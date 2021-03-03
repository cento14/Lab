
import numpy as np
import scipy, scipy.optimize, scipy.ndimage

from astropy.io import fits
from astropy.table import Table 
from astropy.coordinates import SkyCoord 
from astropy import units as u
#from astropy.modeling import models
import scipy.optimize as optimization

import matplotlib.pyplot as plt

import scipy.integrate as integrate

def func_pulsed(x, k1_fit, s1_fit, C=0, m1_fit=0):
            return C + k1_fit * np.exp((-(x-m1_fit)**2)/(2*s1_fit**2))  # 1 gauss

def gaucost(x, k, sig, mu, C):
    return C + k * np.exp((-(x-mu)**2)/(2*sig**2))   # 1 gauss


def PsfGauss(evt, width=5., bin=.01):  
  x=evt['DETX']
  y=evt['DETY']
  #a1, b1 = np.histogram(x,bins=np.arange(-1,1,.05))
  a1, b1, cc = plt.hist(x,bins=np.arange(-width,width,bin))
  b2 = (b1 + np.roll(b1,1))/2.0
  b2 = np.delete(b2, 0)
  err_a1 = a1**0.5
  err_a1 = [1e6 if x==0 else x for x in err_a1]
  #guess = np.array([100, 0.20,0.,0.])
  guess = np.array([10000., 0.1, 0.0, 1000.])
  fitParams, fitCovariances = scipy.optimize.curve_fit(gaucost, b2, a1, guess, sigma=err_a1, maxfev=100000)
  sigma = fitParams, fitCovariances 
  return sigma

def PsfGauss2D(evt, width=5., bin=.01):
    x=evt['DETX']
    y=evt['DETY']
    #a1, b1 = np.histogram(x,bins=np.arange(-1,1,.05))
    th2 = (x**2 + y**2)
    #print max(th2)
    a1, b1, cc = plt.hist(th2,bins=np.arange(-width,width,bin))
    #plt.xlim(0,0.7)
    #plt.show()
    b2 = (b1 + np.roll(b1,1))/2.0
    b2 = np.delete(b2, 0)
    err_a1 = a1**0.5
    err_a1 = [1e6 if x==0 else x for x in err_a1]
    guess = np.array([10000., 0.1, 0.0, 1000.])
    fitParams, fitCovariances = scipy.optimize.curve_fit(gaucost, b2, a1, guess, sigma=err_a1, maxfev=100000)
    sigma = fitParams, fitCovariances 
    return sigma


def PsfRCont(evt, dir_src,p=0.68):  
  dir_ph=SkyCoord(ra=evt['RA']*u.degree,dec=evt['DEC']*u.degree,frame='fk5')
  d = dir_ph.separation(dir_src)
  ds=np.sort(d)
  rcont = ds.deg[int(len(ds)*p)]
  return rcont


def PsfSmooth(psf):
  for i in np.arange(len(psf['SIGMA_1'][0])):
    cn=psf['SIGMA_1'][0][i]
    psf_new=psf
    if sum(cn) != 0 :
      w=np.where(cn != 0)
      cn1=cn[w]
      x=np.arange(len(cn1))
      #p=polyfit(x,cn1,3)
      #c=p[3]+x*p[2]+x**2.*p[1]+x**3.*p[0]
      p=np.polyfit(x,cn1,2)
      c=p[2]+x*p[1]+x**2.*p[0]
      c[c<0] = 0.
      psf[0][5][i][w] = c
  return psf_new


def AeffSmooth(aeff):  
  for i in np.arange(len(aeff['EFFAREA'][0])):
    cn=aeff['EFFAREA'][0][i]
    aeff_new=aeff
    if sum(cn) != 0 :
        w=np.where(cn != 0)
        cn1=cn[w]
        x=np.arange(len(cn1))
        #p=polyfit(x,cn1,3)
        #c=p[3]+x*p[2]+x**2.*p[1]+x**3.*p[0]
        p=np.polyfit(x,cn1,2)
        c=p[2]+x*p[1]+x**2.*p[0]
        c[c<0] = 0.
        aeff[0][4][i][w]= c
  return aeff_new


def EdispSmooth(edisp,smoothing = 2.):
  edisp_new=edisp
  for j_theta in np.arange(len(edisp['THETA_LO'][0])):
    for i_elo in np.arange(len(edisp['ETRUE_LO'][0])):
      guess = np.array([0, 100, 1, 0.1])
      q = edisp['MATRIX'][0][j_theta]
      q_smooth = scipy.ndimage.filters.gaussian_filter(q[:,i_elo],smoothing,order=0)
      emigra=(edisp['MIGRA_LO'][0]+edisp['MIGRA_HI'][0])/2.
      fitParams, fitCovariances = scipy.optimize.curve_fit(func_pulsed, emigra, q_smooth, guess, 0.01, maxfev=100000)
      edisp['MATRIX'][0][j_theta][:,i_elo] = func_pulsed(emigra,fitParams[0],fitParams[1],fitParams[2],fitParams[3])
  
  return edisp_new


def BkgSmooth(bkg,smoothing = 3.):
  bkg_new=bkg
  for i in np.arange(len(bkg_new['BGD'][0])):
    bkg_new[0][6][i] = scipy.ndimage.filters.gaussian_filter(bkg[0][6][i],smoothing,order=0)
  return bkg_new

def BkgSmooth_MC(bkg,smoothing = 3.):               #per leggere il background da MonteCarlo
  bkg_new=bkg
  for i in np.arange(len(bkg_new)):
    bkg_new = scipy.ndimage.filters.gaussian_filter(bkg,smoothing,order=0)
    #plt.show()
  return bkg_new

def IrfUpdate(table,irf_file):
  tempfile='.irftemp'
  table.write(tempfile,format='fits', overwrite=True)
  table2 = fits.open(tempfile)
  irf = fits.open(irf_file, mode='update')
  #irf[table2[1].header['EXTNAME']]=table2[1]
  fits.update(irf_file,table2[1].data,table2[1].header,table2[1].header['EXTNAME'])
  irf.close()
  return table2[1].header['EXTNAME']+' updated in '+irf_file



def mkAeffTable(emin,emax,detmin,detmax,a):
  sp=Table()
  sp['ENERG_LO']= [emin]
  sp['ENERG_LO'].unit='TeV'

  sp['ENERG_HI']= [emax]
  sp['ENERG_HI'].unit='TeV'

  sp['THETA_LO']= [detmin]
  sp['THETA_LO'].unit='deg'

  sp['THETA_HI']= [detmax]
  sp['THETA_HI'].unit='deg'

  sp['EFFAREA']=[[a]*len(detmin)]
  sp['EFFAREA'].unit='m2'

  sp['EFFAREA_RECO']=[[a]*len(detmin)]
  sp['EFFAREA_RECO'].unit='m2'

  sp.meta['EXTNAME']='EFFECTIVE AREA'
  
  return sp


def mkPsfTable(emin,emax,detmin,detmax,sigmas):
  sp=Table()
  sp['ENERG_LO']= [emin]
  sp['ENERG_LO'].unit='TeV'

  sp['ENERG_HI']= [emax]
  sp['ENERG_HI'].unit='TeV'

  sp['THETA_LO']= [detmin]
  sp['THETA_LO'].unit='deg'

  sp['THETA_HI']= [detmax]
  sp['THETA_HI'].unit='deg'

  sp['SCALE']=[[1./(2*np.pi*np.array(sigmas)**2.)]*len(detmin)]

  sp['SIGMA_1']=[[sigmas]*len(detmin)]
  sp['SIGMA_1'].unit='deg'

  sp['AMPL_2']=[[sigmas-sigmas]*len(detmin)]
 
  sp['SIGMA_2']=[[sigmas-sigmas]*len(detmin)]
  sp['SIGMA_2'].unit='deg'
  
  sp['AMPL_3']=[[sigmas-sigmas]*len(detmin)]
  
  sp['SIGMA_3']=[[sigmas-sigmas]*len(detmin)]
  sp['SIGMA_3'].unit='deg'

  sp.meta['EXTNAME']='POINT SPREAD FUNCTION'
  
  return sp

def mkbkgTable(DET_L,DET_H,emin,emax,bkg):
    background=Table()
    #background = Table(bkg_irf)
    background['DETX_LO'] = np.float32([DET_L])
    background['DETX_LO'].unit = 'deg'

    background['DETX_HI'] = np.float32([DET_H])
    background['DETX_HI'].unit = 'deg'

    background['DETY_LO'] = np.float32([DET_L])
    background['DETY_LO'].unit = 'deg'

    background['DETY_HI'] = np.float32([DET_H])
    background['DETY_HI'].unit = 'deg'

    background['ENERG_LO'] = np.float32([emin])
    background['ENERG_LO'].unit = 'TeV'

    background['ENERG_HI'] = np.float32([emax])
    background['ENERG_HI'].unit = 'TeV'

    background['BGD'] = np.float32([bkg])
    background['BGD'].unit = '1/s/MeV/sr'

    background.meta['EXTNAME'] = 'BACKGROUND'

    return background


def mkEdispTable(emin,emax,miglo,mighi,rmf):
  sp=Table()
  sp['ETRUE_LO']= [emin]
  sp['ETRUE_LO'].unit='TeV'

  sp['ETRUE_HI']= [emax]
  sp['ETRUE_HI'].unit='TeV'

  sp['MIGRA_LO']= [miglo]

  sp['MIGRA_HI']= [mighi]

  sp['THETA_LO']= [[0.0]]
  sp['THETA_LO'].unit='deg'

  sp['THETA_HI']= [[5.0]]
  sp['THETA_HI'].unit='deg'

  sp['MATRIX']=[[rmf]]

  sp.meta['EXTNAME']='ENERGY DISPERSION'
  
  return sp





############### Per fit sulle distribuzioni degli eventi

def func_pulsedgau(x, k11, s11): # for 1 gauss
    return k11/np.sqrt(2*np.pi*s11**2)  * np.exp((-(x)**2)/(2*s11**2))  # 1 gauss     ...

def func_pulsed2gau(x, k11, s11,k22, x0): # for 2 gauss
    return k11/np.sqrt(2*np.pi*s11**2)  * np.exp((-(x-x0)**2)/(2*s11**2)) + k22/np.sqrt(2*np.pi*sigma_bkg[i]**2)  * np.exp((-(x)**2)/(2*sigma_bkg[i]**2))   # 3 gauss (1 data, 2 bkg)     ...

def fit_hist_bkg(ene_d, ene_u, evt2, ene_evt2, delta_integr, range, guess=np.array([10000, 0.1]),func=func_pulsedgau):
    print "Estimate the shape of the distribution: ..."
    psf = np.zeros(len(ene_d),float)
    norm_bkg = np.zeros((3,len(ene_d)),float)
    sigma_bkg = np.zeros((3,len(ene_d)),float)
    err_norm_bkg = np.zeros((3,len(ene_d)),float)
    err_sigma_bkg = np.zeros((3,len(ene_d)),float)
    for i in np.arange(0,len(ene_d)):
        bb2, aa1, bb3, nevt = make_hist(ene_d, ene_u, evt2, ene_evt2, delta_integr, i, range)
        #bb2 = bb2[aa1>0]
        #aa1 = aa1[aa1>0]
        err_aa1 = aa1**0.5
        err_aa1 = [1e6 if x==0 else x for x in err_aa1]
        fitParams, fitCovariances = optimization.curve_fit(func, bb2, aa1, guess, sigma=err_aa1, maxfev=100000)
        norm_bkg[0,i] = fitParams[0]
        sigma_bkg[0,i] = fitParams[1]
        err_norm_bkg[0,i] = fitCovariances[0,0]**0.5
        err_sigma_bkg[0,i] = fitCovariances[1,1]**0.5
        #        if fitParams[4]:
        #            norm_bkg[1,i] = fitParams[4]
        #            sigma_bkg[1,i] = fitParams[5]
        #            err_norm_bkg[1,i] = fitCovariances[4,4]**0.5
        #            err_sigma_bkg[1,i] = fitCovariances[5,5]**0.5
        #        if fitParams[7]:
        #            norm_bkg[2,i] = fitParams[7]
        #            sigma_bkg[2,i] = fitParams[8]
        #            err_norm_bkg[2,i] = fitCovariances[7,7]**0.5
        #            err_sigma_bkg[2,i] = fitCovariances[8,8]**0.5
        print "sigma1:",fitParams[1], '+\-', fitCovariances[1,1]**0.5, "; Norm1:",fitParams[0], '+\-', fitCovariances[0,0]**0.5
        #print "mean1:",fitParams[3], '+\-', fitCovariances[3,3]**0.5, "; mean2:",fitParams[6], '+\-', fitCovariances[6,6]**0.5, "; mean3:",fitParams[9], '+\-', fitCovariances[9,9]**0.5
        plt.ylim(0,max(aa1)+10.)
        plt.errorbar(bb2, aa1, yerr = aa1**0.5, fmt = 'red')
        plt.plot(bb2,func(bb2,*fitParams))
        plt.show()
    return norm_bkg, err_norm_bkg, sigma_bkg, err_sigma_bkg


def fit_hist_theta2(ene_d, ene_u, evt2, ene_evt2, delta_integr, range, guess=np.array([10000, 0.1]),func=func_pulsedgau):
    print "Estimate the shape of the distribution: ..."
    psf = np.zeros(len(ene_d),float)
    norm_bkg = np.zeros((3,len(ene_d)),float)
    sigma_bkg = np.zeros((3,len(ene_d)),float)
    err_norm_bkg = np.zeros((3,len(ene_d)),float)
    err_sigma_bkg = np.zeros((3,len(ene_d)),float)
    for i in np.arange(0,len(ene_d)):
        bb2, aa1, bb3, nevt = make_hist_theta2(ene_d, ene_u, evt2, ene_evt2, delta_integr, i, range)
        #bb2 = bb2[aa1>0]
        #aa1 = aa1[aa1>0]
        err_aa1 = aa1**0.5
        err_aa1 = [1e6 if x==0 else x for x in err_aa1]
        fitParams, fitCovariances = optimization.curve_fit(func, bb2, aa1, guess, sigma=err_aa1, maxfev=100000)
        norm_bkg[0,i] = fitParams[0]
        sigma_bkg[0,i] = fitParams[1]
        err_norm_bkg[0,i] = fitCovariances[0,0]**0.5
        err_sigma_bkg[0,i] = fitCovariances[1,1]**0.5
#        if fitParams[2]:
#                    norm_bkg[1,i] = fitParams[2]
#                    sigma_bkg[1,i] = fitParams[3]
#                    err_norm_bkg[1,i] = fitCovariances[2,2]**0.5
#                    err_sigma_bkg[1,i] = fitCovariances[3,3]**0.5
        #        if fitParams[7]:
        #            norm_bkg[2,i] = fitParams[7]
        #            sigma_bkg[2,i] = fitParams[8]
        #            err_norm_bkg[2,i] = fitCovariances[7,7]**0.5
        #            err_sigma_bkg[2,i] = fitCovariances[8,8]**0.5
        print "sigma1:",fitParams[1], '+\-', fitCovariances[1,1]**0.5, "; Norm1:",fitParams[0], '+\-', fitCovariances[0,0]**0.5, "Cost = ", fitParams[2]
        #print "mean1:",fitParams[3], '+\-', fitCovariances[3,3]**0.5, "; mean2:",fitParams[6], '+\-', fitCovariances[6,6]**0.5, "; mean3:",fitParams[9], '+\-', fitCovariances[9,9]**0.5
        plt.ylim(0,max(aa1)+10.)
        plt.errorbar(bb2, aa1, yerr = aa1**0.5, fmt = 'red')
        plt.plot(bb2,func(bb2,*fitParams))
        plt.show()
    return norm_bkg, err_norm_bkg, sigma_bkg, err_sigma_bkg


def fit_hist_src_bkg(ene_d, func_pulsed2gau, evt1, ene_evt1, delta_integr, offset):
    print "\nEstimate the shape of the distribution: ..."
    norm_src = np.zeros((len(ene_d)),float)
    norm_bkg2 = np.zeros((len(ene_d)),float)
    sigma_src = np.zeros((len(ene_d)),float)
    global sigma_bkg, i
    for i in np.arange(0,len(ene_d)):
        bb2, aa1, bb3 = make_hist(ene_d, evt1, ene_evt1, delta_integr, i)
        #bb2 = bb2[aa1>0]
        #aa1 = aa1[aa1>0]
        err_aa1 = aa1**0.5
        err_aa1 = [1e2 if x==0 else x for x in err_aa1]
        guess = np.array([10000, 0.10,100, offset])
        fitParams, fitCovariances = optimization.curve_fit(func_pulsed2gau, bb2, aa1, guess, sigma=err_aa1, maxfev=100000)
        print "sigma :",fitParams[1], '+\-', fitCovariances[1,1]**0.5
        plt.ylim(0,max(aa1)+10.)
        plt.errorbar(bb2, aa1, yerr = aa1**0.5, fmt = 'red')
        plt.plot(bb2,func_pulsed2gau(bb2,fitParams[0],fitParams[1],fitParams[2],fitParams[3]))
        plt.show()
        norm_src[i] = fitParams[0]
        sigma_src[i] = (fitParams[1]**2.)**0.5
        norm_bkg2[i] = fitParams[2]
    return norm_src, sigma_src, norm_bkg2 ,bb3

def make_hist(ene_d, ene_u, evt2, ene_evt2, delta_integr, i, range):
    index2 = np.where((ene_evt2>=ene_d[i]) & (ene_evt2<ene_u[i]))
    x2=evt2['DETX'].data
    x2 = x2[index2]
    y2=evt2['DETY'].data
    y2 = y2[index2]
    aa = plt.hist(x2,bins=np.arange(-1.*range,range,delta_integr))
    aa1=aa[0]
    bb1=aa[1]
    bb2 = (bb1 + np.roll(bb1,1))/2.0
    bb2 = np.delete(bb2, 0)
    bb3 = bb2
    nevt = len(index2[0])
    print nevt
    #print "Number of events:", len(index2)
    return bb2, aa1, bb3, nevt

def make_hist_2d(ene_d, ene_u, evt2, ene_evt2, delta_integr, i, range):
    index2 = np.where((ene_evt2>=ene_d[i]) & (ene_evt2<ene_u[i]))
    x2=evt2['DETX'].data
    x2 = x2[index2]
    y2=evt2['DETY'].data
    y2 = y2[index2]
    aa = plt.hist2d(x2,y2,bins=np.arange(-1.*range,range+delta_integr,delta_integr))         #conteggi
    #plt.show()
    aa1=aa[0]
    aaE1 = aa[0]/(ene_u[i]-ene_d[i])      #normalizzato all'intervallo energetico
    return aa1, aaE1


def make_hist_theta2(ene_d, ene_u, evt2, ene_evt2, delta_integr, i, range):
    index2 = np.where((ene_evt2>=ene_d[i]) & (ene_evt2<ene_u[i]))
    x2=evt2['DETX'].data
    x2 = x2[index2]
    y2=evt2['DETY'].data
    y2 = y2[index2]
    theta2 = (x2**2 + y2**2)
    aa = plt.hist(theta2,bins=np.arange(-1.*range,range,delta_integr))
    aa1=aa[0]
    bb1=aa[1]
    bb2 = (bb1 + np.roll(bb1,1))/2.0
    bb2 = np.delete(bb2, 0)
    bb3 = bb2
    nevt = len(index2[0])
    print nevt
    #print "Number of events:", len(index2)
    return bb2, aa1, bb3, nevt


def calc_effective_area_irf(energy_mean, delta_energy_mean, aeff, tab, ene_d, bb3, norm_src, sigma_src, delta_integr, exposure):
    crab_model = astropy.modeling.powerlaws.LogParabola1D(amplitude=3.23e-11, x_0=1, alpha=2.47, beta=0.24)
    evt_src = np.zeros((len(energy_mean)),float)
    Aeff = np.zeros((len(energy_mean)),float)
    for j in np.arange(len((aeff['THETA_LO'][0]))):
        for i in np.arange(len(ene_d)):
            evt_src[i] = np.sum(func_pulsedgau(bb3[(bb3>-1.*sigma_src[i]) & (bb3<1.*sigma_src[i])],norm_src[i],sigma_src[i])*delta_integr)
            Aeff[i] = evt_src[i] / exposure / crab_model(energy_mean[i]) * 1e-4 / delta_energy_mean[i]
            aeff[0][4][j][i] = Aeff[i]
            tab[0][5][j][i]= sigma_src[i]
    return aeff, tab

def calc_bkg_irf(bkg_irf, emin, DET_L, norm_bkg, sigma_bkg, delta_energy_mean, exposure, ddeg):
    #m=np.zeros((len(emin), len(DET_L), len(DET_H)),float)
    #bkg_irf deve essere definita in modo da essere una tabella ottenuta da astropy.fits
    idx = 0
    while idx < len(emin):
        if len(DET_L) % 2 == 0:
            n = np.fromfunction(lambda i, j: norm_bkg[0][idx]/(2*np.pi*sigma_bkg[0][idx]**2)
                                * np.exp((-(((i+0.5-len(DET_L)/2.)*ddeg)**2+((j+0.5-len(DET_L)/2.)*ddeg)**2))/(2*sigma_bkg[0][idx]**2))
                                #/ ((np.pi*(( (bkg_irf['DETX_HI'][0][i] * np.pi/180.)**2 - (bkg_irf['DETX_LO'][0][i] * np.pi/180.)**2)))**2)**0.5
                                #/ (np.pi * (np.pi/180.)**2 * (( ((i+0.5-len(DET_L)/2.)*0.2)**2 - ((j+0.5+0.5-len(DET_L)/2.)*0.2)**2 )**2)**0.5)
                                / (ddeg*ddeg * (np.pi/180.)**2)
                                / (delta_energy_mean[idx]*1e6) / exposure , (50, 50), dtype=int)
        
        else:
            n = np.fromfunction(lambda i, j: norm_bkg[0][idx]/(2*np.pi*sigma_bkg[0][idx]**2)
                                * np.exp((-(((i-len(DET_L)/2.)*ddeg)**2+((j-len(DET_L)/2.)*ddeg)**2))/(2*sigma_bkg[0][idx]**2))
                                #/ ((np.pi*(( (bkg_irf['DETX_HI'][0][i] * np.pi/180.)**2 - (bkg_irf['DETX_LO'][0][i] * np.pi/180.)**2)))**2)**0.5
                                #/ (np.pi * (np.pi/180.)**2 * (( ((i+0.5-len(DET_L)/2.)*0.2)**2 - ((j+0.5+0.5-len(DET_L)/2.)*0.2)**2 )**2)**0.5)
                                / (ddeg*ddeg * (np.pi/180.)**2)
                                / (delta_energy_mean[idx]*1e6) / exposure , (50, 50), dtype=int)
        bkg_irf[0][6][idx] = n
        idx += 1
    return bkg_irf


#def calc_bkg_irf(bkg_irf, ene_d, norm_bkg, sigma_bkg, delta_energy_mean, exposure):
#    m,n,i = 0,0,0
#    c = len(bkg_irf['DETX_LO'][0])/ 2
#    while i < len(ene_d):
#        while m < len(bkg_irf['DETX_LO'][0]):
#            while n < len(bkg_irf['DETX_LO'][0]):
#                n, m = int(n), int(m)
#                n1 = n + 0.5
#                m1 = m + 0.5
#                #TO DO: correggere le unita' di misura e verificare la procedura
#                bkg_new = norm_bkg[i] * 1. / np.sqrt(2*np.pi*sigma_bkg[i]**2.) * np.exp(- (((c-m1)**2 + (c-n1)**2)* ((bkg_irf['DETX_HI'][0] - bkg_irf['DETX_LO'][0]))[0]) / (2 * sigma_bkg[i]**2.))
#                #                if m == 0 or m==11:
#                #                    print bkg_new, n, m, norm_bkg[i], n1, m1, c, ((c-m1)**2 + (c-n1)**2), ((bkg_irf['DETX_HI'][0] - bkg_irf['DETX_LO'][0]))[0]
#                #ster = (np.pi*(( (bkg_irf['DETX_HI'][0] * np.pi/180)**2 - (bkg_irf['DETX_LO'][0] * np.pi/180)**2)))**2.
#                ster = (np.pi*(( (bkg_irf['DETX_HI'][0] * np.pi/180.)**2 - (bkg_irf['DETX_LO'][0] * np.pi/180.)**2)))**2.
#                bkg_irf[0][6][i][m][n] = bkg_new / (delta_energy_mean[i]*1e6) / ster[m]**0.5/ exposure
#                n += 1
#            n=0
#            m+=1
#        m = 0
#        i += 1
#    return bkg_irf


def fill_aeff_irf(aeff,aeff_refill,aeff_reco_refill):
    for j in np.arange(0,len((aeff['THETA_LO'][0]))):
        aeff[0][4][j][:] = aeff_refill #aeff[0][4][0][:]
        aeff[0][5][j][:] = aeff_reco_refill #aeff[0][5][0][:]
    return aeff

def fill_psf_irf(tab,psf1_refill,psf2_refill,psf3_refill,psf4_refill):
    #for j in np.arange(0,len((tab['THETA_LO'][0]))):                #per riempire anche i bin OFF-AXIS
    for j in range(0,1):             # nel caso vogliamo riempire solo i bin ON-AXIS
        tab[0][4][j][:]= (1./(2.*np.pi*psf1_refill**2.)) #tab[0][4][0][:] #aggiorna la colonna scale
        tab[0][5][j][:]= (psf1_refill**2.)**0.5 #tab[0][5][0][:] #aggiorna la prima sigma
        tab[0][7][j][:]= (psf2_refill**2.)**0.5 #tab[0][7][0][:] #aggiorna la seconda sigma
        tab[0][9][j][:]= (psf3_refill**2.)**0.5 #tab[0][9][0][:] #aggiorna la terza sigma
        tab[0][6][j][:]= (psf4_refill**2.)**0.5 #tab[0][6][0][:] #aggiorna la second norm
    return tab

def fit_hist_MC(ene_d, ene_u, evt2, ene_evt2, delta_integr, range, guess=np.array([10000, 0.1]),func=func_pulsedgau):
    print "Estimate the shape of the distribution: ..."
    norm_bkg = np.zeros((3,len(ene_d)),float)
    sigma_bkg = np.zeros((3,len(ene_d)),float)
    err_norm_bkg = np.zeros((3,len(ene_d)),float)
    err_sigma_bkg = np.zeros((3,len(ene_d)),float)
    aeff_corr = np.zeros((3,len(ene_d)),float)
    for i in np.arange(0,len(ene_d)):
        bb2, aa1, bb3, nevt = make_hist(ene_d, ene_u, evt2, ene_evt2, delta_integr, i, range)
        #bb2 = bb2[aa1>0]
        #aa1 = aa1[aa1>0]
        err_aa1 = aa1**0.5
        err_aa1 = [1e6 if x==0 else x for x in err_aa1]
        fitParams, fitCovariances = optimization.curve_fit(func, bb2, aa1, guess, sigma=err_aa1, maxfev=100000)
        norm_bkg[0,i] = fitParams[0]
        sigma_bkg[0,i] = fitParams[1]
        err_norm_bkg[0,i] = fitCovariances[0,0]**0.5
        err_sigma_bkg[0,i] = fitCovariances[1,1]**0.5
        #aeff_corr[0,i] = (np.sqrt(2*np.pi)*norm_bkg[0,i]*sigma_bkg[0,i]) / nevt
        #print "Number of events under the gaussian", aeff_corr[0,i]
        #print "sigma1:",np.round(fitParams[1],3), '+\-', np.round(fitCovariances[1,1]**0.5,5), "; Norm1:",np.round(fitParams[0],2), '+\-', np.round(fitCovariances[0,0]**0.5,2), "\n"
        if len(fitParams)==5:
            norm_bkg[1,i] = fitParams[3]
            sigma_bkg[1,i] = fitParams[4]
            err_norm_bkg[1,i] = fitCovariances[3,3]**0.5
            err_sigma_bkg[1,i] = fitCovariances[4,4]**0.5
            aeff_corr[0,i] = ((np.sqrt(2*np.pi)*norm_bkg[0,i]*np.sqrt((sigma_bkg[0,i])**2))/delta_integr + (np.sqrt(2*np.pi)*norm_bkg[1,i]*np.sqrt((sigma_bkg[1,i])**2))/delta_integr) / nevt
            print "Number of events under the gaussian", aeff_corr[0,i], nevt, ((np.sqrt(2*np.pi)*norm_bkg[0,i]*np.sqrt((sigma_bkg[0,i])**2)) + (np.sqrt(2*np.pi)*norm_bkg[1,i]*np.sqrt((sigma_bkg[1,i])**2)))/delta_integr
            print "sigma1:",np.round(fitParams[1],3), '+\-', np.round(fitCovariances[1,1]**0.5,5), "; Norm1:",np.round(fitParams[0],2), '+\-', np.round(fitCovariances[0,0]**0.5,2), "\nsigma2:",np.round(fitParams[4],3), '+\-', np.round(fitCovariances[4,4]**0.5,4), "; Norm2:",np.round(fitParams[3],2), '+\-', np.round(fitCovariances[3,3]**0.5,2)
        if len(fitParams)==7:
            norm_bkg[2,i] = fitParams[5]
            sigma_bkg[2,i] = fitParams[6]
            err_norm_bkg[2,i] = fitCovariances[5,5]**0.5
            err_sigma_bkg[2,i] = fitCovariances[6,6]**0.5
            print "sigma1:",np.round(fitParams[1],3), '+\-', np.round(fitCovariances[1,1]**0.5,5), "; Norm1:",np.round(fitParams[0],2), '+\-', np.round(fitCovariances[0,0]**0.5,2), "\nsigma2:",np.round(fitParams[4],3), '+\-', np.round(fitCovariances[4,4]**0.5,4), "; Norm2:",np.round(fitParams[3],2), '+\-', np.round(fitCovariances[3,3]**0.5,2), "\nsigma3:", np.round(fitParams[6],3), '+\-', np.round(fitCovariances[6,6]**0.5,4), "; Norm3:",np.round(fitParams[5],2), '+\-', np.round(fitCovariances[5,5]**0.5,2)
        plt.ylim(0,max(aa1)+10.)
        plt.errorbar(bb2, aa1, yerr = aa1**0.5, fmt = 'red')
        plt.plot(bb2,func(bb2,*fitParams))
        plt.show()
    return norm_bkg, err_norm_bkg, sigma_bkg, err_sigma_bkg, aeff_corr


