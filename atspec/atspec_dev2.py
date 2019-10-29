#/pollux/giuliani/anaconda2/bin/python

import sys,os
from astropy import units as u
from astropy.coordinates import SkyCoord as coords
from astropy.io import fits, ascii  
from astropy.table import Table 
import numpy
import matplotlib.pyplot as plt 

class pointsource:
 def __init__(self,evtfile="",insfile="",l=0.,b=-100., name=''):
  

  if evtfile == '' :
    evtfile = raw_input('Enter events file : ')
    
  if insfile == '' :
    insfile = raw_input('Enter response file : ')
    
  if b == -100 :
      l =  raw_input('Enter G. longitude : ')
      b =  raw_input('Enter G. latitude : ')
      
  if name != '' :
    dir_src=coords.from_name(name)
    l=dir_src.galactic.l.value
    b=dir_src.galactic.b.value

  self.l_src=l
  self.b_src=b

  ### Events

  file_p = fits.open(evtfile)
  evt=file_p[1].data

  self.tsec=file_p[1].header['LIVETIME']
  #self.tsec=t*3600.
  
  ra_pnt=file_p[1].header['RA_PNT']
  dec_pnt=file_p[1].header['DEC_PNT']
 
  self.enetot=evt.field("ENERGY")
  ra=evt.field("RA")
  dec=evt.field("DEC")

  dir_ph=coords(ra=ra*u.degree,dec=dec*u.degree,frame='fk5')
    
  
  #if b_onaxis == -100 :
    #b_onaxis=b
    #l_onaxis=l

  ### Ext. direction

  dir_src=coords(l=self.l_src*u.degree,b=self.b_src*u.degree,frame='galactic')
  dir_onaxis=coords(ra=ra_pnt*u.degree,dec=dec_pnt*u.degree,frame='fk5')
  
  self.d = dir_ph.separation(dir_src)
  oa=dir_src.separation(dir_onaxis) 

  print 'Exposure time : ', self.tsec/3600., ' hrs'
  print 'Extraction region :' , dir_src.galactic
  
  
  ### IRF

  check=open(insfile)
  inizio=check.readline()
  check.close()
  
  #if inizio[0:6] == 'SIMPLE' :
  try :
      
    irf=fits.open(insfile)

    #irf=fits.open(insfile)
  
    aeff2d=irf[1].data
    lene_c2= (numpy.log10(aeff2d.field('ENERG_LO')[0]) + numpy.log10(aeff2d.field('ENERG_HI')[0]) ) *.5
    aeff_c = aeff2d.field('EFFAREA')[0][0]
  
    psf2d=irf[2].data
    lene_c=  (numpy.log10(psf2d.field('ENERG_LO')[0]) + numpy.log10(psf2d.field('ENERG_LO')[0]) ) *.5
    r80_c=psf2d.field('SIGMA_1')[0][0]*1.8  
    r68_c=psf2d.field('SIGMA_1')[0][0]*1.51
    
  #else:
  except:
      
    fov=3.2 *u.degree
    irf = open(insfile, 'r')

    lene_c=[0]
    aeff_c=[0]
    r68_c=[0]
    r80_c=[0]

    for li in irf: 
      ll=li.split()
      if len(ll) == 7 : 
      #print ll[3]
        lene_c=lene_c+[float(ll[0])]
        aeff_c=aeff_c+[float(ll[1])]
        r68_c=r68_c+[float(ll[2])]
        r80_c=r80_c+[float(ll[3])]
    
    lene_c2=lene_c
    irf.close()


  self.lene_c=numpy.array(lene_c[1:])
  self.lene_c2=numpy.array(lene_c2[1:])
  self.aeff_c=numpy.array(aeff_c[1:]) #*numpy.exp(-((oa/fov)**4.)/2.)
  self.r68_c=numpy.array(r68_c[1:])
  self.r80_c=numpy.array(r80_c[1:])

  #print oa, numpy.exp(-((oa/fov)**4.)/2.)
  
  
  ### S.Extraction


  self.exrad=numpy.interp(self.enetot,10.**self.lene_c,self.r80_c)*u.degree
  self.corr=0.8
  #self.corr=0.75
  #self.corr=0.68
  

  ww=numpy.where(self.d < .5*u.degree) 
  plt.hist2d(dir_ph[ww].galactic.l,dir_ph[ww].galactic.b,bins=200)
  
  acic=numpy.arange(0,2*numpy.pi,.01)
  radius=numpy.mean(self.exrad.value)
  #print radius,l,b
  plt.plot(radius*numpy.cos(acic)+l, radius*numpy.sin(acic)+b,'w-')
  plt.plot(2*radius*numpy.cos(acic)+l, 2*radius*numpy.sin(acic)+b,'w--')
  plt.plot(3*radius*numpy.cos(acic)+l, 3*radius*numpy.sin(acic)+b,'w--')
  
  plt.show()

 
 def getspec(self, emin=1., emax=100., nch=10, outfile='', sigmas=3., verbose=False):

    lemin=numpy.log10(emin)
    lemax=numpy.log10(emax)

    self.sb = ( lemax - lemin) / nch
    self.canali_min=  numpy.array(range(int(nch))) *self.sb  +lemin

    fflux=[0]
    eerr=[0]
    
    if verbose :
      print 'Ener. ', 'Non, ' , 'Noff, ' , 'Sigmas, ' ,  'Detection?'

    for lene_sig in self.canali_min:
      
      ene_sig=10.0**lene_sig

      won =  numpy.where(  (self.d < self.exrad) * (self.enetot >= ene_sig) * (self.enetot < ene_sig*10.**self.sb) )
      woff = numpy.where(  (self.d > 3*self.exrad) * (self.d < 4*self.exrad) * (self.enetot > ene_sig) * (self.enetot < ene_sig*10.**self.sb) )
   
      alpha=1./7.

      n_on=float(numpy.size(won))
      n_off=float(numpy.size(woff))
    
      signal= ( n_on-n_off*alpha  ) 
      
      if (n_off != 0)*( n_on != 0.) : 
        sigm=numpy.sqrt(2. *( n_on *numpy.log( (1+alpha)/alpha * (n_on / (n_on +n_off))  ) + n_off *numpy.log( (1+alpha) *( n_off /(n_on +n_off) ) ) ) )
      elif n_off == 0. :
        sigm=numpy.sqrt(n_on)
      else: 
        sigm=0.
    

      aeff=numpy.interp(lene_sig+self.sb/2.,self.lene_c2,self.aeff_c)*1e4   

  #print aeff

      if aeff > 0 : 
        flux=signal/aeff/self.tsec  /(ene_sig*(10.**self.sb - 1.)) /self.corr
        err= numpy.maximum((numpy.sqrt(n_on)+numpy.sqrt(n_off)*alpha)  ,  1.) /aeff/self.tsec/(ene_sig*(10.**self.sb - 1.))/self.corr
      else:
        flux=0.0
        err=0.0

      nosig=(sigm < sigmas)

      if nosig : 
        err=flux+err
        flux=0 

      #ene_ch=10.**(lene_sig+self.sb/2.)
      fflux=fflux+[flux]
      eerr=eerr+[err]
      
      if verbose :
        print numpy.round(10.**(lene_sig+self.sb/2.),3), ' ', long(n_on),' ', numpy.round(n_off*alpha,2),' ', numpy.round(sigm,2) ,' ', sigm > sigmas

      #print   signal, sigm
      #print
  
  
  
    ene_ch=10.**(self.canali_min+self.sb/2.)
    fflux=fflux[1:]              
    eerr=eerr[1:]
    
    
    ### Output

    out=Table([ene_ch*1e12, fflux*ene_ch**2.*1.602, eerr*ene_ch**2.*1.602, eerr*ene_ch**2.*1.602],names=['Energy','Flux','Flux_Error_Min','Flux_Error_Max'])
      
    #out=Table
    #out['ENERGY']=ene_ch*1e12
    #out['FLUX']=fflux*ene_ch**2.*1.6
    #out['FLUX_ERROR_MIN']=eerr*ene_ch**2.*1.6
    #out['FLUX_ERROR_MAX']=eerr*ene_ch**2.*1.6

    out.meta['EXTNAME']='DATA'
    #modello.meta['']=''

    if verbose :
      out.pprint()
    
    if outfile != '' : 
      out.write(outfile, format='fits', overwrite='True')

 
    return fflux, eerr, ene_ch
  
    
def script():
  
  print 'pippo'
  s=pointsource()



#script()




