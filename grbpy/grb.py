#!/usr/bin/python

###########################################################
# Authors : S. Cutini (ASDC) and A. Giuliani (IASF Milano)
# Last-Modified: 12/02/2018
###########################################################


import string
import os

#import pyfits
import astropy.io.fits as pyfits
from astropy import units as u
from astropy.coordinates import SkyCoord as sky

from numpy import * 
import sys
#from math import *
#import coords
#import pydl


####### FIND #######


def find(GRB_evt,GRB_log,paramfile):

    #GRB_evt = sys.argv[1]
    #GRB_log = sys.argv[2]

    parfile=open(paramfile,"r")



    # Sempre utili 

    dreq =  3.141592/180.0
    dreq1 = 1./(3.141592/180.0)

    # Leggo par file

    ciccia=parfile.readline()
    GRB_time = float(parfile.readline()) 
    GRB_ra = float(parfile.readline()) 
    GRB_dec = float(parfile.readline()) 
    raggio=float(parfile.readline())
    t1s=float(parfile.readline())
    t2s=float(parfile.readline())
    t1b=float(parfile.readline())
    t2b=float(parfile.readline())

    fov=float(parfile.readline())
    #ea_th=80.
    ea_th=float(parfile.readline())

    print ''
    print 'GRB T0 :',GRB_time
    print 'GRB (Ra, Dec) :',GRB_ra, GRB_dec
    print 'Ricerca eventi (Raggio):',raggio
    print 'Ricerca eventi (Tmin, Tmax):',t1s, t2s
    print 'Background (Tmin, Tmax):',t1b, t2b
    print 'F.O.V.  :',fov
    print 'Albedo cut  :',ea_th
    print ''

    # Leggo log file

    hdulist_log = pyfits.open(GRB_log)
    tbdata_log = hdulist_log[1].data
    TIME_log = tbdata_log.field('TIME')
    RA_earth= tbdata_log.field('EARTH_RA')
    DEC_earth= tbdata_log.field('EARTH_DEC')
    livetime= tbdata_log.field('LIVETIME')
    ra_punt=tbdata_log.field('ATTITUDE_RA_Y')
    dec_punt=tbdata_log.field('ATTITUDE_DEC_Y')
    phase_log=tbdata_log.field('PHASE')

    TIMEnew_log=TIME_log-GRB_time


    # Leggo evt file

    hdulist = pyfits.open(GRB_evt)
    tbdata = hdulist[1].data
    RAcol = tbdata.field('RA')
    DECcol = tbdata.field('DEC')   
    TIMEcol = tbdata.field('TIME') 
    PH_col = tbdata.field('PH_EARTH')
    EVcol = tbdata.field('EVSTATUS')
    Enecol = tbdata.field('ENERGY')
    THETAcol = tbdata.field('THETA')
    PHASEcol = tbdata.field('PHASE')
    
    TIMEnew=TIMEcol-GRB_time
    tempi=TIMEnew-TIMEnew+8888.

    dec = DECcol*dreq
    ra = RAcol*dreq
    dec_ea = DEC_earth*dreq
    ra_ea = RA_earth*dreq

    #Vx=zeros(len(dec),float)
    #Vy=zeros(len(dec),float)
    #Vz=zeros(len(dec),float)

    #argx=zeros(len(dec),float)
    #argy=zeros(len(dec),float)
    #argz=zeros(len(dec),float)
    #arg=zeros(len(dec),float)

    #DELTA=zeros(len(dec),float)

    Vx_ea=zeros(len(dec_ea),float)
    Vy_ea=zeros(len(dec_ea),float)
    Vz_ea=zeros(len(dec_ea),float)

    argx_ea=zeros(len(dec_ea),float)
    argy_ea=zeros(len(dec_ea),float)
    argz_ea=zeros(len(dec_ea),float)
    arg_ea=zeros(len(dec_ea),float)

    DELTA_ea=zeros(len(dec_ea),float)

    expo=zeros(len(dec_ea),float)
    expo1=zeros(len(dec_ea),float)
    expo2=zeros(len(dec_ea),float)
    exposure=zeros(len(dec_ea),float)

    #TIMEcol=zeros(len(dec),float)
    #TIMEnew=zeros(len(dec),float)
    #TIMEnew_log=zeros(len(dec_ea),float)
    
    
    t1b=max(min(TIMEnew),t1b)
    t2b=min(max(TIMEnew),t2b)
    

    t1b=max(min(TIME_log-GRB_time),t1b)
    t2b=min(max(TIME_log-GRB_time),t2b)
    #print min(TIME_log-GRB_time),t1b
    print t1b,t2b


    # GRB in coord. cartesiane 

    Vxgrb = cos(GRB_dec*dreq)*cos(GRB_ra*dreq)
    Vygrb = cos(GRB_dec*dreq)*sin(GRB_ra*dreq)
    Vzgrb = sin(GRB_dec*dreq)

    # Distanza angolare grb-eventi

    Vx = cos(dec)*cos(ra)
    Vy = cos(dec)*sin(ra)
    Vz = sin(dec)

    argx = Vxgrb*Vx
    argy = Vygrb*Vy
    argz = Vzgrb*Vz
    arg= argx+argy+argz


    DELTA = arccos(arg)*dreq1


    # Off-axis angle del GRB

    Vx_punt = cos(dec_punt*dreq)*cos(ra_punt*dreq)
    Vy_punt = cos(dec_punt*dreq)*sin(ra_punt*dreq)
    Vz_punt = sin(dec_punt*dreq)

    argx_punt = Vxgrb*Vx_punt
    argy_punt = Vygrb*Vy_punt
    argz_punt = Vzgrb*Vz_punt
    arg_punt= argx_punt+argy_punt+argz_punt

    OFF = arccos(arg_punt)*dreq1


    # Calcolo segnale e background


    #ws=where((Enecol > 0)*(PH_col>90)*(DELTA<raggio)*(TIMEnew>t1s)*(TIMEnew<t2s))
    #print size(ws)
    #wb=where((Enecol > 0)*(PH_col>90)*(DELTA<raggio)*(((TIMEnew>t1b)*(TIMEnew<t1s))+((TIMEnew>t2s)*(TIMEnew<t2b)))  )
    #print size(wb)

    source=0
    bkg=0


    for k in range(len(dec)):    

        
	if DELTA[k] <raggio:  
	  if Enecol[k] > 0.:
            if PH_col[k] >ea_th:
              if PHASEcol[k] != 1:
                  if THETAcol[k]<fov: 
                    tempi[k]=TIMEnew[k] 
                    if TIMEnew[k] > t1s:
                        if TIMEnew[k] < t2s: 
                            source=source+1
                            print '      trovato (t-T0) : ', TIMEnew[k]
                    if (TIMEnew[k] > t1b) and (TIMEnew[k] <t1s) or (TIMEnew[k] > t2s) and (TIMEnew[k] <t2b): 

                            bkg=bkg+1

    print ''
    print "Source :", source
    print "Bkg :", bkg
    print ''

    # calcolo delle significativita con il metodo di Li&Ma
    summ_src=0
    ntot_src=0
    summ_bkg=0
    ntot_bkg=0


    for i in range(len(dec_ea)):

        Vx_ea[i] = cos(dec_ea[i])*cos(ra_ea[i])
        Vy_ea[i] = cos(dec_ea[i])*sin(ra_ea[i])
        Vz_ea[i] = sin(dec_ea[i])
        #Vxgrb_ea = cos(GRB_dec*dreq)*cos(GRB_ra*dreq)
        #Vygrb_ea = cos(GRB_dec*dreq)*sin(GRB_ra*dreq)
        #Vzgrb_ea = sin(GRB_dec*dreq)
        argx_ea[i] = Vxgrb*Vx_ea[i]
        argy_ea[i] = Vygrb*Vy_ea[i]
        argz_ea[i] = Vzgrb*Vz_ea[i]
        arg_ea[i] = argx_ea[i]+argy_ea[i]+argz_ea[i]


        DELTA_ea[i] = arccos(arg_ea[i])*dreq1
        if (DELTA_ea[i]-ea_th+raggio)/(2.*raggio) < 0.:
            expo[i]=0.
        elif (DELTA_ea[i]-ea_th+raggio)/(2.*raggio) > 1.:
            expo[i]=1.

        else:
            expo[i]=(DELTA_ea[i]-ea_th+raggio)/(2*raggio)

        if (fov-OFF[i]+raggio)/(2.*raggio) < 0.:
            expo2[i]=0.
        elif (fov-OFF[i]+raggio)/(2.*raggio) > 1.:
            expo2[i]=1.
        else:
            expo2[i]=(fov-OFF[i]+raggio)/(2*raggio) 

        if (((ra_punt[i] < 0) == 0)*((ra_punt[i]>0) == 0) ) :
            expo2[i]=0        

        if (livetime[i] != 0) and (phase_log[i] != 1):
            expo1[i]=livetime[i]/100.
        else:
            expo1[i]=0.

        exposure[i] = expo[i]*expo1[i]*expo2[i]

        if (TIMEnew_log[i] >t1s) and  (TIMEnew_log[i] < t2s):
            summ_src=exposure[i]+summ_src
            ntot_src=ntot_src+1  

        if (TIMEnew_log[i] > t1b) and (TIMEnew_log[i] <t1s) or (TIMEnew_log[i] > t2s) and (TIMEnew_log[i] <t2b):
            summ_bkg=exposure[i]+summ_bkg
            ntot_bkg=ntot_bkg+1 

    print "mean src occulted ", summ_src/ntot_src
    print "mean bkg occulted", summ_bkg/ntot_bkg
    mean_src=summ_src/ntot_src
    mean_bkg=summ_bkg/ntot_bkg
    
    if (t1s >= t2b) or (t2s <= t1b):
      tback=t2b-t1b
    else:
      if (t1s >= t1b) and (t2s <= t2b):
	tback=t2b-t1b - (t2s-t1s)
      else: 
        print "   !!!! Scegli un altro intervallo di background !!!!"
        tback=0
    
    print tback
    
    #alp = ((t2s-t1s)*mean_src)/((t2b-t1b-(t2s-t1s))*mean_bkg)
    alp = ((t2s-t1s)*mean_src)/(tback*mean_bkg)
    alp1=alp/(1+alp)
    alp2=alp+1
    print "source", source
    print "bkg", bkg/((tback)*mean_bkg)*((t2s-t1s)*mean_src), "(", bkg, ")"
    source1=float(source)
    bkg1=float(bkg)

    if source>0 :
       L1 = math.pow(((source1+bkg1)/source1)*alp1,source1)
       L2 = math.pow(((bkg1+source1)/bkg1)/alp2,bkg1)
       L=L1*L2
       #print "L", alp2
       S=math.sqrt(-2.*math.log(L))
       print "Li&Ma sigma", S


    hdulist.close()
    
    
    #tempi=tempi[where(abs(tempi) < 6000. )]
    #pydl.pist(tempi,5)
    #zoom=float(raw_input('Zoommiamo un po ?? '))
    #tempi=tempi[where(abs(tempi) < zoom )]
    #pydl.pist(tempi,5)
    #raw_input('Please press return to finish...\n')
    
    
    #pyfits.writeto('exp.fits',exposure)
    #pydl.plotta(TIMEnew_log,exposure)



########### GET ######

def get(filename):

   #filename=sys.argv[1]
   fileletto=open(filename,"r")

   riga1=fileletto.readline()  # UTC
   riga2=fileletto.readline()  # OBT
   riga3=fileletto.readline()  # Long
   riga4=fileletto.readline()  # Lat
   riga5=fileletto.readline()  # Theta
   riga6=fileletto.readline()  # Phi
   riga7=fileletto.readline()  # Orbit
   riga8=fileletto.readline()  # Phase

   t0_obt=double(riga2[16:])
   long=double(riga3[16:])
   lat=double(riga4[16:])
   theta=double(riga5[19:])
   phi=double(riga6[19:])
   orbit=int(riga7[19:])

   #t0_obt=double(riga2)
   #long=double(riga3)
   #lat=double(riga4)
   #theta=double(riga5)
   #phi=double(riga6)
   #orbit=int(riga7)

   #grb=coords.Position((long,lat),system='galactic')
   grb=sky(l=long*u.degree,b=lat*u.degree,frame='galactic')

   #radec=grb.j2000()
   #ra=radec[0]
   #dec=radec[1]

   ra=grb.fk5.ra
   dec=grb.fk5.dec


   filepar=filename+".par"
   scritto=open(filepar,"w")
   scritto.write(str(orbit)+"\n")
   scritto.write(str(t0_obt)+" \n")
   scritto.write(str(ra)+" \n")
   scritto.write(str(dec)+" \n")
   scritto.write("15 \n")
   scritto.write("-20 \n")
   scritto.write("340 \n")
   scritto.write("-6000 \n")
   scritto.write("6000 \n")
   scritto.write("60 \n")
   scritto.write("80 \n")
   scritto.close()

   os.system("scp giuliani@agileatc2.iasf-roma.inaf.it:/AGILE_PROC3/F4_2/EVT/0"+str(orbit)+".evt.gz .")
   os.system("scp giuliani@agileatc2.iasf-roma.inaf.it:/AGILE_PROC3/DATA/LOG/0"+str(orbit)+".log.gz .")
   os.system("gunzip 0"+str(orbit)+".*gz")

   return str(orbit),filepar


#######  Script


if len(sys.argv) <= 2:

   files=get(sys.argv[1])

   orbit=files[0]
   filepar=files[1]

   find('0'+str(orbit)+'.evt','0'+str(orbit)+'.log',filepar)

else:
   find(sys.argv[1],sys.argv[2],sys.argv[3])

