#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 14:17:25 2023

@author: giuliani
"""

import sys, os
from numpy import *
from matplotlib.pyplot import *
from astropy.table import Table 
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord

from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun, get_moon

#from specsim import atmosphere


# In[]

def visInfo(source = SkyCoord(0.0,0.0,unit='deg',frame='galactic'), 
            sito = EarthLocation(lat=28.3012*u.deg, lon=-16.5082*u.deg, height=2000*u.m),
            time=Time('2024-1-1 00:00:00') ) :
    
    
    # sky at a given time  
    
    sito_locCoord = AltAz(obstime=time, location=sito)
    
    source_alt = source.transform_to(sito_locCoord).alt  
    
    sun = get_sun(time).transform_to(sito_locCoord)
    altezza_sole = sun.alt  
    
    moon = get_moon(time, location=sito)
    
    altezza_luna = moon.transform_to(sito_locCoord).alt  
    
    fase_luna = moon.separation(sun)
    dist_luna = moon.separation(source)
    
    rr = {'sun_alt'   : altezza_sole,  
          'source_alt': source_alt,  
          'moon_alt'  : altezza_luna,
          'moon_phase': fase_luna,
          'moon_sep'  : dist_luna 
          }
    
    return rr
    

def day2vis(source = SkyCoord(0.0,0.0,unit='deg',frame='galactic'), 
            sito = EarthLocation(lat=28.3012*u.deg, lon=-16.5082*u.deg, height=2000*u.m),
            day=Time('2024-1-1 00:00:00'), deltah=.25,
            za=40.*u.deg, moonThr=0.0*u.deg ):    
    
    # Estimates the visibility at a given day 
    
    ora = arange(0,24,deltah)

    time = day  + ora*u.hour
    sito_locCoord = AltAz(obstime=time, location=sito)
    
    source_alt = source.transform_to(sito_locCoord).alt  
    
    sun = get_sun(time).transform_to(sito_locCoord)
    altezza_sole = sun.alt  
    
    if moonThr < 90.*u.deg :
        moon = get_moon(time, location=sito)
        altezza_luna = moon.transform_to(sito_locCoord).alt  
        #fase_luna = moon.separation(sun)
        #print(fase_luna.deg)
    else :
        altezza_luna = 0.

    #plot(ora,source_alt,'g.',label='Source')
    #plot(ora,altezza_sole,'y.',label='Sun')
    #plot(ora,altezza_luna,'b.',label='Moon')    
    
    vis = (source_alt   >  (90. *u.deg -za)) *    \
          (altezza_sole < -12. *u.degree  ) *    \
          (altezza_luna <  moonThr  )

    orebuone = sum(vis) * deltah

    return orebuone*u.hr
    

site= {
       'Teide' :  EarthLocation(lat = 28.3012*u.deg, lon=-16.5082*u.deg, height=2360*u.m) ,
       'SLN'   :  EarthLocation(lat = 37.6933*u.deg, lon= 14.9747*u.deg, height=1730*u.m)
#      'SLN'   :  EarthLocation(lat = 37.6930*u.deg, lon= 14.9745*u.deg, height=1725*u.m)
}



#%%



































