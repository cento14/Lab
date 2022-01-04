#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 20 00:22:14 2021

@author: giuliani
"""

import sys, os
from numpy import *
from matplotlib.pyplot import *
from astropy.table import Table 
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord

#%%    GammaCat

def fromGammaCat(coord):
  cat2=Table.read(os.environ['HOME']+'/git/gamma-cat/output/gammacat.fits.gz')
  cat2_coord = SkyCoord(cat2['ra'], cat2['dec'],frame='icrs')
  idx, d2d, d3d = coord.match_to_catalog_sky(cat2_coord)
  return cat2[idx],d2d


#%%  LHAASO


# t = Table.read('lhaaso.txt', format='ascii')

# print('Name, Distance [deg], Association , Class')
# for s in t:
#     sc = SkyCoord(s['RA'],s['dec'],frame='icrs',unit='deg')
#     p,d = fromGammaCat(sc)
#     if d < 1*u.deg:
#       print( s['name'],d.deg , p['common_name'], p['classes'] )
      
      
      
#%%  Fermi

pth = os.environ['HOME']+'/inafCloud/databases/catalogs/'

def fermiCat(sourceName='4FGL J0000.3-7355'):
    
    fgl = Table.read( pth+ 'fgl4.fits', 1)
    s = ''
    
    for source in fgl:
        
        if source['Source_Name'].replace(' ','') == sourceName.replace(' ','') :
            s = source
    
    return s
    
