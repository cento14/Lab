#!/usr/bin/env python3
"""
Created on Fri May  1 23:04:06 2020

@author: Andrea  ( andrea.giuliani@inaf.it )
"""

import sys, os
from numpy import *
from matplotlib.pyplot import *
from astropy.table import Table 
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord


irf = fits.open('irf_file.fits')

bkg=irf['BACKGROUND']

bkg.data[0][6] = bkg.data[0][6]  *1e-4

fits.update('irf_prova.fits',bkg.data,bkg.header,bkg.header['EXTNAME'])