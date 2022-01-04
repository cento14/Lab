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

import regions

      
#%%   cat e' una table con una colonna 'SkyDir' di tipo skycoord, una 'Name' e una 'Size'


class cat(): 

    def __init__(self, cd):
        self.table = cd

    
    def toReg(self, regfile='pippo.reg'):
    
        cc = self.table
        regs=[]
        
        for source in cc : 
                        
            try :  
              radius = u.Quantity( nan_to_num(source['Size']), cc['Size'].unit)
              radius = max(radius ,  0.01 * u.deg)
            except:
              radius= 0.01 * u.deg
    
            reg= regions.CircleSkyRegion( source['SkyDir'], radius)

            reg.visual['color']='Yellow'
        
            if source['Name'] != '' :
                reg.meta['label']= source['Name']
            else : 
                print('No Name')
        
            regs = regs + [reg] 
             
        print ('  Writing file : ',regfile) 
        regions.write_ds9(regs,regfile) 
        
    
    
    def selectSources(self,   center = SkyCoord( 0,0, unit='deg', frame='galactic') , radius = 10*u.deg, 
                      fromName='' , 
                      fromDir=0
                      ):
        
        
        
        tab = self.table
        lt = len(tab)
        ii=[]

        if fromName != '' :           #  macthes from 1 name
            for i  in arange(lt):
                if tab[i]['Name'].replace(' ','') == fromName.replace(' ','') :
                    ii = ii + [i]

        elif fromDir != 0 :           #  macths from a list of directions
            ii, d2d, d3d = fromDir.match_to_catalog_sky(tab['SkyDir'])
    
        else :                         # mactes from ROI
            b = center.separation(tab['SkyDir'])
            ii = where(b < radius)
    
        newtab = tab[ii]
        newccz = cat(newtab)
        
        return newccz
        
                                                                                                                                                                 
  

      
#%%  Fermi

pth = os.environ['HOME']+'/inafCloud/databases/catalogs/'

def fermiCat(sourceName='4FGL J0000.3-7355'):
    
    fgl = Table.read( pth+ 'fgl4.fits', 1)
    s = ''
    
    for source in fgl:
        
        if source['Source_Name'].replace(' ','') == sourceName.replace(' ','') :
            s = source
    
    return s
    

def fromFermi( file = pth + 'fgl4.fits' ):
    
    fgl = Table.read( file, 1)
    
    ss = SkyCoord( fgl['RAJ2000'], fgl['DEJ2000'], frame='icrs'  )
    
    cdc = fgl
    cdc['Name'] = fgl['Source_Name']
    cdc['SkyDir'] = ss
    cdc['Size'] = fgl['Conf_68_SemiMajor']
    
    ccz = cat(cdc)
    return ccz
    



#%%    GammaCat

catfile = os.environ['HOME']+'/git/gamma-cat/output/gammacat.fits.gz'

def fromGammaCat(file = catfile):
  
    cdc=Table.read(file)
    
    cdc['Name'] = cdc['common_name']
    cdc['SkyDir'] = SkyCoord(cdc['ra'], cdc['dec'], frame='icrs')
    cdc['Size'] = cdc['morph_sigma']
    
    ccz = cat(cdc)    

    return ccz


# def fromGammaCat(coord):
#   cat2=Table.read(os.environ['HOME']+'/git/gamma-cat/output/gammacat.fits.gz')
#   cat2_coord = SkyCoord(cat2['ra'], cat2['dec'],frame='icrs')
#   idx, d2d, d3d = coord.match_to_catalog_sky(cat2_coord)
#   return cat2[idx],d2d




#%%  LHAASO


# t = Table.read('lhaaso.txt', format='ascii')

# print('Name, Distance [deg], Association , Class')
# for s in t:
#     sc = SkyCoord(s['RA'],s['dec'],frame='icrs',unit='deg')
#     p,d = fromGammaCat(sc)
#     if d < 1*u.deg:
#       print( s['name'],d.deg , p['common_name'], p['classes'] )
      
      












