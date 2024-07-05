#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 18:03:19 2023

@author: andrea
"""

import sys, os
from numpy import *
from matplotlib.pyplot import *
from astropy.table import Table 
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord

from cats import fromGreen,fromManitoba, cat

pth = '/home/andrea/inafCloud/databases/catalogs/'


#%%


def fromFermiSNR(file = pth+'snr/1SC_catalog_v02.fits'):
    
    cdc = Table.read( file )
    cdc['SkyDir'] = SkyCoord(cdc['GLON'], cdc['GLAT'],  frame='galactic')
    cdc['Size'] = cdc['RADIUS']
    cdc['Name'] = cdc['SNR_Name']
    
    vere = cdc['Classification'] == 'Classified            '
    
    ccz = cat(cdc)    

    if 1:   ccz = cat(cdc[vere]) 


    return ccz



class SNR():
    
    def __init__(self, gname='G299.2-02.9'):
        
        self.gname=gname
        self.green = fromGreen().selectSources(fromName=gname)
        self.manitoba = fromManitoba().selectSources(fromName=gname)
       
        print(' Green matches : ',len(self.green.table))
        print(' Mani. matches : ',len(self.manitoba.table))
        
    
    def age(self, giveme='average',catalog='Manitoba'):
        
        
        manit={'minmax' : [self.manitoba.table['age_min (yr)'][0]*u.yr,
                          self.manitoba.table['age_max (yr)'][0]*u.yr],
              'average': (self.manitoba.table['age_min (yr)'][0]*u.yr+
                         self.manitoba.table['age_max (yr)'][0]*u.yr)/2.
              }
        
        ages = { 'Manitoba' : manit }
        
        
        return ages[catalog][giveme]
        
    
    
    def distance(self , giveme='average',catalog='Manitoba'):
        
        
        manit={'minmax' : [self.manitoba.table['distance_min (kpc)'][0]*u.kpc,
                          self.manitoba.table['distance_max (kpc)'][0]*u.kpc],
              'average': (self.manitoba.table['distance_min (kpc)'][0]*u.kpc+
                         self.manitoba.table['distance_max (kpc)'][0]*u.kpc)/2.
              }
        
        dist = { 'Manitoba' : manit }
        
        try : 
            dd = dist[catalog][giveme]
        except :
            print(    'Prova  ',  dist.keys())
            for ctlg in dist.keys():               
                print('      con  ', dist[ctlg].keys()  )
            dd = -1
        
        return dd
    
    
    def SkyDir(self):
        
        try :
            ss = self.green.table['SkyDir']
        except:
            ss = self.manitoba.table['SkyDir']

        
        return ss
        

class SNRset():
    
    def __init__(self, snrlist):
        
        self.snrlist = snrlist
        
    def Ages(self):
        
        self.ages = [] 
        
        for snr in self.snrlist:

            self.ages = self.ages + [ snr.age('average') ]
            
        return
    
    
    

        
#%%

G006 = SNR( gname='G006.4-00.1' )
G008 = SNR( gname='G008.7-00.1' )

W28 = SNR(gname='G006.4-00.1' )


sl = SNRset([ G006, G008])


sl.Ages()

            
#%%

superNR = Table()
ff = fromFermiSNR()
gnames = []

for fsnr in ff.table['Name']:
    
    gname='G'+fsnr[3:]
    print( gname )
    
    gnames = gnames + [gname]
    
    #q = SNR( gname=gname )
    
    # if len(q.manitoba.table) > 0:
    #     errorbar(q.gname[:4],
    #              q.age('average'),
    #              q.age('average')-q.age('minmax')[0] 
    #              )
    #     print( '  ',q.age('minmax') )
    
#yscale('log')

superNR['Name'] = gnames
superNR['Common Name'] = 'none        '
superNR['GeV'] = 1 == 1
superNR['TeV'] = 1==0

# superNR.write('superNR.ecsv' )


#%%

su = Table.read('superNR.ecsv')

#su['Age'] = W28.age(giveme='minmax')

#l = []


for s in su: 
    print(s['Name'])
    
    sn = SNR(s['Name'])

    plot( sn.age(), sn.distance(),'b.' )

    #l.append(SNR(s['Name']).age(giveme='minmax'))
    #s['Age'] = SNR(s['Name']).age()
      
        
for isnr in isnrs:     
    sn = SNR(isnr)
    plot( sn.age(), sn.distance(),'r.' )

    
xscale('log')


# su.write('superNR.ecsv' , overwrite = True)


#%%

su = Table.read('superNR.ecsv')


caption = 'List of SNRs firmly identified with a gamma-ray source. The Coloumns GeV, TeV and PeV indicate the detection respetively in the bands 0.1-100 GeV, 0.1-100 TeV and $>$ 0.1 PeV'

su.write('su.tex', overwrite=True  ,
         latexdict={'caption': caption } )






#%%


isnrs = ['G006.4-00.1','G008.7-00.1',
         'G023.3-00.3','G034.7-00.4',
         'G049.2-00.7','G089.0+04.7','G043.3-00.2',
         'G189.1+03.0','G348.5+00.1',
         'G349.7+00.2','G357.7-00.1'
         ]

for isnr in isnrs:
    
    q = SNR(gname=isnr)
    errorbar(q.gname[:4],
             q.age('average'),
             q.age('average')-q.age('minmax')[0]
             )
    print( q.age('minmax') )
    
    
    
yscale('log')
    
    
#%%



shell_snrs = ['G266.2-01.2',
              'G315.4-02.3',
              'G347.3-00.5',
              'G327.6+14.6']






SNR(shell_snrs[3]).age('minmax')















