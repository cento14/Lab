
###  
import sys, os
from numpy import *
from matplotlib.pyplot import *
from astropy.table import Table 
from astropy import units as u
from astropy.io import fits
###

ns=int(sys.argv[1])

obs=Table.read('obs-index.fits.gz')
obsfiles=obs['EVENTS_FILENAME']


def removeSource(filin, source='', ns=0., output='evt_cleaned.fits'):
  evt=fits.open(filin)
  t=evt[1].data
  tn=t[where(t['MC_ID'] != ns )[0]]
  evt[1].data=tn
  evt.writeto(output)
  evt.close()
  return


for evtfile in obsfiles:
  output=evtfile.replace('.fits','_senza'+str(ns)+'.fits')
  print 'Removing source ',ns,'from ',evtfile,'   -->   writing ',output
  removeSource(evtfile,ns=ns,output=output)

