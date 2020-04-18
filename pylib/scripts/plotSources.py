###
import sys, os
from numpy import *
from matplotlib.pyplot import *
from astropy.table import Table 
from astropy import units as u
from astropy.io import fits
###

from astropy.coordinates import SkyCoord, Angle
from astropy.table import  vstack



def plotSources(obsfile='', evtfile='', xmlfile='', ni=1, ns=100.):

  if obsfile != '': 
    
    obs=Table.read(obsfile,1)
    files=obs['EVENTS_FILENAME']
  
    tt=[]
    for file in files:
      tt=tt+[Table.read(file,1)]
    t=vstack(tt,metadata_conflicts='silent')
  else:
    t=Table.read(evtfile,1)
  
   
  phs= SkyCoord(ra=t['RA'], dec=t['DEC'] , unit='deg', frame='icrs')


  #target= SkyCoord(ra=270.12, dec=-24.03, unit='deg', frame='icrs')
  #sep=phs.separation(target)
  #t_on =  t[where(sep < 0.4*u.deg)]
  #hist(t_on['MC_ID'])


  for i in arange(ni,ns+1):
    if xmlfile != '':
      sourcename=inmodels[i-1].name()
    else:
      sourcename=str(i)
    w=where(t['MC_ID'] == i ) ; plot(phs.galactic.l[w],phs.galactic.b[w],'.',label=sourcename,markersize=2)
    print (len(w[0]),' events from source ',sourcename)
  legend()
  show()


######## Main 


if len(sys.argv) == 1:
  print
  print ("  -->  usalo cosi' : python plotSources.py NUltimaSorgente [NPrimaSorgente [modelli.xml]]")
  print

ns=int(sys.argv[1])

if len(sys.argv) > 2:
  ni=int(sys.argv[2])
else:
  ni=1


if len(sys.argv) > 3:
  import gammalib
  xmlfile=sys.argv[3]
  inmodels = gammalib.GModels(xmlfile)
  ns=min(len(inmodels),ns)  
  
print ('Plot components from',ni,'to',ns)

# Obsfile
#obsfile='obs-index.fits.gz'
#plotSources(obsfile=obsfile, xmlfile='', ni=ni, ns=ns)


# Evtfile 
evtfile='events.fits'
plotSources(evtfile=evtfile, xmlfile='', ni=ni, ns=ns)



