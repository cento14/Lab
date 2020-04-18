
import os
from astropy import units as u
from astropy.coordinates import SkyCoord as coords
from astropy.io import fits
from ctoolsfordummy import createObs

os.system('ls evt*fits > fitsfiles')

lista=open('fitsfiles')
righe=lista.readlines()

code= ['<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n',
       '<observation_list title="observation list"> \n'] 

for riga in righe:
  
  evtfile=riga[0:-1]
  file_p = fits.open(evtfile)

  ra_pnt=file_p[1].header['RA_PNT'] 
  dec_pnt=file_p[1].header['DEC_PNT']

  dir_onaxis=coords(ra=ra_pnt*u.degree,dec=dec_pnt*u.degree,frame='fk5')

  lon=dir_onaxis.galactic.l.value
  lat=dir_onaxis.galactic.b.value

  tmin= file_p[1].header['TSTART']
  tmax= file_p[1].header['TSTOP']

  print 'Aggiungo file : '+evtfile,'  puntato in (l,b):',lon,lat

  database= file_p[1].header['CALDB']
  response= file_p[1].header['IRF']
  
  obsid= file_p[1].header['OBS_ID']
  
  cc=createObs(name=evtfile[0:-5], evtfile=evtfile,lon=lon, lat=lat, tmin=tmin, tmax=tmax, id=str(obsid), outfile='.cancellami', database=database, response=response)

  code=code+cc[2:-2]

code=code+['</observation_list>\n','\n']


print ('Writing file : '+'index.xml')
indx=open('index.xml','w')
indx.writelines(code)
indx.close()
