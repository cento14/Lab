

import os, sys
from astropy.table import Table 
from gammapyfordummy import get_events_file_info

#irfpath='../../../caldb/data/cta/prod2/bcf/Astri_50h_June/'
irfpath=os.path.relpath(os.path.expandvars('$CALDB/data/cta/prod2/bcf/Astri_50h_June/'))

### Main ##

os.system('ls ev*fits > fitsfiles')
filename='fitsfiles'
#filename=sys.argv[1]

lista=open(filename)
righe=lista.readlines()

#t=Table.read('hdu-index.fits.gz')
#t.dtype

q=[('OBS_ID', '>i8'), ('RA_PNT', '>f8'), ('DEC_PNT', '>f8'), ('GLON_PNT', '>f8'), ('GLAT_PNT', '>f8'), ('ZEN_PNT', '>f8'), ('ALT_PNT', '>f8'), ('AZ_PNT', '>i8'), ('ONTIME', '>f8'), ('LIVETIME', '>f8'), ('DEADC', '>f8'), ('TSTART', '>f8'), ('TSTOP', '>f8'), ('DATE_OBS', 'S10'), ('TIME_OBS', 'S8'), ('DATE_END', 'S10'), ('TIME_END', 'S8'), ('N_TELS', '>i8'), ('OBJECT', 'S21'), ('CALDB', 'S6'), ('IRF', 'S13'), ('EVENTS_FILENAME', 'S50'), ('EVENT_COUNT', '>i8')]

q2=[('OBS_ID', '>i8'), ('HDU_TYPE', 'S6'), ('HDU_CLASS', 'S10'), ('FILE_DIR', 'S142'), ('FILE_NAME', 'S50'), ('HDU_NAME', 'S21')]

obsIndex=Table(dtype=q)
hduIndex=Table(dtype=q2)

for riga in righe:
  
  evtfile=riga[0:-1]
  
  info=get_events_file_info(evtfile, irfpath=irfpath)
  
  obsIndex.add_row()
  obsIndex['OBS_ID'][-1]   =info['OBS_ID'] 
  obsIndex['RA_PNT'][-1]   =info['RA_PNT']
  obsIndex['DEC_PNT'][-1]  =info['DEC_PNT']
  obsIndex['GLON_PNT'][-1] =info['GLON_PNT']
  obsIndex['GLAT_PNT'][-1] =info['GLAT_PNT']
  obsIndex['ZEN_PNT'][-1]  =info['ZEN_PNT']
  obsIndex['ALT_PNT'][-1]  =info['ALT_PNT']
  obsIndex['AZ_PNT'][-1]   =info['AZ_PNT']
  obsIndex['ONTIME'][-1]   =info['ONTIME']
  obsIndex['LIVETIME'][-1]  =info['LIVETIME']
  obsIndex['DEADC'][-1]  =info['DEADC']
  obsIndex['TSTART'][-1]  =info['TSTART']
  obsIndex['TSTOP'][-1]  =info['TSTOP']
  #obsIndex['DATE_OBS'][-1]  =info['DATE_OBS']
  #obsIndex['TIME_OBS'][-1]  =info['TIME_OBS']
  #obsIndex['DATE_END'][-1]  =info['DATE_END']
  #obsIndex['TIME_END'][-1]  =info['TIME_END']
  obsIndex['N_TELS'][-1]  =info['N_TELS']
  obsIndex['OBJECT'][-1]  =info['OBJECT']
  obsIndex['CALDB'][-1]  =info['CALDB']
  obsIndex['IRF'][-1]  =info['IRF']
  obsIndex['EVENTS_FILENAME'][-1]  =info['EVENTS_FILENAME']
  obsIndex['EVENT_COUNT'][-1]  =info['EVENT_COUNT']
  
  print ('Aggiungo file : '+evtfile,'  puntato in (l,b):',info['GLON_PNT'],info['GLAT_PNT'])

  hduIndex.add_row()    #  Events
  hduIndex['OBS_ID'][-1]   =info['OBS_ID']
  hduIndex['HDU_TYPE'][-1] ='events'
  hduIndex['HDU_CLASS'][-1]='events'
  hduIndex['FILE_DIR'][-1] =info['EVENTS_PATH']
  hduIndex['FILE_NAME'][-1]=info['EVENTS_FILENAME']
  hduIndex['HDU_NAME'][-1]='EVENTS'
  
  hduIndex.add_row()    #  GTI
  hduIndex['OBS_ID'][-1]   =info['OBS_ID']
  hduIndex['HDU_TYPE'][-1] ='gti'
  hduIndex['HDU_CLASS'][-1]='gti'
  hduIndex['FILE_DIR'][-1] =info['EVENTS_PATH']
  hduIndex['FILE_NAME'][-1]=info['EVENTS_FILENAME']
  hduIndex['HDU_NAME'][-1]='GTI'
  
  hduIndex.add_row()    #  Aeff
  hduIndex['OBS_ID'][-1]    =info['OBS_ID']
  hduIndex['HDU_TYPE'][-1]  ='aeff'
  hduIndex['HDU_CLASS'][-1] ='aeff_2d'
  hduIndex['FILE_DIR'][-1]  = info['IRF_PATH']
  hduIndex['FILE_NAME'][-1] = info['IRF_FILENAME']
  hduIndex['HDU_NAME'][-1]  ='EFFECTIVE AREA'
  
  hduIndex.add_row()    #  EDisp
  hduIndex['OBS_ID'][-1]    =info['OBS_ID']
  hduIndex['HDU_TYPE'][-1]  ='edisp'
  hduIndex['HDU_CLASS'][-1] ='edisp_2d'
  hduIndex['FILE_DIR'][-1]  = info['IRF_PATH']
  hduIndex['FILE_NAME'][-1] = info['IRF_FILENAME']
  hduIndex['HDU_NAME'][-1]  ='ENERGY DISPERSION'
  
  hduIndex.add_row()    #  Psf
  hduIndex['OBS_ID'][-1]    =info['OBS_ID']
  hduIndex['HDU_TYPE'][-1]  ='psf'
  hduIndex['HDU_CLASS'][-1] ='psf_3gauss'
  hduIndex['FILE_DIR'][-1]  = info['IRF_PATH']
  hduIndex['FILE_NAME'][-1] = info['IRF_FILENAME']
  hduIndex['HDU_NAME'][-1]  ='POINT SPREAD FUNCTION'
 
  hduIndex.add_row()    #  Background
  hduIndex['OBS_ID'][-1]    =info['OBS_ID']
  hduIndex['HDU_TYPE'][-1]  ='bkg'
  hduIndex['HDU_CLASS'][-1] ='bkg_3d'
  hduIndex['FILE_DIR'][-1]  = info['IRF_PATH']
  hduIndex['FILE_NAME'][-1] = info['IRF_FILENAME']
  hduIndex['HDU_NAME'][-1]  ='BACKGROUND'

  

obsIndex.meta={'MJDREFI': 51544,'MJDREFF': 0.5}
obsIndex.write('obs-index.fits',format='fits')
hduIndex.write('hdu-index.fits',format='fits')

os.system('gzip *-index.fits ')























