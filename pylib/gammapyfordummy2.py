from numpy import *
from matplotlib.pyplot import *
from gammapy.irf import Background3D                     
from astropy.table import vstack as table_vstack
from astropy import units as u
from astropy.table import Table 
from astropy.coordinates import SkyCoord #as coords
from astropy.io import fits
from collections import OrderedDict


def showbkg(irf_filename='irf_file.fits'):

  bkg = Background3D.read(irf_filename, hdu="BACKGROUND")
  tb=bkg.to_table() 

  for i  in arange(len(tb['ENERG_HI'][0])):
  
    #print 'Energy range : ',tb['ENERG_LO'][0][i],tb['ENERG_HI'][0][i]
    title(str(round(tb['ENERG_LO'][0][i],2))+' - '+str(round(tb['ENERG_HI'][0][i],2))+'  TeV')
    imshow(tb['BKG'][0][i])
    pause(1)
  
  return 1

def stackEvents(obsindex='obs-index.fits.gz'):

  inda= Table.read(obsindex)

  eev=[]
  for file in inda['EVENTS_FILENAME']:
    eev=eev+[Table.read(file)]
  
  evv=table_vstack(eev, metadata_conflicts="silent")
  
  return evv


def get_events_file_info(evtfile,evtpath='./',irffile='irf_file.fits',irfpath='$CALDB/'):
    #log.debug(f'Reading {filename}')
    header = fits.open(evtfile)['EVENTS'].header
    events = fits.open(evtfile)['EVENTS'].data
    
    gti = fits.open(evtfile)['GTI'].data
    
    info = OrderedDict()
    #
    info['OBS_ID'] = header['OBS_ID']
    info['RA_PNT'] = header['RA_PNT']
    info['DEC_PNT'] = header['DEC_PNT']
    #
    pos = SkyCoord(info['RA_PNT'], info['DEC_PNT'], unit='deg').galactic
    info['GLON_PNT'] = pos.l.deg
    info['GLAT_PNT'] = pos.b.deg
    #
    info['ZEN_PNT'] = 90 - float(header['ALT_PNT'])
    info['ALT_PNT'] = header['ALT_PNT']
    info['AZ_PNT'] = header['AZ_PNT']
    info['ONTIME'] = header['ONTIME']
    info['LIVETIME'] = header['LIVETIME']
    info['DEADC'] = header['DEADC']
    info['TSTART'] = header['TSTART']
    info['TSTOP'] = header['TSTOP']
    
    try :
      info['DATE_OBS'] = header['DATE_OBS']
      info['DATE_END'] = header['DATE_END']
    except:
      print ' DATE not defined'

    try :
     info['TIME_OBS'] = header['TIME_OBS']
     info['TIME_END'] = header['TIME_END']
    except:
      print ' TIME_OBS not defined'
      
    
    info['N_TELS'] = header['N_TELS']
    info['OBJECT'] = header['OBJECT']
    info['CALDB'] = header['CALDB']
    info['IRF'] = header['IRF']
    #
    # Not part of the spec, but good to know from which file the info comes
    info['IRF_FILENAME']='irf_file.fits'
    info['IRF_PATH']=irfpath
    info['EVENTS_FILENAME'] = evtfile
    info['EVENTS_PATH']='./'
    #
    # import IPython; IPython.embed(); 1/0
    info['EVENT_COUNT'] = header['NAXIS2']
    #
    info['EVENT_TIME_MIN'] = events['TIME'].min()
    info['EVENT_TIME_MAX'] = events['TIME'].max()
    info['EVENT_ENERGY_MIN'] = events['ENERGY'].min()
    info['EVENT_ENERGY_MAX'] = events['ENERGY'].max()
    #
    # gti = Table.read(filename, hdu='GTI')
    info['GTI_START'] = gti['START'][0]
    info['GTI_STOP'] = gti['STOP'][0]
    #
    return info






