
from numpy import *
from matplotlib.pyplot import *
from astropy.table import Table 
from astropy import units as u

import gammalib, ctools, cscripts  #, gammapy

def golike(obs='obs.xml',model="model.xml",likeout='results.xml',data='data.xml'):
    
  sim=ctools.ctobssim()
  sim["inobs"]=obs     
  sim["inmodel"]=model
  sim["outevents"] = data
  sim['rad']=20.
  sim['seed']= int(random.random(1)[0]*1e5)
  sim["debug"]= True
  sim.execute()

  like  =  ctools.ctlike()
  like['inobs']    = data
  like['inmodel']  = model
  like['outmodel'] = likeout
  like["debug"]= True
  like.execute()

  return like.opt().value()

def createObs(outfile='obs.xml', name ='o1', id='0', lon=0.0, lat=0.0, emin=0.1e6, emax=150.0e6, tmin=0.0, tmax=3600.0, database='1dc', response='South_z20_50h', evtfile=''):
  code=['<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n',
  '<observation_list title="observation list"> \n',
  '   <observation name="iSNR_0" id="asdfg" instrument="CTA"> \n',
  '       <parameter name="Pointing" lon="888" lat="777" /> \n',
  '       <parameter name="EnergyBoundaries" emin="1e6" emax="160e6" /> \n',
  '       <parameter name="GoodTimeIntervals" tmin="0.0" tmax="5400" /> \n',
  '       <parameter name="TimeReference" mjdrefi="51544" mjdreff="0.5" timeunit="s" timesys="TT" timeref="LOCAL" /> \n',
  '       <parameter name="Deadtime" deadc="1.0" /> \n',
  '       <parameter name="Calibration" database="1dc" response="South_z20_50h" /> \n',
  '    </observation> \n',
  '</observation_list>\n',
  '\n']
  
  code[2] = code[2].replace('iSNR_0',name)
  code[2] = code[2].replace('asdfg',id)
  code[3] = code[3].replace('888',str(lon))
  code[3] = code[3].replace('777',str(lat))
  code[4] = code[4].replace('1e6'  ,str(emin))
  code[4] = code[4].replace('160e6',str(emax))
  code[5] = code[5].replace('0.0' ,str(tmin))
  code[5] = code[5].replace('5400',str(tmax))
  code[8] = code[8].replace('1dc',database)
  code[8] = code[8].replace('South_z20_50h',response)
  
  
  if evtfile != '' :
      evtlist='       <parameter name="EventList" file="'+evtfile+'" /> \n'
      code.insert(3,evtlist)
  
  print ('Writing file : '+outfile)
  obs = open(outfile, 'w')
  obs.writelines(code)
  obs.close()
  
  return code


def createModel(outfile='model.xml',specfile='specfile.txt', mapfile='mapfile.fits'):
  code=['<source_library title="source library"> \n',
  '   <source name="source0" type="DiffuseSource"  tscalc="1">\n',
  '       <spectrum file="model0.txt" type="FileFunction">\n',
  '           <parameter free="0" max="1000.0" min="0.0" name="Normalization" scale="1.0" value="1.0" />\n',
  '       </spectrum>\n',
  '       <spatialModel type="DiffuseMap" file="skymap0.fits">\n',
  '           <parameter name="Prefactor" value="1" scale="1" min="0.001" max="1000" free="0" />\n',
  '       </spatialModel>\n',
  '   </source>\n',
  '   <source name="CTABackgroundModel" type="CTAIrfBackground" instrument="CTA">\n',
  '    <spectrum type="PowerLaw">\n',
  '      <parameter name="Prefactor"   scale="1.0"  value="1.0"  min="1e-3" max="1e+3"   free="0"/>\n',
  '      <parameter name="Index"       scale="1.0"  value="0.0"  min="-5.0" max="+5.0"   free="0"/>\n',
  '      <parameter name="PivotEnergy" scale="1e6"  value="1.0"  min="0.01" max="1000.0" free="0"/>\n',
  '    </spectrum>\n',
  '</source>   \n',
  '</source_library>\n']

  code[2] = code[2].replace('model0.txt',specfile)
  code[5] = code[5].replace('skymap0.fits',mapfile)

  print ('Writing file : '+outfile)
  mod = open(outfile, 'w')
  mod.writelines(code)
  mod.close()
  
  return code


def createModel2(outfile='model.xml',specfile='specfile.txt', lon=0., lat=0.):
  code=['<source_library title="source library"> \n',
  '   <source name="source0" type="ExtendedSource""  tscalc="1">\n',
  '       <spectrum file="model0.txt" type="FileFunction">\n',
  '           <parameter free="0" max="1000.0" min="0.0" name="Normalization" scale="1.0" value="1.0" />\n',
  '       </spectrum>\n',
  '       <spatialModel type="PointSource">\n',
  '          <parameter name="LON"  scale="1.0" value="8888" min="-360" max="360" free="1"/>     \n'
  '          <parameter name="LAT" scale="1.0" value="8888" min="-90"  max="90"  free="1"/>     \n',
  '       </spatialModel>\n',
  '   </source>\n',
  '   <source name="CTABackgroundModel" type="CTAIrfBackground" instrument="CTA">\n',
  '    <spectrum type="PowerLaw">\n',
  '      <parameter name="Prefactor"   scale="1.0"  value="1.0e-3"  min="1e-3" max="1e+3"   free="0"/>\n',
  '      <parameter name="Index"       scale="1.0"  value="0.0"  min="-5.0" max="+5.0"   free="0"/>\n',
  '      <parameter name="PivotEnergy" scale="1e6"  value="1.0"  min="0.01" max="1000.0" free="0"/>\n',
  '    </spectrum>\n',
  '</source>   \n',
  '</source_library>\n']

  code[2] = code[2].replace('model0.txt',specfile)
  code[6] = code[6].replace('8888',str(lon))
  code[7] = code[7].replace('8888',str(lat))

  print ('Writing file : '+outfile)
  mod = open(outfile, 'w')
  mod.writelines(code)
  mod.close()
  
  return code


def createModel3(outfile='model.xml',lon=0.0,lat=0.0,specfile='specfile.txt', mapfile='mapfile.fits'):  
  code=['<source_library title="source library"> \n',
  '   <source name="source0" type="DiffuseSource"  tscalc="1">\n',
  '       <spectrum file="model0.txt" type="FileFunction">\n',
  '           <parameter free="1" max="2.0" min="0.5" name="Normalization" scale="1.0" value="1.0" />\n',
  '       </spectrum>\n',
  '      <spatialModel type="RadialGaussian">\n',
  ' 		 <parameter name="GLON"    scale="1.0" value="888" min="-360" max="360" free="0"/>\n',
  ' 		 <parameter name="GLAT"   scale="1.0" value="777" min="-90"  max="90"  free="0"/>\n',
  '  	 	 <parameter name="Sigma" scale="1.0" value="0.20"    min="0.01" max="10"  free="1"/>\n',
  '	</spatialModel>\n',  
  '   </source>\n',
  '   <source name="CTABackgroundModel" type="CTAIrfBackground" instrument="CTA">\n',
  '    <spectrum type="PowerLaw">\n',
  '      <parameter name="Prefactor"   scale="1.0"  value="1.0"  min="1e-3" max="1e+3"   free="1"/>\n',
  '      <parameter name="Index"       scale="1.0"  value="0.0"  min="-5.0" max="+5.0"   free="0"/>\n',
  '      <parameter name="PivotEnergy" scale="1e6"  value="1.0"  min="0.01" max="1000.0" free="0"/>\n',
  '    </spectrum>\n',
  '</source>   \n',
  '</source_library>\n']

  code[2] = code[2].replace('model0.txt',specfile)
  code[6] = code[6].replace('888',str(lon))
  code[7] = code[7].replace('777',str(lat))
  #code[5] = code[5].replace('skymap0.fits',mapfile)

  print ('Writing file : '+outfile)
  mod = open(outfile, 'w')
  mod.writelines(code)
  mod.close()
  
  return code


class xmlModel:                                                                                                                                                         
  def __init__(self, file_name):                                                                                                                                        
    self.models = gammalib.GModels(file_name)                                                                                                                          
                                                                   
  def PL(self, sname, plotta=True):
     k     = self.models[sname]["Prefactor"].value()                                                                                                                         
     index = self.models[sname]["Index"].value()                                                                                                                             
     e_piv = self.models[sname]["PivotEnergy"].value()                                                                                                                       
                                                                                                                                                            
     en = 10**( arange(110)/20. ) * 1e3                # MeV                                                                                                         
     f = k * (en / e_piv)**(index)                        # ph / cm2 s MeV                                                                                             
                                                              
     if plotta: 
       loglog()                                                                                                                                                       
       plot(en*1e6,f*en**2.*1.6e-6)                                                                                                                               
                                                              
     return en*1e6*u.eV,f*en**2.*1.6e-6*u.TeV                                                                                                                                                    
                                                              
  def PLC(self, sname, plotta=True):
     k     = self.models[sname]["Prefactor"].value()                                                                                                                         
     index = self.models[sname]["Index"].value()                                                                                                                             
     e_piv = self.models[sname]["PivotEnergy"].value()                                                                                                                       
     ecut  = self.models[sname]["CutoffEnergy"].value()                                                                                                                      
                                                                                                                                                            
     en = 10**( arange(110)/20. ) * 1e3                # MeV                                                                                                         
     f = k * (en / e_piv)**(index) * exp(-en / ecut)    # ph / cm2 s MeV                                                                                             
        
     if plotta: 
       loglog()                                                                                                                                                       
       plot(en*1e6,f*en**2.*1.6e-6)                                                                                                                               
                                                              
     #fig.show()                                                                                                                                                        
     return en*1e6*u.eV,f*en**2.*1.6e-6*u.TeV                                                      


def stackedPipeline(obsfile='index.xml', l=0.01, b=0.01, emin=0.1, emax=100.0,
                 enumbins=20, nxpix=200, nypix=200, binsz=0.02,
                 coordsys='CEL', proj='CAR', caldb='prod2', irf='acdc1a',debug=False,inmodel='Crab.xml'):

    # da esempi ctools 1.4.2

    """
    Simulation and stacked analysis pipeline

    Parameters
    ----------
    obs : `~gammalib.GObservations`
        Observation container
    ra : float, optional
        Right Ascension of counts cube centre (deg)
    dec : float, optional
        Declination of Region of counts cube centre (deg)
    emin : float, optional
        Minimum energy (TeV)
    emax : float, optional
        Maximum energy (TeV)
    enumbins : int, optional
        Number of energy bins
    nxpix : int, optional
        Number of pixels in X axis
    nypix : int, optional
        Number of pixels in Y axis
    binsz : float, optional
        Pixel size (deg)
    coordsys : str, optional
        Coordinate system
    proj : str, optional
        Coordinate projection
    debug : bool, optional
        Debug function
    """

    # Bin events into counts map
    bin = ctools.ctbin()
    bin['inobs']= obsfile
    bin['ebinalg']  = 'LOG'
    bin['emin']     = emin
    bin['emax']     = emax
    bin['enumbins'] = enumbins
    bin['nxpix']    = nxpix
    bin['nypix']    = nypix
    bin['binsz']    = binsz
    bin['coordsys'] = coordsys
    bin['proj']     = proj
    bin['xref']     = l
    bin['yref']     = b
    bin['debug']    = debug
    bin.execute()
    print 'Datacube : done!'


    # Create exposure cube
    expcube = ctools.ctexpcube()
    #expcube['incube']=bin.obs()
    expcube['inobs']= obsfile
    expcube['incube']   = 'NONE'
    expcube['ebinalg']  = 'LOG'
    expcube['emin']     = emin
    expcube['emax']     = emax
    expcube['enumbins'] = enumbins
    expcube['nxpix']    = nxpix
    expcube['nypix']    = nypix                                                                                  
    expcube['binsz']    = binsz                                                                                  
    expcube['coordsys'] = coordsys                                                                               
    expcube['proj']     = proj                                                                                   
    expcube['xref']     = l                                                                                    
    expcube['yref']     = b                                                                                    
    expcube['debug']    = debug                                                                                  
    expcube.execute()                                                                                                
    print 'Expcube : done!'  
 
    # Create PSF cube                                                                                            
    psfcube = ctools.ctpsfcube()
    psfcube['inobs']= obsfile
    psfcube['incube']   = 'NONE'                                                                                 
    psfcube['ebinalg']  = 'LOG'                                                                                  
    psfcube['emin']     = emin                                                                                   
    psfcube['emax']     = emax                                                                                   
    psfcube['enumbins'] = enumbins                                                                               
    psfcube['nxpix']    = 10                                                                                     
    psfcube['nypix']    = 10                                                                                     
    psfcube['binsz']    = 1.0                                                                                    
    psfcube['coordsys'] = coordsys                                                                               
    psfcube['proj']     = proj                                                                                   
    psfcube['xref']     = l                                                                                     
    psfcube['yref']     = b                                                                                    
    psfcube['debug']    = debug                                                                                  
    psfcube.execute()                                                                                                
    print 'Psfcube : done!'
    
    edispcube = ctools.ctedispcube()
    edispcube['inobs']    =  obsfile
    edispcube['caldb']    = caldb
    edispcube['irf']      = irf
    edispcube['incube']   = 'NONE'
    edispcube['xref']     = l
    edispcube['yref']     = b
    edispcube['proj']     = proj
    edispcube['coordsys'] = coordsys
    edispcube['binsz']    = 1.0   # deg/bin; the energy dispersion only varies slowly
    edispcube['nxpix']    = 10
    edispcube['nypix']    = 10
    edispcube['emin']     = emin
    edispcube['emax']     = emax
    edispcube['enumbins'] = enumbins
    #edispcube['outcube']  = 'edispcube.fits'
    edispcube.execute()
    print 'Edispcube : done!'

    # Create background cube                                                                                     
    bkgcube = ctools.ctbkgcube()
    bkgcube['inobs']= obsfile
    bkgcube['incube']   = 'NONE'
    bkgcube['ebinalg']  = 'LOG'
    bkgcube['emin']     = emin
    bkgcube['emax']     = emax
    bkgcube['enumbins'] = enumbins
    bkgcube['nxpix']    = 10
    bkgcube['nypix']    = 10
    bkgcube['binsz']    = 1.0
    bkgcube['coordsys'] = coordsys
    bkgcube['proj']     = proj
    bkgcube['xref']     = l
    bkgcube['yref']     = b
    bkgcube['debug']    = debug
    bkgcube['inmodel'] = inmodel
    bkgcube.execute()
    print 'Bkgcube : done!'
    
    # Attach background model to observation container
    bin.obs().models(bkgcube.models())

    # Set Exposure and Psf cube for first CTA observation
    # (ctbin will create an observation with a single container)
    bin.obs()[0].response(expcube.expcube(), psfcube.psfcube(), edispcube.edispcube(), bkgcube.bkgcube())

    # Perform maximum likelihood fitting
    like = ctools.ctlike(bin.obs())
    like['inmodel']=inmodel
    like['edisp'] = True   
    like['debug'] = True # Switch this always on for results in console
    #like['statistic']='CSTAT'
    like.run()

    # Return
    return



