#! /usr/bin/env python
# filename: gsed.py
# authors: Andrea Giuliani <giuliani@iasf-milano.inaf.it>
#          mario <mario@piffio.org>
# version: v3
# date: 2018-02-16
# description:
#   This library implements basic classes for manipulating
#   FITS catalogs of TeV sources for CTA.
# changelog:
#   v0 - entry version (classes catalog and source)
#   v1 - dump to XML file
#   v2 - prepend the index file path and background component
#   v3 - added mkTable() method 

import sys, os
import xml.etree.ElementTree as ET
import astropy.coordinates as apc
from astropy.table import Table
from astropy.io import fits as pf
import matplotlib.pyplot as fig
from numpy import *
import numpy as np

__version__ = '2020-02-16'

class catalog():
  def __init__(self, file_name = 'index.fits'):
    self.index = Table.read(file_name)
    self.path = os.path.dirname(os.path.abspath(file_name))
  
  
  def getSource(self, n):
    src=self.index[n]
    srcfile = (self.path+'/'+src['PATH']+'/'+src['Data_File']).replace(' ','')
    s=source(file_name=srcfile)
    return s
  
  def find(self,patt):
    for i in arange(len(self.index)):
      if (self.index['SNR_Name'][i].find(patt) != -1)+(self.index['Other_Names'][i].find(patt) != -1) :
        print( i, self.index['SNR_Name'][i], self.index['Other_Names'][i])
        
  
  
  # you can pass either celestial or Galactic coordinates, as
  # named arguments, or celestial, if positional
  def ToXML(self, ra = -361., dec = 0., roi = 2., l = -361., b = 0., \
      skydir = 0., XML_file_name='model.xml', add_background=True):
    if ra != -361.:
      skydir = apc.SkyCoord(ra, dec, unit='deg', frame='icrs')
    elif l != -361.:
      skydir = apc.SkyCoord(l, b, unit='deg', frame='galactic')

    if skydir == 0.:
      print('Error: Invalid direction')
      return False
    print('Select sources within %.4f deg from (RA, Dec) = (%.6f, %.6f)' % \
        (roi, skydir.icrs.ra.value, skydir.icrs.dec.value))

    # init the XML tree
    src_lib = ET.Element('source_library', {'title': 'source library'})
    for src_entry in self.index:
      src_file_name = '%s/%s/%s' % (self.path, \
          src_entry['PATH'].replace(' ', ''), \
          src_entry['Data_File'].replace(' ', ''))
      src_name = src_entry['TeVCat'].replace(' ', '')
      Src = source(file_name = src_file_name, name = src_name, \
        ra = src_entry['ra'], dec = src_entry['dec'])
      angsep_deg = Src.skydir.separation(skydir).value
      if angsep_deg <= roi:
        src_string = 'source %s at %.4f deg from (RA, Dec) = (%.6f, %.6f)' % \
            (Src.name, angsep_deg, skydir.icrs.ra.value, skydir.icrs.dec.value)
        src_model = src_entry['MODEL'].replace(' ', '')
        # skip sources without a defined model
        if not src_model:
          print('Warning: missing model for %s' % src_string)
          continue
        # skip sources for which either the model cannot be read
        # or the file cannot be written
        if not Src.ToXML(src_lib, ext_model = src_model):
          print('Warning: problems adding %s' % src_string) 
          continue
        # yey!! a new source
        print('Adding %s' % src_string)
    if add_background:
      print('Adding the background')
      background('IRFs').ToXML(src_lib)
      
    indent(src_lib)
    #print ET.tostring(src_lib, method='html')
    ET.ElementTree(src_lib).write(XML_file_name)

  def mkhtml(self,regions=''):
    lista=self.index.copy()
 
    bb=lista['b']
    ll=lista['l']
    #lista2=lista
 
    lista['l']= ll*(ll <= 180) + (ll-360.)*(ll > 180)
    lista.sort('l')
 
    ll=lista['l']
 
    lista['Fits Data_File']='                                                                               '
    for s in arange(len(ll)):
       if lista['Data_File'][s] != 'none':
         lista['Fits Data_File'][s]='*(*a href="data/'+lista['Data_File'][s]+'"*)* Download *(*/a*)*'
       else:
          lista['Fits Data_File'][s]='Not yet available'
 
 
    lista['Plot']='                                                                               '
    for s in arange(len(ll)):
       if lista['Data_File'][s] != 'none':
         lista['Plot'][s]='*(*a href="data/'+lista['Data_File'][s].replace(' ','')+'.png'+'" target="_new" *)* Show *(*/a*)*'
       else:
          lista['Plot'][s]='Not yet available'
 
 
    ### Galactic regions
    #  Scutum   ~30   (40,15)
    #  GC        0     (15,-10)
    #  3 kpc    ~337   (-10,-35)
    #  Crux     ~317   (-35,-60)
    #  Carina   ~285   (-60,-85)
    #  Vela     ~265   (-85,-110)
 
 
   
    mw=regions
    #mw={'Scutum Arm':[40,15], 'Galactic Center':[15,-10], '3-kpc Arm':[-10,-35],'Crux Arm':[-35,-60],'Carina Arm':[-60,-85],'Vela Arm':[-85,-110]}
    #mw={'Southern Milky Way':[40,-110]}
 
    colonnedascrivere=['TeVCat','l','b','type','Other_Names','Fits Data_File','Plot']
 
    wfatte=where(lista['Data_File'] != 'none')
    fig.plot(ll,bb,'ro',ll[wfatte],bb[wfatte],'go',markersize=8)
 
    for s in arange(len(ll)):
      fig.text(ll[s]-.5,bb[s]*1.2,lista['TeVCat'][s],fontsize=8)
 
    for region in mw:
      w=where((ll > mw[region][1]) * (ll < mw[region][0]))
      htmlfile=region.replace(' ','')+'.html'
      lista[colonnedascrivere][w].write(htmlfile,format='html',htmldict={'cssfiles':['inuso.css']})
 
      of=open(htmlfile,'r')
      codice=of.readlines()
 
      for riga in arange(len(codice)):
        codice[riga]=codice[riga].replace('*(*','<')
        codice[riga]=codice[riga].replace('*)*','>')
      nf=open(htmlfile,'w')
      nf.writelines(codice)
      nf.close()
 
      fig.xlabel('Gal. Latitude')
      fig.ylabel('Gal. Longitude')
      #fig.title(region)
      fig.axis([mw[region][0],mw[region][1],-3.7,1.5])
      fig.savefig(region.replace(' ','')+'.png')
    fig.close()

class background():
  def __init__(self, model_type='IRFs'):
    self.model = model_type
    return

  def ToXML(self, src_lib):
    if self.model == 'IRFs':
      # build an XML tree for the background
      src = ET.SubElement(src_lib, 'source', \
          {'name': 'Background', 'type': 'CTAIrfBackground', 'instrument': 'CTA'})
      # spectral sub-section
      src_spec = ET.SubElement(src, 'spectrum', {'type': 'Constant'})
      ET.SubElement(src_spec, 'parameter', \
          {'name' : 'Normalization', 'free' : '0', \
          'scale' : '1.0', 'value' : '1.0', \
          'min' : '0.0', 'max' : '1000.0'})
    else:
      print('Warning: Background model not supported: %s' % self.model)
    return

class source():
  def __init__(self, file_name='TeVJ1616-508.fits', ra=-361., dec=0, name=''):
    self.exts = pf.open(file_name)
    header = self.exts[0].header
    # parse for source name
    if name:
      self.name = name
    elif 'SOURCE' in header:
      self.name = header['SOURCE']
    else :
      self.name = 'GenericSource'
    # parse for source coordinates
    if ra != -361.:
      self.skydir = apc.SkyCoord(ra, dec, unit='deg', frame='icrs')
    elif 'L' in header and 'B' in header:
      self.skydir = apc.SkyCoord( \
          header['L'], header['B'], unit='deg', frame='galactic')
    elif 'RA' in header and 'DEC' in header:
      self.skydir = apc.SkyCoord( \
          header['RA'], header['DEC'], unit='deg', frame='icrs')
    else:
      print ('Error: missing source coordinates')
      return

  def plot(self, vediamo=0, outfile=0, over=0, point='o',xlim=[1e10,1e15],ylim=[1e-15,1e-10],nsig=1.):
    if vediamo == 0:
      vediamo=arange(len(self.exts)-1)+1
    #else :
      #print vediamo
      #vediamo=arange(len(exts)-1)+1
    #print self.name


    col = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', \
           'blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', \
           'blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black'  ]

    for ii in vediamo:
      i=int(ii)
      ext=self.exts[i]
      
      head = ext.header
          
    
      if outfile == 0 :
        print(i, ext.name, end = '')
        
        if 'PAPER' in head:
            print(', ',head['PAPER'], end = '')
        else:
            print('', end = ' ')
            
        if 'COMMENT' in head:
            print(', ',head['COMMENT'])
        else:
            print('')
          
          
      if not ext.is_image  :
      #print ext.name

        if ext.name[0:5] == 'MODEL' :
          fig.plot(ext.data['Energy'],ext.data['Flux'],'-',color=col[i],label=ext.name)
        else :
          w=where(ext.data['Flux'] != 0.)
          fig.errorbar(ext.data['Energy'][w], ext.data['Flux'][w],
                       yerr=[ext.data['FLUX_ERROR_MIN'][w]*nsig,ext.data['FLUX_ERROR_MAX'][w]*nsig],
                       fmt=point,color=col[i],label=ext.name)
          wul=where(ext.data['Flux'] == 0.)
          fig.errorbar(ext.data['Energy'][wul], ext.data['FLUX_ERROR_MAX'][wul],
                       ext.data['FLUX_ERROR_MAX'][wul]*.2,
                       uplims=True,fmt=point,color=col[i])

    if over != 0 : 
      return
    
    fig.ylim(ylim)
    fig.xlim(xlim)

    #matplotlib.rcParams.update({'font.size': 22})

    fig.xlabel('Energy  [ $eV$ ]')
    fig.ylabel('Flux  [ $erg / cm^2 s$ ]')

    fig.title(self.name)
    fig.legend(loc='best')
    fig.loglog()

    if outfile == 0 :
      fig.show()
    else :
      fig.savefig(outfile)
      fig.close()

  def ListModels(self):
    for ext in self.exts[1:]:
      
      #print('model: %s' % ext.name)
      print(ext.name, end = '')     
      head = ext.header
      
      if 'PAPER' in head:
          print(', ',head['PAPER'], end = '')
      else:
          print('', end = ' ')
           
      if 'COMMENT' in head:
          print(', ',head['COMMENT'])
      else:
          print('')
 
      
      

  def ToXML(self, src_lib, ext_model = 'MODEL'):
    spec_file_name = '%s.txt' % self.name.replace(' ','')
    # check that the model extension actually exists
    if ext_model not in self.exts:
      print('Error: missing extension %s' % ext_model)
      return False

    ext = self.exts[ext_model]
    energy = ext.data['Energy']  /1e6                    #  In MeV
    energy , ui = unique(energy, return_index=True) 
    
    flux = ext.data['Flux'][ui] /1.602e-6   / energy**2      #  In ph / cm2 sec MeV
    
    try:
      positive = np.where(flux>0)
      np.savetxt(fname = spec_file_name, \
          X = list(zip(energy[positive], flux[positive])), fmt = '%.6g %.6g')
    except IOError:
      print('Error: cannot write into file %s' % spec_file_name)
      return False

    # build an XML tree for the source
    src = ET.SubElement(src_lib, 'source', \
        {'name': self.name, 'type': 'PointSource'})
    # spectral sub-section
    src_spec = ET.SubElement(src, 'spectrum', \
        {'type': 'FileFunction', 'file': spec_file_name})
    ET.SubElement(src_spec, 'parameter', \
        {'name' : 'Normalization', 'free' : '1', \
        'scale' : '1.0', 'value' : '1.0', \
        'min' : '0.0', 'max' : '1000.0'})
    # spatial sub-section
    src_coords = ET.SubElement(src, 'spatialModel', {'type': 'SkyDirFunction'})
    ET.SubElement(src_coords, 'parameter', \
        {'name' : 'RA', 'free' : '0', \
        'scale' : '1.0', 'value' : '%.6f' % self.skydir.icrs.ra.value, \
        'min' : '0.0', 'max' : '360.0'})
    ET.SubElement(src_coords, 'parameter', \
        {'name' : 'DEC', 'free' : '0', \
        'scale' : '1.0', 'value' : '%.6f' % self.skydir.icrs.dec.value, \
        'min' : '-90.0', 'max' : '90.0'})
    return True

  def ctools_model(self, ext_model = 'MODEL', model_file_name='model.xml'):
    spec_file_name = '%s.txt' % self.name.replace(' ','')
    # check that the model extension actually exists
    if ext_model not in self.exts:
      print('Error: missing extension %s' % ext_model)
      return

    ext = self.exts[ext_model]
    energy = ext.data['Energy']/1e6                    #  In MeV
    flux = ext.data['Flux'] /1.602e-6   / energy**2      #  In ph / cm2 sec MeV
    try:
      np.savetxt(fname = spec_file_name, \
          X = list(zip(energy, flux)), fmt = '%.6g %.6g')
    except IOError:
      print('Error: cannot write into file %s' % spec_file_name)
      return
    # build an XML tree for the source
    src_lib = ET.Element('source_library', {'title': 'source library'})
    src = ET.SubElement(src_lib, 'source', \
        {'name': self.name, 'type': 'PointSource'})
    # spectral sub-section
    src_spec = ET.SubElement(src, 'spectrum', \
        {'type': 'FileFunction', 'file': spec_file_name})
    ET.SubElement(src_spec, 'parameter', \
        {'name' : 'Normalization', 'free' : '1', \
        'scale' : '1.0', 'value' : '1.0', \
        'min' : '0.0', 'max' : '1000.0'})
    # spatial sub-section
    src_coords = ET.SubElement(src, 'spatialModel', {'type': 'SkyDirFunction'})
    ET.SubElement(src_coords, 'parameter', \
        {'name' : 'RA', 'free' : '0', \
        'scale' : '1.0', 'value' : '%.6f' % self.skydir.icrs.ra.value, \
        'min' : '0.0', 'max' : '360.0'})
    ET.SubElement(src_coords, 'parameter', \
        {'name' : 'DEC', 'free' : '0', \
        'scale' : '1.0', 'value' : '%.6f' % self.skydir.icrs.dec.value, \
        'min' : '-90.0', 'max' : '90.0'})
    indent(src_lib)
    #print ET.tostring(src_lib, method='html')
    ET.ElementTree(src_lib).write(model_file_name)



def mkTable(en,flux,err=0,err2=0,name='DATA',paper='',comment=''):
    print('*mkTable* Energies must be in eV, Fluxes in erg / cm2 sec')
    sp=Table()
    sp['ENERGY']=en
    sp['ENERGY'].unit='eV'
    sp['FLUX']=flux
    sp['FLUX'].unit='erg cm-2 s-1'

    if err is not 0 :
      sp['FLUX_ERROR_MIN']=err
      sp['FLUX_ERROR_MIN'].unit='erg cm-2 s-1'

      sp['FLUX_ERROR_MAX']=err.copy()
      sp['FLUX_ERROR_MAX'].unit='erg cm-2 s-1'
       
    if err2 is not 0 :
      sp['FLUX_ERROR_MAX']=err2
      sp['FLUX_ERROR_MAX'].unit='erg cm-2 s-1'
    
    sp.meta['EXTNAME']=name
    sp.meta['PAPER']=paper
    sp.meta['COMMENT']=comment
    
    return sp


def ctools2gsed(specfile='spectrum.fits',name='CTools',comment=''):
  sp= Table.read(specfile) 
  en=sp['Energy'].data*1e12
  #en.convert_unit_to('eV')
  flux=sp['Flux'].data 
  err=sp['e_Flux'].data 
  err2=sp['e_Flux'].data.copy()
  w=where(sp['TS'] < 2)[0]
  err2[w]=sp['UpperLimit'][w].data
  flux[w]=0.
  err[w]=0.
  gt=mkTable(en,flux,err=err,err2=err2,name=name,comment=comment)
  return gt


def write(gs,file='spec.gsed'):
  d=pf.table_to_hdu(gs)
  pf.append(file,d.data,d.header)
  return




def indent(elem, level=0):
  i = "\n" + level * "  "
  if len(elem):
    if not elem.text or not elem.text.strip():
      elem.text = i + "  "
    if not elem.tail or not elem.tail.strip():
      elem.tail = i
    for elem in elem:
      indent(elem, level+1)
    if not elem.tail or not elem.tail.strip():
      elem.tail = i
  else:
    if level and (not elem.tail or not elem.tail.strip()):
      elem.tail = i

