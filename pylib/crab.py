

#import sys, os
from numpy import *
#from matplotlib.pyplot import *
#from astropy.table import Table 
from astropy import units as u
#from astropy.io import fits



def pl(en,k=1.,index=2,e_piv=1.):
  
  flux = k * (en / e_piv)**(index) 
  return flux



def crab(energy=[.1,10],model='hegra', giveme='intflux'):

  try : 
    print('Input energy units :',energy.unit) 
  except: 
    energy=energy*u.TeV 
    print('Energy units set to Tev') 

  if giveme == 'intFlux' :

    emin=energy[0].to('TeV')
    emax=energy[1].to('TeV')

    dd=log10(emax/emin)

    energy_b = 10**( arange(dd*1000.+1.)/1000. )* emin      # TeV           
    de = energy_b[1:]-energy_b[:-1]

    en = sqrt(energy_b[1:]*energy_b[:-1])

  else :
    
    en=energy    
    

  def hegra(en):   # HEGRA

    k = 2.83e-11  *(u.cm**2 *u.s *u.TeV)**-1  # ph / cm2 s TeV 
    index = -2.62 
    e_piv = 1. * u.TeV
    flux = pl(en,k=k, index=index, e_piv=e_piv)
  
    return flux


  def amenomori(en):      # Amenomori et al. (PL)

    k = 1.49e-15  *(u.cm**2 *u.s *u.TeV)**-1  # ph / cm2 s TeV 
    index = -2.91 
    e_piv = 40. * u.TeV
    flux = pl(en,k=k, index=index, e_piv=e_piv)

    return flux


  def hawc(en):          # Hawc model (Log.Par)

    k = 2.35e-13  *(u.cm**2 *u.s*u.TeV)**-1  # typo nel paper ! -->  ph / TeV cm2 s  
    alpha = 2.79
    beta = 0.06
    e_piv = 7. *u.TeV
    flux = k * (en / e_piv)**(-alpha-beta*log(en/e_piv))   

    return flux


  crabFlux={'hegra':hegra, 'amenomori':amenomori, 'hacw':hawc}

  try:
    flux=crabFlux[model](en)
  except: 
    print('Available models : ',crabFlux.keys())
    flux=0.

  
  sed=flux*en*en.to('erg')

  intFlux  = sum(flux*de )
  intFluxE = sum(flux*de *en.to('erg'))

  
    

  return intFlux,intFluxE

