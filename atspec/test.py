
from numpy import *
from matplotlib.pyplot import *
import gsed
import atspec_dev as at
from astropy import units as u

evtfile='EventFiles/astrevents.fits'
evtfile='EventFiles/astrevents_amenomori.fits'

insfile='IRFs/Astri_50h_June/irf_file.fits'
outfile='out/pippo.fits'

l = 184.557448
b = -5.784362

so=at.pointsource(evtfile=evtfile, insfile=insfile, l=l, b=b, name='Crab')    
sp=so.getspec(outfile=outfile, verbose=True,emax=150., emin=3., nch=7.)

ss=gsed.source(outfile)


# Model 1

k = 5.7e-17  *(u.cm**2 *u.s *u.MeV)**-1  # ph / cm2 s MeV 
index = -2.48 
e_piv = 1e6 * u.MeV

en = 10**( arange(110)/20. )* 1e3 *u.MeV     # MeV           
f = k * (en / e_piv)**(index)    
vfv = f * en**2.
label='Ctools (PL)'

#plot(en.to_value('eV'), vfv.to_value('erg / (cm2 s)'),label=label )          

# Model 2

k = 1.49e-15  *(u.cm**2 *u.s *u.TeV)**-1  # ph / cm2 s TeV 
index = -2.91 
e_piv = 40. * u.TeV

en = 10**( arange(110)/20. )* 1e-3 *u.TeV     # TeV           
f = k * (en / e_piv)**(index)    
vfv = f * en**2.
label='Amenomori et al. (PL)'

plot(en.to_value('eV'), vfv.to_value('erg / (cm2 s)'),label=label )          

# Model Hegra

k = 2.83e-11  *(u.cm**2 *u.s *u.TeV)**-1  # ph / cm2 s TeV 
index = -2.62 
e_piv = 1. * u.TeV

en = 10**( arange(110)/20. )* 1e-3 *u.TeV     # TeV           
f = k * (en / e_piv)**(index)    
vfv = f * en**2.
label='HEGRA (PL)'

plot(en.to_value('eV'), vfv.to_value('erg / (cm2 s)'),label=label )          




# Model Hawc

k = 2.35e-13  *(u.cm**2 *u.s*u.TeV)**-1  # typo nel paper ! -->  ph / TeV cm2 s  
alpha = 2.79
beta = 0.06
e_piv = 7. *u.TeV

en = 10**( arange(110)/20. )* 1e-3 *u.TeV     # TeV           
f = k * (en / e_piv)**(-alpha-beta*log(en/e_piv))    
vfv = f * en**2.
label='Hawc spectrum (Log Par.)'






# Plot

plot(en.to_value('eV'), vfv.to_value('erg / (cm2 s)'),label=label )          # Model
ss.plot(nsig=2,xlim=[3e12,210.e12],ylim=[1e-13,1e-9],point='o')





