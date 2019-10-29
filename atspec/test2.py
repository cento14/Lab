
import numpy 
from matplotlib.pyplot import *
from astropy.io import fits

insfile='IRFs/Astri_50h_June/irf_file.fits'
irf=fits.open(insfile)
  
aeff2d=irf[1].data
lene_c2= (numpy.log10(aeff2d.field('ENERG_LO')[0]) + numpy.log10(aeff2d.field('ENERG_HI')[0]) ) *.5
aeff_c = aeff2d.field('EFFAREA')[0][0]

xlabel('Energy [TeV]') ; ylabel('Eff. Area [m2]')
plot(10**lene_c2, aeff_c,'-o')
loglog()
show()  

psf2d=irf[2].data
lene_c=  (numpy.log10(psf2d.field('ENERG_LO')[0]) + numpy.log10(psf2d.field('ENERG_LO')[0]) ) *.5
r80_c=psf2d.field('SIGMA_1')[0][0]*1.8  
r68_c=psf2d.field('SIGMA_1')[0][0]*1.51

plot(10**lene_c2, r68_c,'-o')
xlabel('Energy [TeV]') ; ylabel('68% c.r. [deg]')
loglog()
show()
