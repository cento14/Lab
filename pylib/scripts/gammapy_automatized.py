import matplotlib.pyplot as plt

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle
from astropy.convolution import Gaussian2DKernel
from regions import CircleSkyRegion
from gammapy.utils.energy import EnergyBounds
from gammapy.data import DataStore
from gammapy.spectrum import (
                              SpectrumExtraction,
                              SpectrumFit,
                              SpectrumResult,
                              models,
                              SpectrumEnergyGroupMaker,
                              FluxPointEstimator,
                              )
from gammapy.maps import Map, MapAxis, WcsNDMap, WcsGeom
from gammapy.cube import MapMaker
from gammapy.background import ReflectedRegionsBackgroundEstimator
from gammapy.detect import TSMapEstimator, find_peaks

import logging

logging.basicConfig()
log = logging.getLogger("gammapy.spectrum")
log.setLevel(logging.ERROR)


## SELECT EVENTS
print("#################")
print("Starting the analysis of the dataset...")
print("#################")

data_store = DataStore.from_dir("./")


ra_src=248.04  # DA SETTARE A MANO PER INDIVIDUARE LA SORGENTE PRINCIPALE DEL CAMPO
dec_src=-47.82  # DA SETTARE A MANO PER INDIVIDUARE LA SORGENTE PRINCIPALE DEL CAMPO
src_pos_galactic = (SkyCoord(ra_src*u.degree, dec_src*u.degree, frame='icrs')).galactic


#Just as a reminder: this is how to select observations
table = data_store.obs_table
pos_obs = SkyCoord(table['GLON_PNT'], table['GLAT_PNT'], frame='galactic', unit='deg')
pos_target = src_pos_galactic
offset = pos_target.separation(pos_obs).deg
mask = (-8 < offset) & (offset < 8)
table = table[mask]
#table.show_in_browser(jsviewer=True)


obs_id = table['OBS_ID']
obs_list = data_store.obs_list(obs_id)

obs_cols = ['OBS_ID', 'RA_PNT', 'DEC_PNT', 'GLON_PNT', 'GLAT_PNT', 'LIVETIME']
data_store.obs_table.select_obs_id(obs_id)[obs_cols]



## MAPS
print("#################")
print("These are the images of your event lists...")
print("#################")

axis = MapAxis.from_edges(
                          np.logspace(0.045, 1.97, 12), unit="TeV", name="energy", interp="log"
                          )
geom = WcsGeom.create(
                      skydir=(ra_src, dec_src), npix=(250, 250), binsz=0.02, coordsys="CEL", axes=[axis]
                      )    # DA SETTARE A MANO PER CREARE MAPPE DI DIMENSIONI DIVERSE
#geom



#exclusion_mask = geom.to_image().region_mask([on_region], inside=False)
#exclusion_mask = WcsNDMap(geom.to_image(), exclusion_mask)
#exclusion_mask.plot();

maker = MapMaker(geom, offset_max="5 deg")   # DA SETTARE A MANO PER CREARE MAPPE DI DIMENSIONI DIVERSE
maps = maker.run(obs_list)
print(maps.keys())


images = maker.make_images()

excess = images["counts"].copy()
excess.data -= images["background"].data
images["excess"] = excess


images["counts"].smooth(2).plot(vmax=40);
name="detected_sources_counts.eps"
plt.savefig(name, format='eps', dpi=100)
plt.gcf().clear()
plt.close

images["background"].plot(vmax=25);
name="detected_sources_background.eps"
plt.savefig(name, format='eps', dpi=100)
plt.gcf().clear()
plt.close

images["excess"].smooth(3).plot(vmax=20);
plt.show()


### SOURCE DETECTION
#print("#################")
#print("Starting the source detection...")
#print("#################")
#
#
#kernel = Gaussian2DKernel(1, mode="oversample").array
#plt.imshow(kernel);
#
#ts_image_estimator = TSMapEstimator()
#images_ts = ts_image_estimator.run(images, kernel)
#print(images_ts.keys())
#
#sources = find_peaks(images_ts["sqrt_ts"], threshold=11) # DA SETTARE A MANO PER SELEZIONARE LA SIGMA DI DETECTION
#sources
#
#source_pos = SkyCoord(sources["ra"], sources["dec"])
#source_pos
#
#
#images_ts["sqrt_ts"].plot(add_cbar=True)
#
#plt.gca().scatter(
#                  source_pos.ra.deg,
#                  source_pos.dec.deg,
#                  transform=plt.gca().get_transform("icrs"),
#                  color="none",
#                  edgecolor="white",
#                  marker="o",
#                  s=200,
#                  lw=1.5,
#                  );
#
#name="detected_sources_ts.eps"
#plt.savefig(name, format='eps', dpi=100)
#plt.gcf().clear()
#plt.close
#
#print("#################")
#print("....These are the detected sources at "+str(thr)+"sigma...")
#print("#################")

########
#Sources

src_sim_ra = np.array([248.04, 248.74, 244.091682, 244.724776, 246.753409, 250.144976, 250.165452, 250.259430, 251.979034])
src_sim_dec = np.array([-47.82, -47.27, -50.915090, -50.659715, -49.173316, -46.654457, -46.536048, -46.311721, -46.258415])


# SPECTRAL ANALYSIS
print("#################")
print("Getting the spectrum of each source...")
print("#################")

#thresh = 12 # SOGLIA IN SIGMA, DA SETTARE A MANO PER DECIDERE QUALI SORGENTI SI VOGLIONO ANALIZZARE E TAGLIARE TUTTE QUELLE SPURIE
for i in np.arange(0,len(src_sim_ra)):
#    if sources['value'][i] > thresh:
        sz = 0.4
        pos = SkyCoord(src_sim_ra[i], src_sim_dec[i], unit='deg').galactic
        pos_rest = SkyCoord(src_sim_ra, src_sim_dec, unit='deg').galactic
        radius = Angle(sz, 'deg')
        on_region = CircleSkyRegion(center=pos, radius=radius)

        sep = pos.separation(pos_rest)
        idx = np.where(sep > 0.2 *u.deg)
        other_src = pos_rest[idx]
        src_excluded = CircleSkyRegion(center=other_src, radius=radius)
        
#        src_excluded=[]
#        for j in np.delete(np.arange(0,len(src_sim_ra)),i):
#            if sep[j] > (0.2 * u.deg):
#                src2_pos_galactic = (SkyCoord(src_sim_ra[j], src_sim_dec[j], unit='deg', frame='icrs')).galactic
#                src2_radius = 0.2 * u.deg
#                src2_region = CircleSkyRegion(center=src2_pos_galactic, radius=src2_radius)
#                src_excluded.append(src2_region)


        exclusion_mask = geom.to_image().region_mask(src_excluded, inside=False)
        exclusion_mask = WcsNDMap(geom.to_image(), exclusion_mask)
        exclusion_mask.plot();


        bkg_estimator = ReflectedRegionsBackgroundEstimator(
                                                            obs_list=obs_list, on_region=on_region, exclusion_mask=exclusion_mask
                                                            )
        bkg_estimator.run()
        bkg_estimate = bkg_estimator.result
        #bkg_estimator.plot();


        extract = SpectrumExtraction(obs_list=obs_list, bkg_estimate=bkg_estimate)
        extract.run()
        observations = extract.observations



        model = models.PowerLaw(
                                   index = 2.12,
                                   amplitude = 5e-12 * u.Unit('cm-2 s-1 TeV-1'),
                                   reference = 1 * u.TeV,
                                   )
        fit = SpectrumFit(observations, model, fit_range=[1.5,95.0]*u.TeV)
        fit.run()
        print(fit.result[0])


        stacked_obs = extract.observations.stack()
        print(stacked_obs)


        ebounds = EnergyBounds.equal_log_spacing(3.0, 95, 10, unit=u.TeV)

        seg = SpectrumEnergyGroupMaker(obs=stacked_obs)
        seg.compute_groups_fixed(ebounds=ebounds)

        fpe = FluxPointEstimator(
                                 obs=stacked_obs, groups=seg.groups, model=fit.result[0].model
                                 )
        fpe.compute_points()
        fpe.flux_points.table


        #spec_simul = SpectrumResult(
        #                              model=model, points=fpe.flux_points
        #                              )

        #model.plot(energy_range=[1, 100] * u.TeV,energy_power=2);
        #plt.show()

        total_result = SpectrumResult(
                                      model=fit.result[0].model, points=fpe.flux_points
                                      )
        ax0, ax1 = total_result.plot(
                          energy_range=[1.1, 90] * u.TeV,
                          energy_power=2, flux_unit='erg-1 cm-2 s-1',
                          fig_kwargs=dict(figsize=(8, 8)),
                          point_kwargs=dict(color="green"),
                          );

        #ax0.set_xlim(1.0, 100)
        #ax0.set_ylim(1e-12, 1e-11)

        #opts = {
        #    "energy_range": [1.1, 100] * u.TeV,
        #    "energy_power": 2,
        #    "flux_unit": "erg-1 cm-2 s-1",
        #}
        #CrabSpectrum('magic_lp').model.plot(ax=ax0, **opts)

        #total_mod_sim = SpectrumResult(
        #                              model=model, points=fpe.flux_points
        #                              )
        #total_mod_sim.plot(
        #                   energy_range=[1.1, 90] * u.TeV,
        #                   energy_power=2, flux_unit='erg-1 cm-2 s-1',
        #                   fig_kwargs=dict(figsize=(8, 8)),
        #                   point_kwargs=dict(color="red"),
        #                   );

        #plt.show()

        name="spectrum_source_powerlaw"+"_"+str(i)+".eps"
        plt.savefig(name, format='eps', dpi=100)
        plt.gcf().clear()
        plt.close


        np.savetxt("gammapy_powerlaw"+"_"+str(i)+".txt", np.c_[fpe.flux_points.table['e_ref']*1e12, fpe.flux_points.table['dnde']*fpe.flux_points.table['e_ref']**2.*1.60217330, fpe.flux_points.table['dnde_err']*fpe.flux_points.table['e_ref']**2.*1.60217330, fpe.flux_points.table['dnde_err']*fpe.flux_points.table['e_ref']**2.*1.60217330, fpe.flux_points.table['dnde_ul']*fpe.flux_points.table['e_ref']**2.*1.60217330, fpe.flux_points.table['sqrt_ts']],delimiter=' ')
