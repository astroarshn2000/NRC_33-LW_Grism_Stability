""" NIRCam CAR-3: TS Stability (HD 140982 - GRISMR) Simulation
             ~ Arsh R. Nadkarni (UArizona), 2021 """  

# ----------   Set Environment Variables   ----------

import os
os.environ["MIRAGE_DATA"] = "/home/anadkarni/JWST_Programs/mirage_reference_files/mirage_data/"
os.environ["CRDS_PATH"] = os.path.expandvars("$HOME/crds_cache")
os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"
os.environ["PYSYN_CDBS"] = "/home/anadkarni/JWST_Programs/stsynphot_reference_files/grp/hst/cdbs/"
os.environ["WEBBPSF_PATH"] = "/home/anadkarni/JWST_Programs/webbpsf_reference_files/webbpsf-data/"


#  ----------   Import Libraries  ----------

from astropy.io import fits, ascii
from astropy.table import Table
from astropy.visualization import simple_norm, imshow_norm
from astropy import units as u
import batman
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams; rcParams["figure.dpi"] = 300
from matplotlib.ticker import (AutoMinorLocator)
from matplotlib import colors
import matplotlib.cm as cmx
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-small')
import stsynphot as stsyn
from synphot import SourceSpectrum, SpectralElement
from synphot import units
import yaml
from mirage.catalogs.hdf5_catalog import save_tso
from mirage.catalogs.catalog_generator import GrismTSOCatalog, ImagingTSOCatalog, PointSourceCatalog
from mirage.catalogs.catalog_generator import TSO_GRISM_INDEX
from mirage.grism_tso_simulator import GrismTSO
from mirage.imaging_simulator import ImgSim
from mirage.seed_image.catalog_seed_image import Catalog_seed
from mirage.utils.utils import ensure_dir_exists
from mirage.yaml import yaml_generator


# ----------   Define Paths to organize Inputs and Outputs  ----------

input_data_path = os.path.abspath('/home/anadkarni/NRC_33-LW_Grism_Stability/NCam033_pid1442_grism_stability_apt_data/')
output_dir = '/home/anadkarni/NRC_33-LW_Grism_Stability/'
output_yaml_dir = os.path.abspath('/home/anadkarni/NRC_33-LW_Grism_Stability/NCam033_pid1442_grism_stability_yaml_data/')
ensure_dir_exists(output_yaml_dir)
output_data_dir = os.path.abspath('/home/anadkarni/NRC_33-LW_Grism_Stability/NCam033_pid1442_grism_stability_sim_data/')
ensure_dir_exists(output_data_dir)
tsdir = '/home/anadkarni/NRC_33-LW_Grism_Stability/HD140982_Params/'

# ----------  Prepare Inputs  ----------

xml_file = os.path.join(input_data_path, 'NCam033_pid1442_grism_stability.xml')
pointing_file = xml_file.replace('.xml', '.pointing')


# ----------  Stellar Spectrum  ----------

t_eff = 6400  # surface temperature
metallicity = 0 # Fe/H 
log_g = 4.3  # log(surface gravity)
sp = stsyn.grid_to_spec('ck04models', t_eff, metallicity, log_g)
bp = SpectralElement.from_filter('johnson_k')
vega = SourceSpectrum.from_vega()
sp_norm = sp.normalize(7.60 * units.VEGAMAG, bp, vegaspec=vega)
wavelengths = sp_norm.waveset.to(u.micron)
fluxes = sp_norm(wavelengths, flux_unit='flam')
wavelength_units = 'microns'
flux_units = 'flam'
sed_file = os.path.join(output_dir, 'NCam033_pid1442_grism_stability_stellar_spectrum.hdf5')
fluxes = [fluxes]
wavelengths = [wavelengths]
with h5py.File(sed_file, "w") as file_obj:
    for i in range(len(fluxes)):
        dset = file_obj.create_dataset(str(i+TSO_GRISM_INDEX), data=[wavelengths[i].value, fluxes[i].value],
                                       dtype='f', compression="gzip", compression_opts=9)
        dset.attrs[u'wavelength_units'] = wavelength_units
        dset.attrs[u'flux_units'] = flux_units


# ----------  BATMAN Parameters and Light-Curve File  ----------

params = batman.TransitParams()       # object to store transit parameters
params.t0 = 2458651.789013            # epoch of eclipse center
params.per = 1.2696257                # orbital period
params.rp = 0.12845                   # planet/EB radius (in units of stellar radii)
params.a = 6.5                        # semi-major axis (in units of stellar radii)
params.inc = 90                       # orbital inclination (in degrees)
params.ecc = 0                        # eccentricity
params.w = 90                         # longitude of periastron (in degrees)
params.fp = 0                         # planet-to-star flux ratio
params.t_secondary = 0                # central eclipse time
params.limb_dark = "quadratic"        # limb darkening model
params.u = [0.038056000, 0.25594800]  # limb darkening coefficients (ExoFAST)
nint = 1336
times = np.linspace(2458651.688213, 2458651.8898130003, nint) # Total = 4.84 hrs. So, 2.42 hrs = 0.1008 days on either side of t0
m = batman.TransitModel(params, times)
flux = m.light_curve(params)
lightcurve_file = os.path.join(output_dir, 'HD-140982_lightcurve.hdf5')
contents = {}
contents['1'] = {'times': times,
                 'fluxes': flux}
save_tso(contents, lightcurve_file, time_unit='days')

# ----------  Plot Light Curve  ----------

f, a = plt.subplots()
a.plot(times, flux, ls='None', marker='.', ms=0.5, zorder=1, color='black')
a.set_xlabel('Time (BJD)')
a.set_ylabel('Normalized Flux')
a.set_title('HD-140982 Eclipse', fontsize='large', fontweight='bold')
a.xaxis.set_minor_locator(AutoMinorLocator()) 
a.yaxis.set_minor_locator(AutoMinorLocator())
a.tick_params(which='minor', length=2.5, color='k')
plt.savefig('HD-140982_LC.png')
plt.show()


# ----------  Transmission Spectrum  ----------

tran_spec_file = os.path.join(tsdir,'HD140982_tran_spec.txt')


# ---------- Create Grism TSO Catalog   ----------

grism_tso_catalog = os.path.join(output_dir,'tso_grism_source.cat')
object_ra = +235.8768367
object_dec = +59.11375
object_f322w2_mag = 7.588659701753237
grism_cat = GrismTSOCatalog(ra=[object_ra], dec=[object_dec], semimajor_axis=[params.a],
                            orbital_inclination=[params.inc], eccentricity=[params.ecc],
                            orbital_period=[params.per], longitude_of_periastron=[params.w],
                            limb_dark_model=[params.limb_dark], limb_dark_coeffs=[params.u],
                            time_units=['days'], start_time=[np.min(times)],
                            end_time=[np.max(times)], inferior_conj=[params.t0],
                            transmission_spectrum=[tran_spec_file])
grism_cat.add_magnitude_column([object_f322w2_mag], magnitude_system='vegamag',
                               instrument='nircam', filter_name='f322w2')
grism_cat.save(grism_tso_catalog)


# ---------- Create Imaging TSO Catalog   ----------

object_f210m_mag = 7.600444438235966
imaging_tso_catalog = os.path.join(output_dir, 'tso_imaging_source.cat')
tsimg_cat = ImagingTSOCatalog(ra=[object_ra], dec=[object_dec], lightcurve_file=[lightcurve_file])
tsimg_cat.add_magnitude_column([object_f210m_mag], magnitude_system='vegamag',
                               instrument='nircam', filter_name='wlp8')
tsimg_cat.save(imaging_tso_catalog)


# ---------- Create Input YAML Files for Mirage   ----------

catalogs = {'HD-140982': {'point_source': None,
                         'tso_imaging_catalog': imaging_tso_catalog,
                         'tso_grism_catalog': grism_tso_catalog,
                        }
           }
background = 'medium'
pav3 = 0
yam = yaml_generator.SimInput(xml_file, pointing_file, catalogs=catalogs, verbose=True,
                              output_dir=output_yaml_dir, simdata_output_dir=output_data_dir,
                              background=background, roll_angle=pav3,
                              dates='2022-3-1', datatype='linear, raw', dateobs_for_background=True,
                              reffile_defaults='crds')

yam.use_linearized_darks = True
yam.create_inputs()

# ---------- SIMULATE F322W2 GRISMR TSO ---------- 

gr_tso_yaml_file = os.path.join(output_yaml_dir, 'jw01442001001_01101_00002_nrca5.yaml')
gr_f322w2 = GrismTSO(gr_tso_yaml_file, SED_file=sed_file, SED_normalizing_catalog_column=None,
                    final_SED_file=None, save_dispersed_seed=True, source_stamps_file=None,
                    extrapolate_SED=True, override_dark=None, disp_seed_filename=None,
                    orders=["+1", "+2"])
gr_f322w2.create()

# ---------- SIMULATE WLP8+F210M IMAGING TSO ---------- 

img_tso_sw_A1_yaml = os.path.join(output_yaml_dir, 'jw01442001001_01101_00001_nrca1.yaml')
img_tso_sw_A3_yaml = os.path.join(output_yaml_dir, 'jw01442001001_01101_00001_nrca3.yaml')
img_tso_A1 = ImgSim(paramfile=img_tso_sw_A1_yaml)
img_tso_A1.create()
img_tso_A3 = ImgSim(paramfile=img_tso_sw_A3_yaml)
img_tso_A3.create()