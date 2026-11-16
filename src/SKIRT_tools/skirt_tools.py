import h5py
import shutil
import numpy as np
from lxml import etree
from pathlib import Path
import astropy.units as u
import sys
from astropy.io import fits
import matplotlib.pyplot as plt

import pts.band
import pts.do
import pts.simulation as sm
from pts.storedtable.io import writeStoredColumns

#sys.path.append("/export/nfs0home/jennydw/FIRE_studio2/src") # Replace with where your firestudio folder is installed

from firestudio.utils.stellar_utils import load_stellar_hsml, make_threeband_image
from . import observation_utils

if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid

SKIRT_tools_dir = Path(__file__).resolve().parent

class Galaxy:
    """
    A representation of a galaxy containing the data required to run a SKIRT simulation.

    Attributes
    ----------
    name : str
        The prefix that is used to name the run directory, SKIRT parameter file, and output files. Ex: object_50_s600
    star_cols : dict
        A dictionary containing the all star particle data required to run the SKIRT simulation, in the form
        
        star_cols = {
            "coord"        : (N, 3) numpy array, units: kpc
            "initial_mass" : (N,) numpy array, units: Msun
            "metallicity"  : (N,) numpy array, units: 1
            "age"          : (N,) numpy array, units: Gyr
            "velocity"     : (N, 3) numpy array, units: km/s
                (OPTIONAL -- Includes Doppler shift from velocity relative to the model coordinate system)
        }
    gas_cols : dict
        A dictionary containing the all gas particle data required to run the SKIRT simulation, in the form  
        
        gas_cols = {
            "coord"        : (N, 3) numpy array, units: kpc
            "mass"         : (N,) numpy array, units: Msun
            "metallicity"  : numpy array, units: 1
            "hsml"         : (N,) numpy array, units: kpc
            "temperature"  : (N,) numpy array, units: K
                (OPTIONAL -- Provides the option to ignore gas particles above a maximum temperature value)
            "velocity"     : (N, 3) numpy array, units: km/s
                (OPTIONAL -- Includes Doppler shift from velocity relative to the model coordinate system)
        }        
    """
    
    def __init__(self, name, star_cols=None, gas_cols=None, hdf5_file=None, snapshot=None):
        self.name = name
        self.hdf5_file = hdf5_file
        self.snapshot = snapshot      

        # loading from hdf5 file
        if self.hdf5_file:
            self.star_cols, self.gas_cols = {}, {}
            try: 
                with h5py.File(self.hdf5_file, "r") as f:
                    for star_key in f['star'].keys():
                        self.star_cols[star_key] = f['star'][star_key][()]
                    for gas_key in f['gas'].keys():
                        self.gas_cols[gas_key] = f['gas'][gas_key][()]
                try: snapshot = f.attrs["snapshot"]
                except Exception: print("No snapshot number found.")
            except Exception as e:
                print(f"Error occured while loading hdf5 file:{e}.")
                if (self.star_cols and self.gas_cols):
                    print("Using star columns and gas columns instead.")
                    self.star_cols = star_cols
                    self.gas_cols = gas_cols
                else: raise ValueError("Data columns are empty.")
        else:
            # loading from column data
            if (star_cols and gas_cols):
                print("Loading data columns...")
                self.star_cols = star_cols
                self.gas_cols = gas_cols
            # empty data
            else: raise ValueError("Must provide column data or hdf5 file.")
            

    def save_to_hdf5(self, save_folder):
        with h5py.File(self.name + ".h5", "w") as f:
            star = f.create_group("star")
            gas = f.create_group("gas")
            for star_key, data in self.star_cols.items():
                star.create_dataset(star_key, data=data)
            for gas_key, data in self.gas_cols.items():
                gas.create_dataset(gas_key, data=data)
            if self.snapshot:
                f.attrs["snapshot"] = self.snapshot

class Sim:
    """
    Child class that inherits attributes of the Galaxy object. Writes input files, parameter files, 
    and runs the SKIRT simulation.

    Attributes
    ----------
    rundir : str
        Directory that stores simulation input and output files
    max_radius : int or float
        The maximum radial distance (Kpc) of a star or gas particle. Particles outside of this radius are not used in the simualtion.
        If max_radius is None, all particles are kept.
        Default: None
    import_velocity : bool
        Determines whether to include Doppler shift from velocity relative to model coordinates
        Default: False
        (NOTE: requires a "velocity" column in gas_cols and star_cols)
    import_gas_temp : bool
        Determines whether to include gas temperature input files, which enables the option to set a maximum gas temperature in the simulation.
        Default: False
        (NOTE: requires a "temperature" column in gas_cols)
    n_random_projections : int
        The number of random projections generated by the simulation, in addition to 
        the default xy, xz, and yz projections. For each random projection, SKIRT will output a fits_file and sed_file.
    n_photon_packets : int
        The number of photon packets launched in the simulation. Increasing the number photon packets makes the simulation more computationally expensive.

    Methods
    ----------
    Pre-simulation
        add_instrument : Add an instrument to the SKIRT parameter file 
        initialize : Write input files and parameter file to <rundir>
        run : Runs the simulation from <rundir>

    Post-simulation
        get_band_flux   : Return dictionary of band fluxes
        get_luminosity  : Return luminosity
        get_band_images : Return dictionary of band images
        get_mock_image  : Return three-band mock image
    """

    def __init__(
        self, 
        Galaxy,
        rundir,
        max_radius=None,
        import_velocity=False, 
        import_gas_temp=True, 
        n_random_projections=0, 
        n_photon_packets=1e7):

        self.galaxy = Galaxy
        self.rundir = rundir
        self.max_radius = max_radius
        self.import_velocity = import_velocity
        self.import_gas_temp = import_gas_temp
        self.n_random_projections = n_random_projections
        self.n_photon_packets = n_photon_packets
        self.n_star_particles_init = len(self.galaxy.star_cols['coord'][:,0])
        self.n_gas_particles_init = len(self.galaxy.gas_cols['coord'][:,0])

    def _build_scol_files(self, star=True, gas=True):
        if star == False and gas == False:
            return
        star_colnames = ['x', 'y', 'z', 'hsml', 'mass', 'M', 'age', 'vx', 'vy', 'vz']
        gas_colnames = ['x', 'y', 'z', 'hsml', 'mass', 'M', 'temp', 'vx', 'vy', 'vz']
        star_units = ['kpc','kpc','kpc','kpc','Msun','1','Gyr','km/s','km/s','km/s']
        gas_units = ['kpc','kpc','kpc','kpc','Msun','1','K','km/s','km/s','km/s']
        if not self.import_velocity:
            star_colnames = star_colnames[:7]
            gas_colnames = gas_colnames[:7]
            star_units = star_units[:7]
            gas_units = gas_units[:7]
        if not self.import_gas_temp:
            gas_colnames.remove('temp')
            gas_units.remove('K')

        if self.max_radius:
            if star: print("Number of star particles before cut: " + str(self.n_star_particles_init))
            if gas: print("Number of gas particles before cut: " + str(self.n_gas_particles_init))
        else:
            if star: print("Number of star particles: " + str(self.n_star_particles_init))
            if gas: print("Number of gas particles: " + str(self.n_gas_particles_init))

        def get_mask(r_max):
            star_coords = self.galaxy.star_cols['coord']
            gas_coords = self.galaxy.gas_cols['coord']
            if r_max is not None and r_max > 0:
                r_max_sq = r_max**2
                mask_star = np.sum(star_coords**2, axis=1) <= r_max_sq
                mask_gas = np.sum(gas_coords**2, axis=1) <= r_max_sq
            else:
                mask_star = np.ones(len(star_coords), dtype=bool)
                mask_gas = np.ones(len(gas_coords), dtype=bool)
            return mask_star, mask_gas

        mask_star, mask_gas = get_mask(self.max_radius)
        if self.max_radius:
            print("Applying " + str(self.max_radius) + " kpc cut...")
            if star:
                n_star = len(self.galaxy.star_cols['coord'][:,0][mask_star])
                print("Number of star particles after cut: " + str(n_star))
            if gas:
                n_gas = len(self.galaxy.gas_cols['coord'][:,0][mask_gas])
                print("Number of gas particles after cut: " + str(n_gas))

        if star:
            print("Importing star data...")
            # Calculate star particle SKIRT smoothing lengths with FIRE_studio
            star_hsml = load_stellar_hsml.get_particle_hsml(
                self.galaxy.star_cols['coord'][:,0], 
                self.galaxy.star_cols['coord'][:,1], 
                self.galaxy.star_cols['coord'][:,2])
            star_cols = [
                self.galaxy.star_cols['coord'][:,0][mask_star],
                self.galaxy.star_cols['coord'][:,1][mask_star],
                self.galaxy.star_cols['coord'][:,2][mask_star],
                star_hsml[mask_star],
                self.galaxy.star_cols['initial_mass'][mask_star],
                self.galaxy.star_cols['metallicity'][mask_star],
                self.galaxy.star_cols['age'][mask_star]]
            if self.import_velocity: star_cols.append(self.galaxy.star_cols['velocity'][mask_star])
            writeStoredColumns(self.rundir + '/star.scol', star_colnames, star_units, star_cols)

        if gas:
            print("Importing gas data...")
            gas_cols = [
                self.galaxy.gas_cols['coord'][:,0][mask_gas],
                self.galaxy.gas_cols['coord'][:,1][mask_gas],
                self.galaxy.gas_cols['coord'][:,2][mask_gas],
                self.galaxy.gas_cols['hsml'][mask_gas],
                self.galaxy.gas_cols['mass'][mask_gas],
                self.galaxy.gas_cols['metallicity'][mask_gas]]
            if self.import_gas_temp: gas_cols.append(self.galaxy.gas_cols['temperature'][mask_gas])
            if self.import_velocity: gas_cols.append(self.galaxy.gas_cols['velocity'][mask_gas])
            writeStoredColumns(self.rundir + '/gas.scol', gas_colnames, gas_units, gas_cols)

    def add_instrument(
        self,
        instrumentName, 
        inclination, azimuth, roll,
        instrumentType='FullInstrument',
        distance=50,
        FOV_x=60, npixel_x=400, center_x=0,
        FOV_y=60, npixel_y=400, center_y=0,
        recordComponents=True, numScatteringLevels=0, 
        recordPolarization=True, recordStatistics=True):

        """
        Add an instrument to the parameter file with a specified viewing angle, field of view, 
        and additional parameters.
        See SKIRT documentation (https://skirt.ugent.be/skirt9/class_instrument.html) for more details.
        
        Args
        ----------
        instrumentName : str, Prefix for the output files generated by the instrument
        inclination    : int or float between -360 and 360 deg
        azimuth        : int or float between    0 and 180 deg
        roll           : int or float between -360 and 360 deg
        
        instrumentType : "FullIntrument" or "FrameInstrument" or "SEDIntrument" (default: "FullInstrument")
            "FullInstrument"  -- Outputs both FITS files and SED files
            "FrameInstrument" -- Only outputs FITS files
            "SEDIntrument"    -- Only outputs SED files
            NOTE: SED files are required for band flux calculations, and FITS files are required to make three-band images.
        
        distance : int or float, Distance to galaxy in Mpc (default: 50)

        The following are only used if instrumentType is "FullInstrument" or "FrameInstrument":
            FOV_x, FOV_y, FOV_z          : int or float, Full range in Kpc in the x, y, z direction captured by the instrument
            npixel_x, npixel_y, npixel_z : int, Dimensions (pixels) of FITS file images
            center_x, center_y, center_z : int or float, The coordinates in Kpc that specify the center of the image
        
        recordComponents : bool
            True -- outputs images with dust, no dust, and primary scattering only 
            False -- only outputs images with dust 
        """

        ski_file = str(Path(self.rundir, self.galaxy.name + '.ski'))
        tree = etree.parse(ski_file)
        root = tree.getroot()
    
        instruments = root.find(".//instruments")
        instrument = etree.SubElement(instruments, instrumentType)
    
        instrument.set("instrumentName", instrumentName)
        instrument.set("distance", f"{distance} Mpc")
        instrument.set("inclination", f"{inclination} deg")
        instrument.set("azimuth", f"{azimuth} deg")
        instrument.set("roll", f"{roll} deg")
        if instrumentType == "FullInstrument" or instrumentType == "FrameInstrument":
            instrument.set("fieldOfViewX", f"{FOV_x} kpc")
            instrument.set("numPixelsX", f"{npixel_x}")
            instrument.set("centerX", f"{center_x} kpc")
            instrument.set("fieldOfViewY", f"{FOV_y} kpc")
            instrument.set("numPixelsY", f"{npixel_y}")
            instrument.set("centerY", f"{center_y} kpc")
        instrument.set("recordComponents", f"{str(recordComponents).lower()}")
        instrument.set("numScatteringLevels", f"{numScatteringLevels}")
        instrument.set("recordPolarization", f"{str(recordPolarization).lower()}")
        instrument.set("recordStatistics", f"{str(recordStatistics).lower()}")
    
        tree.write(ski_file, encoding="UTF-8", xml_declaration=True, pretty_print=True)

    def _build_ski_file(self):
        ski_default = str(SKIRT_tools_dir / "galaxy.ski")
        ski_path = str(Path(self.rundir, self.galaxy.name + '.ski'))
        shutil.copy(ski_default, ski_path)
        pts.do.initializePTS()
        skifile = sm.SkiFile(ski_path)

        xpath_source = 'MonteCarloSimulation/sourceSystem/SourceSystem/sources/ParticleSource'
        xpath_medium = 'MonteCarloSimulation/mediumSystem/MediumSystem/media/ParticleMedium'
        
        skifile.setStringAttribute(xpath_source, 'filename', str(Path(self.rundir, 'star.scol')))
        skifile.setStringAttribute(xpath_medium, 'filename', str(Path(self.rundir, 'gas.scol')))
        skifile.setStringAttribute('MonteCarloSimulation', 'numPackets', str(self.n_photon_packets))

        # Import optional columns
        if self.import_velocity:
            skifile.setStringAttribute(xpath_source, 'importVelocity', "True")
            skifile.setStringAttribute(xpath_medium, 'importVelocity', "True")
        if self.import_gas_temp:
            skifile.setStringAttribute(xpath_medium, 'importTemperature', "True")

        # Match list of column names to imported columns
        columns_str = {"star":"x,y,z,hsml,mass,M,age", "gas":"x,y,z,hsml,mass,M"}
        if self.import_gas_temp: columns_str["gas"] += ",temp"
        if self.import_velocity:
            columns_str["star"] += ",vx,vy,vz"
            columns_str["gas"]  += ",vx,vy,vz"
        skifile.setStringAttribute(xpath_source, 'useColumns', columns_str["star"])
        skifile.setStringAttribute(xpath_medium, 'useColumns', columns_str["gas"])

        # Add instrument for every random projection
        if self.n_random_projections:
            for n in range(self.n_random_projections):
                inc, azi, roll = observation_utils.random_viewing_angles()
                add_instrument(ski_path, f'rand{n+1}', inc, azi, roll)
            
        skifile.saveTo(ski_path)    

    def build(self):
        """Writes column data files and parameter files to <rundir>"""
        if not Path(self.rundir).exists():
            print(f"Run directory ({rundir}) does not exist.\nCreating run directory...")
            output_dir = Path(self.rundir)
            output_dir.mkdir(parents=True, exist_ok=True)
        star_scol_path = Path(self.rundir, 'star.scol')
        gas_scol_path = Path(self.rundir, 'gas.scol')
        ski_file_path = Path(self.rundir, self.galaxy.name + '.ski')

        if star_scol_path.exists() and gas_scol_path.exists() and ski_file_path.exists():
            print("All parameter and data files already in run directory.")
            return
        
        print(">> Creating .scol files...")
        if not (star_scol_path.exists() and gas_scol_path.exists()): 
            star = True
            gas  = True
        elif star_scol_path.exists() and (not gas_scol_path.exists()):
            star = False
            gas  = True
        elif gas_scol_path.exists() and (not star_scol_path.exists()):
            star = True
            gas  = False
        else:
            print("Run directory already contains star.scol and gas.scol. Skipping this step...")
        self._build_scol_files(star=star, gas=gas)

        print(">> Building .ski file...")
        if not ski_file_path.exists():
            self._build_ski_file()
        else:
            print("Run directory already contains ski file. Skipping this step...")
        print("Done.")
        print(f"Wrote initial simulation files for {self.galaxy.name} in {self.rundir}.")       

    def run(self, nthreads=4, nprocesses=2): # number of threads per process
        """Runs simulation"""
        skirt = sm.Skirt()
        ski_path = Path(self.rundir, self.galaxy.name + '.ski')
        print(ski_path)
        print(self.rundir)
        simulation = skirt.execute(
            skiFilePath = str(Path(self.rundir, self.galaxy.name + '.ski')),
            inDirPath = self.rundir,
            outDirPath = self.rundir,
            #skiRelative = True,
            numThreadsPerProcess = nthreads,
            numProcesses = nprocesses)
        simulation = sm.createSimulation(outDirPath=self.rundir, prefix=self.galaxy.name)

def run_simulation(name, rundir, nthreads=4, nprocesses=2): # number of threads per process
    """Run simulation using the already initialized run directory."""
    skirt = sm.Skirt()
    simulation = skirt.execute(
        skiFilePath = str(Path(rundir, name + '.ski')),
        inDirPath = rundir,
        outDirPath = rundir,
        #skiRelative = True,
        numThreadsPerProcess = nthreads,
        numProcesses = nprocesses)
    simulation = sm.createSimulation(outDirPath=rundir, prefix=name)

# Post-processing functions
        
def _get_output_file_path(simulation_name, rundir, instrument_name, file_type, include_dust=True):
    if file_type == "SED":
        file_path = Path(rundir, simulation_name + '_' + instrument_name + '_sed.dat')
    elif file_type == "FITS":
        if include_dust:
            file_path = Path(rundir, simulation_name + '_' + instrument_name + '_total.fits')
        else: 
            file_path = Path(rundir, simulation_name + '_' + instrument_name + '_transparent.fits') 
    return file_path
            
def _get_lum_file_path(simulation_name, rundir):
    return Path(rundir, simulation_name + '_lum_luminosities.dat')
        
def get_band_flux(simulation_name, rundir, instrument_name, bands):
    """
    Args:
        bands : list of str, List of bands to calculate the flux of
            To see available bands, open README.md
        instrument_name : str, Name of the instrument (e.g. "xy", "rand4", <custom instrument name>)
        
    Returns: 
        flux_dict : dict, Dictionary of band fluxes (Jy) in the form 
            {<band1>: {"total": <flux>, "transparent": <flux>}, 
             <band2>: {"total": <flux>, "transparent": <flux>},
             ...
             }

    Note on components:
        "total" sums the contributions of all components of the simulation, taking dust attenuation into account.
        "transparent" removes the effects of dust attenuation, treating the gas as completely transparent.
        "total" corresponds to data collected with dust, and "transparent" corresponds to data collected without dust.
    """
    sed_file = str(_get_output_file_path(simulation_name, rundir, instrument_name, file_type = "SED"))
    print(sed_file)
    wavelengths, total_fluxes, transparent_fluxes = np.loadtxt(
        sed_file, usecols=(0,1,2), unpack=True, encoding="utf-8")
    wavelengths *= u.micron
    total_fluxes *= u.Jy
    transparent_fluxes *= u.Jy

    flux_dict = {}
    if isinstance(bands, list):
        for band in bands:
            builtinBand = pts.band.broadband.BroadBand.builtinBand(nameSegments=band)
            total_flux = builtinBand.convolve(wavelengths,total_fluxes)
            transparent_flux = builtinBand.convolve(wavelengths,transparent_fluxes)
            flux_dict[band] = {}
            flux_dict[band]['total'] = total_flux
            flux_dict[band]['transparent'] = transparent_flux   
    else:
        builtinBand = pts.band.broadband.BroadBand.builtinBand(nameSegments=band)
        total_flux = builtinBand.convolve(wavelengths,total_fluxes)
        transparent_flux = builtinBand.convolve(wavelengths,transparent_fluxes)
        flux_dict[band]['total'] = total_flux
        flux_dict[band]['transparent'] = transparent_flux

    return flux_dict
    
def get_luminosity(simulation_name, rundir, bands=None):
    """
    Args:
        bands : list of str,
            Bands to extract the luminosity of (will return all bands if not specified)
            To see available bands, open README.md.
        
    Returns:
        lum_dict : list of float,
            A dictionary of luminosities corresponding to each band
    """
    lum_file = str(_get_lum_file_path(simulation_name, rundir))
    wavelengths, Lsun = np.loadtxt(lum_file, usecols=(0,2), unpack=True, encoding="utf-8")
    lum_dict = {}
    for i,wavelength in enumerate(wavelengths):
        band = observation_utils.get_predefined_band_name(wavelength)
        if not bands:
            lum_dict[band] = Lsun[i]
        else:
            if band in bands: lum_dict[band] = Lsun[i]
    return lum_dict

    wavelengths,Lsun = np.loadtxt(probe_file, usecols=(0,2), unpack=True, encoding="utf-8")
    lum_dict = {}
    for i,wavelength in enumerate(wavelengths):
        band = get_predefined_band_name(wavelength)
        if band in bands: lum_dict[band] = Lsun[i]
    return lum_dict if len(bands) > 1 else lum_dict[bands[0]]

def get_band_image(
    simulation_name,
    rundir,
    instrument_name, 
    bands=['SDSS_u','SDSS_g','SDSS_r'],
    include_dust=False
    ):
    """
    Args:
        instrument_name : str, Name of the instrument (e.g. "xy", "rand4", <custom instrument name>)
        bands : list of str, List of bands to get the image of
        include_dust : bool,
            True -- Return image with dust attenuation, 
            False -- Return image without dust attenuation
        To see available bands, open README.md.
            
    Returns:
        [0] img_dict : dict, Dictionary of images for each band
            {<band1>: <image>,
             <band2>: <image>,
             ...,
            }
        [1] header : FITS Header object
    """
        
    def _parse_fits(FITS_file):
        with fits.open(FITS_file) as hdul:
            header=hdul[0].header
            imgcube=hdul[0].data
            z_table=hdul[1].data
        wavelengths=np.array(z_table.field(0))
        return header, imgcube, wavelengths

    img_dict={}
    fits_file = _get_output_file_path(
        simulation_name, 
        rundir, 
        instrument_name, 
        file_type="FITS", 
        include_dust=include_dust
    )
    header, imgcube, wavelengths = _parse_fits(fits_file)
    for i,wavelength in enumerate(wavelengths):
        band=observation_utils.get_predefined_band_name(wavelength)
        if band in bands: img_dict[band]=imgcube[i]

    return img_dict, header

def get_mock_image(
    simulation_name,
    rundir,
    instrument_name,
    bands_rgb=['SDSS_r','SDSS_g','SDSS_u'],
    maxden=None, 
    dynrange=100, 
    weights=[1,1,1],
    include_dust=True,
    ):
    """
    Args:
        instrument_name : str, Name of the instrument (e.g. "xy", "rand4", <custom instrument name>)
        bands_rbg : list of str, Bands in order of RGB (default: ['SDSS_u', 'SDSS_g', 'SDSS_r'])
        maxden : float, Maximum density (default: value in the 99th percentile of "r" band image)
        dynrange : int or float, Dynamic range (default: 100)
        weights : list of str or float, Multipliers on the "r-","b-","g-" band images used to adjust color
            (default: [1,1,1])
        include_dust : bool,
            True -- Return image with dust attenuation
            False -- Return image without dust attenuation
        
        To see available bands, open README.md.
        Also, see the FIRE_studio documentation (https://github.com/agurvich/FIRE_studio) for more information.
        
    Returns:
        image : numpy array, Three-band image of galaxy
    """ 
    fits_file = _get_output_file_path(
        simulation_name, 
        rundir, 
        instrument_name, 
        file_type="FITS", 
        include_dust=include_dust
    )    
    img_dict, _ = get_band_image(simulation_name, rundir, instrument_name, bands=bands_rgb)
    if maxden is None: 
        maxden = np.percentile(img_dict[bands_rgb[0]], 99)
    image, massmap = make_threeband_image.make_threeband_image_process_bandmaps(
        img_dict[bands_rgb[0]]*weights[0],
        img_dict[bands_rgb[1]]*weights[1],
        img_dict[bands_rgb[2]]*weights[2],
        maxden=maxden,dynrange=dynrange)

    return image

def save_mock_image(
    image,
    output_dir='',
    savename='mock.png'):
    fig, ax = plt.subplots()
    ax.imshow(image, origin="lower")
    ax.set_xticks([])
    ax.set_yticks([])
    plt.show()
    if savename:
        fig.savefig(
            str(Path(output_dir, savename)),
            dpi=image.shape[0],bbox_inches="tight",pad_inches=0)
    plt.close(fig)
