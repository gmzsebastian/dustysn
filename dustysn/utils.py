from scipy import interpolate
import numpy as np
from astropy import table
from astropy.cosmology import Planck18 as cosmo
from astropy import units as u
import os


# Get directory with reference data
try:
    # If running as a package
    current_file_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    # If running as a script
    current_file_dir = os.getcwd()
data_dir = os.path.join(current_file_dir, 'ref_data')


def calc_distance(redshift):
    """
    Calculate the luminosity distance in cm for a given redshift.

    Parameters
    ----------
    redshift : float
        Redshift of the object

    Returns
    -------
    distance : float
        Luminosity distance in cm
    """
    distance = cosmo.luminosity_distance(z=redshift).to(u.cm)
    return distance


def import_cenwave_table(data_dir=data_dir):
    """
    Get the Astropy table with the central wavelengths of
    the reference filters.

    Central effective wavelengths are in microns and taken from
    the SVO filter profile service.

    Parameters
    ----------
    data_dir : str
        Directory with the cenwave table

    Returns
    -------
    cenwave_table : astropy.table.Table
        Table with Telescope, instrument, filter, and central wavelength columns
    """
    # Read in the cenwave table
    cenwave_table = table.Table.read(os.path.join(data_dir, 'filters', 'cenwaves.txt'), format='ascii')
    return cenwave_table


def import_coefficients(grain_size=0.1, composition='carbon', interp_grain=False,
                        data_dir=data_dir):
    """
    Import the mass absorption coefficients for a given dust composition.
    The data is stored in text files with two columns: wavelength (in microns)
    and mass absorption coefficient (in cm^2/g). Optionally, if grain size is
    not found in the file, it can be interpolated to the desired size.

    Parameters
    ----------
    grain_size : float, default 0.1
        Grain size in microns (either 0.01, 0.1, or 1.0)
    composition : str, default 'carbon'
        Composition of the dust (Either 'carbon' or 'silicate')
    interp_grain : bool, default False
        If True, the grain size will be interpolated to the desired size
        if it is not found in the file.
    data_dir : str
        Reference data directory where the data files are stored

    Returns
    -------
    wave_kappa : array astropy.units.Quantity
        Wavelength array in microns of the dust opacity data
    kappa : array astropy.units.Quantity
        Dust opacity data in cm^2/g
    """

    # Import absorption coefficients file
    if composition == 'carbon':
        file_path = os.path.join(data_dir, 'dust', 'k_Carbon.dat')
    elif composition == 'silicate':
        file_path = os.path.join(data_dir, 'dust', 'k_Silicate.dat')
    else:
        raise ValueError(f"Invalid composition {composition}. Choose 'carbon' or 'silicate'.")
    kappa_data = table.Table.read(file_path, format='ascii')

    # Extract the relevant data
    if not interp_grain:
        try:
            kappa = kappa_data[f'abs_k(a={grain_size}mic)'] * u.cm**2 / u.g
        except KeyError:
            raise KeyError(f"Grain size {grain_size} microns not found in data file.")
        wave_kappa = kappa_data['wavelength(mic)'] * u.micron
    else:
        # Interpolate to the desired grain size
        grain_sizes = [0.01, 0.1, 1.0]
        kappas = []
        for size in grain_sizes:
            kappas.append(kappa_data[f'abs_k(a={size}mic)'] * u.cm**2 / u.g)
        kappas = np.array(kappas)
        wave_kappa = kappa_data['wavelength(mic)'] * u.micron

        # Interpolate to the desired grain size
        interp_func = interpolate.interp1d(grain_sizes, kappas, axis=0, fill_value='extrapolate')
        kappa = interp_func(grain_size)
        kappa = kappa * u.cm**2 / u.g
        wave_kappa = wave_kappa.to(u.micron)

    return wave_kappa, kappa


def interpolate_kappa(wave_kappa, kappa, wave_rest):
    """
    Interpolate the mass absorption coefficients to a new set of
    wavelengths in microns.

    Parameters
    ----------
    wave_kappa : array
        Wavelength array in microns from reference data
    kappa : array
        Dust opacity data in cm^2/g from reference data
    wave_rest : array, astropy.units.Quantity
        Rest wavelength array in microns

    Returns
    -------
    kappa_interp : array
        Interpolated dust opacity data in cm^2/g
    """

    # Create interpolation function
    f_kappa = interpolate.interp1d(wave_kappa.value, kappa.value,
                                   bounds_error=False, fill_value="extrapolate")

    # Interpolate to rest wavelengths
    kappa_interp = f_kappa(wave_rest.value) * u.cm**2 / u.g

    return kappa_interp


def calc_filter_flux(obs_wave, flux_obs, filt_wave, filt_trans):
    """
    Calculate the flux of a model through a filter.

    Parameters
    ----------
    obs_wave : numpy.ndarray
        Wavelength array of the model in microns
    flux_obs : numpy.ndarray
        Flux density of the model in Jansky
    filt_wave : numpy.ndarray
        Filter wavelength array in microns
    filt_trans : numpy.ndarray
        Filter transmission values (0-1)

    Returns
    -------
    float
        Integrated flux in Jansky through the filter
    """

    # Create interpolation function for the model flux
    flux_interp = interpolate.interp1d(obs_wave, flux_obs, bounds_error=False, fill_value=0.0)

    # Keep only filter wavelengths within the model wavelength range
    min_wave = np.min(obs_wave)
    max_wave = np.max(obs_wave)

    mask = (filt_wave >= min_wave) & (filt_wave <= max_wave)

    if not np.any(mask):
        return 0.0  # Filter outside model wavelength range

    filtered_waves = filt_wave[mask]
    filtered_trans = filt_trans[mask]

    # Interpolate model flux at filter wavelengths
    interp_flux = flux_interp(filtered_waves)

    # Replace any NaN values with zeros (from the -inf fill value)
    interp_flux = np.nan_to_num(interp_flux)

    # Calculate filter-weighted flux (flux * transmission)
    weighted_flux = interp_flux * filtered_trans

    # Integrate using trapezoidal rule
    # For photon-counting detector, we need to weight by wavelength
    numerator = np.trapz(weighted_flux * filtered_waves, filtered_waves)
    denominator = np.trapz(filtered_trans * filtered_waves, filtered_waves)

    # Avoid division by zero
    if denominator == 0:
        return 0.0

    # Return the integrated flux in Jansky
    return numerator / denominator


def import_data(filename, data_dir=data_dir):
    """
    Import photometry data from the specified directory
    and format it for use in the MCMC fitting. If the input file
    does not contain the Telescope or Instrument the default
    values are set to 'JWST' and 'MIRI'.

    Parameters
    ----------
    filename : str
        Name of the input data file
    data_dir : str
        Directory where the data file is located

    Returns
    -------
    obs_wave : array
        Wavelengths in microns
    obs_flux : array
        Flux density in Jy
    obs_flux_err : array
        Uncertainty in flux density in Jy
    obs_limits : array
        Boolean array indicating upper limits
    obs_filters : list of str
        List of filter names used in the observations
    obs_wave_filters : list of arrays
        Wavelengths of the filters used in the observations
    obs_trans_filters : list of arrays
        Transmission of the filters used in the observations
    """

    # Read the input data file
    data_file = os.path.join(filename)
    if not os.path.exists(data_file):
        raise FileNotFoundError(f"Data file {data_file} not found.")

    # Read the data using astropy
    input_data = table.Table.read(data_file, format='ascii')

    # Extract relevant columns
    if 'Cenwave' not in input_data.colnames:
        cenwaves = import_cenwave_table(data_dir)
        obs_wave = np.array([cenwaves[cenwaves['Filter'] == f]['Cenwave'][0] for f in input_data['Filter']]) * u.micron
    else:
        obs_wave = input_data['Cenwave'] * u.micron
    obs_flux = input_data['Flux'] * u.Jy
    obs_flux_err = input_data['Flux_err'] * u.Jy
    obs_limits = input_data['UL'] == 'True'
    obs_filters = input_data['Filter']

    if 'Telescope' in input_data.colnames:
        obs_telescope = input_data['Telescope']
    else:
        obs_telescope = ['JWST'] * len(obs_filters)

    if 'Instrument' in input_data.colnames:
        obs_instrument = input_data['Instrument']
    else:
        obs_instrument = ['MIRI'] * len(obs_filters)

    # Get the filter transmissions curves for MIRI
    obs_wave_filters = []
    obs_trans_filters = []
    for i in range(len(obs_filters)):
        # Get filter, telescope, and instrument names
        filter_name = obs_filters[i].strip().upper()
        telescope_name = obs_telescope[i].strip().upper()
        instrument_name = obs_instrument[i].strip().upper()

        # Get filter filenames
        file_name = os.path.join(data_dir, "filters", f"{telescope_name}_{instrument_name}_{filter_name}.txt")
        filter_data = table.Table.read(file_name, format='ascii')

        obs_wave_filters.append(np.array(filter_data['WAVELENGTH']).astype(float))
        obs_trans_filters.append(np.array(filter_data['THROUGHPUT']).astype(float))

    return obs_wave, obs_flux, obs_flux_err, obs_limits, obs_filters, obs_wave_filters, obs_trans_filters


def compute_rhat(chains):
    """
    Calculate Gelman-Rubin R-hat statistic across multiple independent MCMC chains.

    Parameters
    ----------
    chains : ndarray
        Array of shape (n_chains, n_samples), each row is one independent chain.

    Returns
    -------
    float
        Gelman-Rubin R-hat statistic.
    """

    n_chains, n_samples = chains.shape

    if n_chains < 2:
        return np.nan

    # Calculate within-chain variances
    chain_means = np.mean(chains, axis=1)
    chain_vars = np.var(chains, axis=1, ddof=1)

    # Calculate between-chain variance
    grand_mean = np.mean(chain_means)
    B = n_samples * np.sum((chain_means - grand_mean)**2) / (n_chains - 1)

    # Calculate within-chain variance
    W = np.mean(chain_vars)

    # Estimate of marginal posterior variance
    var_hat = ((n_samples - 1) / n_samples) * W + (B / n_samples)

    # Calculate potential scale reduction factor
    if W > 0:
        R_hat = np.sqrt(var_hat / W)
    else:
        R_hat = np.nan

    return R_hat
