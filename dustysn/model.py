from .utils import calc_distance, calc_filter_flux, import_coefficients, interpolate_kappa, import_data
from .plot import plot_corner, plot_trace
import warnings
import scipy.special as sp
import emcee
import matplotlib.pyplot as plt
from astropy.modeling.models import BlackBody
import numpy as np
import astropy.units as u
from astropy import table
import os
from multiprocessing import Pool
plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'font.family': 'serif'})

# Get directory with reference data
try:
    # If running as a package
    current_file_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    # If running locally
    current_file_dir = os.getcwd()
data_dir = os.path.join(current_file_dir, 'ref_data')

# Define priors for the 2-component model
ref_priors = {
    'log_dust_mass_cold': (-6, 1),
    'temp_cold': (20, 2000),
    'log_dust_mass_hot': (-8, 1),
    'temperature_hot': (20, 3000),
    'radius': (1e14, 1e18)
}


def calc_luminosity(rest_wave, kappa_interp, dust_mass, temperature, radius=None,
                    output_units='nu', dust_type='thin'):
    """
    Calculate the intrinsic luminosity from an optically thin dust emission given a
    temperature, dust mass and composition.

    Parameters
    ----------
    rest_wave : array
        Rest-frame wavelength in microns
    kappa_interp : array
        Dust opacity data in cm^2/g interpolated
        to the rest-frame wavelength
    dust_mass : float or Quantity
        Dust mass in solar masses
    temperature : float or Quantity
        Dust temperature in Kelvin
    radius : float or Quantity, optional
        For the optically thick case, the radius of the dust shell in cm.
    output_units : str, default 'nu'
        Output units for the flux ('nu' or 'lambda')
    dust_type : str, default 'thin'
        Type of dust emission ('thin' for optically thin, 'thick' for optically thick).

    Returns
    -------
    luminosity : array
        Luminosity in erg/s/Hz or erg/s/AA
    """

    # Make sure rest_wave and kappa_interp are the same length
    if len(rest_wave) != len(kappa_interp):
        raise ValueError("rest_wave and kappa_interp must have the same length.")

    # Make sure rest_wave is a Quantity in microns
    if not isinstance(rest_wave, u.Quantity):
        rest_wave = rest_wave * u.micron

    # Make sure temperature is a Quantity in Kelvin
    if not isinstance(temperature, u.Quantity):
        temperature = temperature * u.K

    # Making sure dust mass is a Quantity in solar masses
    if not isinstance(dust_mass, u.Quantity):
        dust_mass = dust_mass * u.Msun

    # Make sure radius is a Quantity in cm if provided
    if dust_type == 'thick':
        if radius is None:
            raise ValueError("For optically thick dust, radius must be provided.")
        if not isinstance(radius, u.Quantity):
            radius = radius * u.cm

    # Create blackbody model
    bb = BlackBody(temperature=temperature)

    # Calculate blackbody spectral radiance at observed wavelengths
    # in units of erg/s/cm^2/sr/Hz
    bb_radiance = bb(rest_wave)

    if dust_type == 'thin':
        # Calculate total luminosity in erg/s/Hz
        luminosity_in = dust_mass * kappa_interp * bb_radiance * (4 * np.pi * u.sr)
    elif dust_type == 'thick':
        # Calculate optical depth
        tau = ((3.0/4.0) * kappa_interp * dust_mass / (np.pi * radius**2)).cgs
        # Calculate the escape probability for a spherical shell
        Pesc = (3.0/(4.0*tau))*(1-1/(2.0*tau**2)+(1/tau + 1/(2*tau**2))*np.exp(-2*tau))
        # Calculate luminosity for optically thick dust
        luminosity_in = Pesc * dust_mass * kappa_interp * bb_radiance * (4 * np.pi * u.sr)
    else:
        raise ValueError("dust_type must be 'thin' or 'thick'.")

    # Convert output to desired units
    if output_units == 'nu':
        luminosity = luminosity_in.to(u.erg / u.s / u.Hz)
    elif output_units == 'lambda':
        luminosity = luminosity_in.to(u.erg / u.s / u.AA, equivalencies=u.spectral_density(rest_wave))
    else:
        raise ValueError(f"Invalid output units {output_units}, must be 'nu', or 'lambda'.")

    return luminosity


def calc_flux(rest_wave, luminosity, distance, redshift, output_units='Jy'):
    """
    Calculate the flux density from the luminosity and distance.

    Parameters
    ----------
    rest_wave : array
        Rest-frame wavelength in microns
    luminosity : array
        Luminosity in erg/s/Hz or erg/s/AA
    distance : float or Quantity
        Luminosity distance in cm
    redshift : float
        Redshift of the object
    output_units : str, default 'nu'
        Output units for the flux ('Jy', 'nu', or 'lambda')

    Returns
    -------
    flux : array
        Flux density in Jy, erg/s/cm^2/Hz, or erg/s/cm^2/AA
    """

    # Make sure rest_wave is a Quantity in microns
    if not isinstance(rest_wave, u.Quantity):
        rest_wave = rest_wave * u.micron

    # Make sure distance is a Quantity in cm
    if not isinstance(distance, u.Quantity):
        distance = distance * u.cm

    # Make sure luminosity is a Quantity, return an error if it's not
    if not isinstance(luminosity, u.Quantity):
        raise ValueError("luminosity must be an Astropy Quantity.")

    # Convert luminosity to flux density if luminosity is in erg/s/Hz
    if luminosity.unit.is_equivalent(u.erg / u.s / u.Hz):
        flux_in = luminosity / (4 * np.pi * distance**2) * (1 + redshift)
    elif luminosity.unit.is_equivalent(u.erg / u.s / u.AA):
        flux_in = luminosity / (4 * np.pi * distance**2) / (1 + redshift)
    else:
        raise ValueError("luminosity must be in erg/s/Hz or erg/s/AA.")

    # Calculate observer-frame wavelength
    obs_wave = rest_wave * (1 + redshift)

    # Convert flux to requested units
    if output_units == 'Jy':
        flux = flux_in.to(u.Jy, equivalencies=u.spectral_density(obs_wave))
    elif output_units == 'nu':
        flux = flux_in.to(u.erg / u.s / u.cm / u.cm / u.Hz, equivalencies=u.spectral_density(obs_wave))
    elif output_units == 'lambda':
        flux = flux_in.to(u.erg / u.s / u.cm / u.cm / u.AA, equivalencies=u.spectral_density(obs_wave))
    else:
        raise ValueError(f"Invalid output units {output_units}, must be 'Jy', 'nu', or 'lambda'.")

    return flux


def calc_model_flux(obs_wave, dust_mass, temperature, redshift, radius=None, distance=None, wave_kappa=None,
                    kappa=None, kappa_interp=None, grain_size=0.1, composition='carbon', interp_grain=False,
                    dust_type='thin'):
    """
    Calculate the model flux for a given dust mass and temperature.

    Parameters
    ----------
    obs_wave : array
        Observer-frame wavelength in microns
    dust_mass : float
        Dust mass in solar masses
    temperature : float
        Dust temperature in Kelvin
    redshift : float
        Redshift of the object
    radius : float, optional
        Radius of the dust shell in cm. Only used if dust_type is 'thick'.
    distance : float, optional
        Luminosity distance in cm. If None, will be calculated.
    wave_kappa : array, optional
        Wavelength array in microns from reference data. If None, it will be imported.
    kappa : array, optional
        Dust opacity data in cm^2/g. If None, will be imported based on composition and grain size.
    kappa_interp : array, optional
        Dust opacity interpoalted to the rest-frame wavelength. If None, will be calculated.
    grain_size : float, default 0.1
        Grain size in microns. Only used if kappa is None.
    composition : str, default 'carbon'
        Composition of the dust. Only used if kappa is None.
    interp_grain : bool, default False
        If True, will interpolate the grain size to the observed wavelengths.
    dust_type : str, default 'thin'
        Type of dust emission ('thin' for optically thin, 'thick' for optically thick).

    Returns
    -------
    flux : array
        Model flux in Jy
    """

    # Calculate rest-frame wavelength
    rest_wave = obs_wave / (1 + redshift)

    # If kappa is None, import the dust opacity data
    if kappa_interp is None:
        if (kappa is None) or (wave_kappa is None):
            wave_kappa, kappa = import_coefficients(grain_size=grain_size, composition=composition,
                                                    interp_grain=interp_grain)

        # Interpolate kappa to the rest-frame wavelength
        kappa_interp = interpolate_kappa(wave_kappa, kappa, rest_wave)

    # Calculate distance
    if distance is None:
        distance = calc_distance(redshift)

    # Calculate luminosity
    luminosity = calc_luminosity(rest_wave, kappa_interp, dust_mass, temperature, radius=radius,
                                 output_units='nu', dust_type=dust_type)

    # Calculate flux density
    flux = calc_flux(rest_wave, luminosity, distance, redshift, output_units='Jy')

    return flux


def model_flux(theta, obs_wave, obs_flux, kappa_interp, redshift, distance, radius=None,
               n_components=1, obs_wave_filters=None, obs_trans_filters=None, kappa_interp_hot=None,
               kappa_interp_cold=None, dust_type='thin'):
    """
    Calculate model flux for given parameters with support for one or two dust components.

    Parameters
    ----------
    theta : tuple
        For one component: (log_dust_mass_cold, temp_cold)
        For two components: (log_dust_mass_cold, temp_cold, log_dust_mass_hot, temperature_hot)
    obs_wave : array or Quantity
        Observer-frame wavelength in microns
    obs_flux : array or Quantity
        Observed flux in Jy
    kappa_interp : array or Quantity
        Pre-interpolated dust opacity data in cm^2/g
    redshift : float
        Redshift of the object
    distance : float or Quantity
        Luminosity distance in cm (with units)
    dust_type : str, default 'thin'
        Type of dust emission ('thin' for optically thin, 'thick' for optically thick).
    n_components : int, default 1
        Number of dust components (1 or 2)
    kappa_interp_hot : array or Quantity, optional
        Pre-interpolated dust opacity data for the hot component in cm^2/g.
    kappa_interp_cold : array or Quantity, optional
        Pre-interpolated dust opacity data for the cold component in cm^2/g.
    obs_wave_filters : list of arrays, optional
        Wavelengths of the filters used in the observations

    Returns
    -------
    flux : array
        Model flux in Jansky
    """

    if n_components == 1:
        log_dust_mass_cold, temp_cold = theta
        dust_mass_cold = 10**log_dust_mass_cold

        # Calculate flux for single component
        flux = calc_model_flux(obs_wave, dust_mass_cold, temp_cold, redshift, distance=distance,
                               kappa_interp=kappa_interp, dust_type=dust_type, radius=radius)

    elif n_components == 2:
        log_dust_mass_cold, temp_cold, log_dust_mass_hot, temp_hot = theta
        dust_mass_cold = 10**log_dust_mass_cold
        dust_mass_hot = 10**log_dust_mass_hot

        # If components have different compositions
        if kappa_interp_cold is None:
            kappa_interp_cold = kappa_interp
        if kappa_interp_hot is None:
            kappa_interp_hot = kappa_interp

        # Calculate flux for each component and add them
        flux1 = calc_model_flux(obs_wave, dust_mass_cold, temp_cold, redshift, distance=distance,
                                kappa_interp=kappa_interp_cold, dust_type=dust_type, radius=radius)
        flux2 = calc_model_flux(obs_wave, dust_mass_hot, temp_hot, redshift, distance=distance,
                                kappa_interp=kappa_interp_hot, dust_type=dust_type, radius=radius)
        flux = flux1 + flux2

    else:
        raise ValueError("n_components must be 1 or 2")

    if obs_wave_filters is not None:
        # Apply filter transmission to the model flux
        flux_model = np.zeros(len(obs_flux))
        for i in range(len(obs_wave_filters)):
            # Get the filter transmission for the current filter
            obs_wave_filter = obs_wave_filters[i]
            obs_trans_filter = obs_trans_filters[i]
            output_flux = calc_filter_flux(obs_wave.value, flux.value, obs_wave_filter, obs_trans_filter)
            flux_model[i] = output_flux
        flux_model = flux_model * u.Jy
    else:
        flux_model = flux

    return flux_model


def log_likelihood(theta, obs_wave, obs_flux, obs_flux_err, obs_limits,
                   kappa_interp, redshift, distance, radius=None, n_components=1,
                   obs_wave_filters=None, obs_trans_filters=None,
                   kappa_interp_hot=None, kappa_interp_cold=None,
                   dust_type='thin'):
    """
    Function to calculate the log likelihood of a model, but accounting
    for upper limits in the data.

    Parameters
    ----------
    theta : tuple
        For one component: (log_dust_mass_cold, temp_cold)
        For two components: (log_dust_mass_cold, temp_cold, log_dust_mass_hot, temperature_hot)
    obs_wave : array
        Observed wavelengths in microns
    obs_flux : array
        Observed flux in Jy
    obs_flux_err : array
        Flux uncertainties in Jy
    obs_limits : array
        Boolean array indicating upper limits
    kappa_interp : array
        Pre-interpolated dust opacity data, must
        be the same length as obs_wave
    redshift : float
        Redshift of the object
    distance : float
        Luminosity distance with units of cm
    n_components : int
        Number of dust components (1 or 2)
    obs_wave_filters : list of arrays, optional
        Wavelengths of the filters used in the observations
    obs_trans_filters : list of arrays, optional
        Transmission of the filters used in the observations
    kappa_interp_hot : array, optional
        Pre-interpolated dust opacity data for the hot component in cm^2/g.
    kappa_interp_cold : array, optional
        Pre-interpolated dust opacity data for the cold component in cm^2/g.
    dust_type : str, default 'thin'
        Type of dust emission ('thin' for optically thin, 'thick' for optically thick).

    Returns
    -------
    ln_like : float
        Log likelihood
    """

    # Calculate model flux
    flux_model = model_flux(theta, obs_wave, obs_flux, kappa_interp, redshift, distance, radius,
                            n_components, obs_wave_filters=obs_wave_filters, obs_trans_filters=obs_trans_filters,
                            kappa_interp_hot=kappa_interp_hot, kappa_interp_cold=kappa_interp_cold,
                            dust_type=dust_type)
    ln_like = 0.0

    # Handle detections as usual
    is_detection = ~obs_limits
    if np.any(is_detection):
        det_error = obs_flux.value[is_detection] - flux_model.value[is_detection]
        det_weight = 1.0 / (obs_flux_err.value[is_detection] ** 2)
        ln_like -= 0.5 * np.sum(det_weight * det_error ** 2)

        # Include normalization term for the detections
        ln_like -= 0.5 * np.sum(np.log(2.0 * np.pi * obs_flux_err.value[is_detection] ** 2))

    # Handle upper limits
    if np.any(obs_limits):
        # For each upper limit, calculate the integral term
        for j in np.where(obs_limits)[0]:
            # Calculate how many sigma the model is from the limit
            z = (flux_model.value[j] - obs_flux.value[j]) / obs_flux_err.value[j]

            # Use the complementary error function to calculate the integral term
            # of Equation 8 in https://arxiv.org/pdf/1210.0285
            prob = 0.5 * (1 + sp.erf(z / np.sqrt(2)))
            # Add the log of this probability to the likelihood
            ln_like += np.log(prob) if prob > 0 else -np.inf

    return ln_like


def log_prior(theta, priors, n_components=1, flat_prior=True):
    """
    Calculate the log-prior for the model parameters.

    Parameters
    ----------
    theta : array
        For one component: (log_dust_mass_cold, temp_cold)
        For two components: (log_dust_mass_cold, temp_cold, log_dust_mass_hot, temperature_hot)
    n_components : int
        Number of dust components (1 or 2)
    flat_prior : bool, default True
        If True, use a flat prior on the parameters.
        If False, use a Gaussian prior on the total dust mass.

    Returns
    -------
    log_prior : float
        Log-prior of the model parameters
    """

    # Implement a gaussian prior on the total dust mass
    mu_mass = -1.0  # center at log(Mass/Msun) = -1
    sigma_mass = 1.0  # standard deviation of 1 dex

    if n_components == 1:
        log_dust_mass_cold, temp_cold = theta

        # Calculate total mass (just the cold component in this case)
        log_total_mass = log_dust_mass_cold

        # Check if the parameters are within the prior bounds
        if (priors['log_dust_mass_cold'][0] < log_dust_mass_cold < priors['log_dust_mass_cold'][1] and
                priors['temp_cold'][0] < temp_cold < priors['temp_cold'][1]):
            # Apply Gaussian prior on total dust mass
            mass_prior = -0.5 * ((log_total_mass - mu_mass) / sigma_mass)**2
            mass_prior = -log_total_mass
            if flat_prior:
                return 0.0
            else:
                return mass_prior
        else:
            return -np.inf

    elif n_components == 2:
        log_dust_mass_cold, temp_cold, log_dust_mass_hot, temperature_hot = theta

        # Calculate total mass (sum of cold and hot components in linear space)
        total_mass = 10**log_dust_mass_cold + 10**log_dust_mass_hot
        log_total_mass = np.log10(total_mass)

        # Check if all parameters are within the prior bounds
        if (priors['log_dust_mass_cold'][0] < log_dust_mass_cold < priors['log_dust_mass_cold'][1] and
                priors['temp_cold'][0] < temp_cold < priors['temp_cold'][1] and
                priors['log_dust_mass_hot'][0] < log_dust_mass_hot < priors['log_dust_mass_hot'][1] and
                priors['temperature_hot'][0] < temperature_hot < priors['temperature_hot'][1]):
            # Add a prior to encourage temperature_hot > temp_cold (order)
            if temperature_hot > temp_cold:
                # Apply Gaussian prior on total dust mass
                mass_prior = -0.5 * ((log_total_mass - mu_mass) / sigma_mass)**2
                mass_prior = -log_total_mass
                if flat_prior:
                    return 0.0
                else:
                    return mass_prior
            else:
                return -np.inf
        else:
            return -np.inf

    else:
        raise ValueError("n_components must be 1 or 2")


def log_probability(theta, obs_wave, obs_flux, obs_flux_err, obs_limits,
                    kappa_interp, redshift, distance, radius=None, n_components=1,
                    obs_wave_filters=None, obs_trans_filters=None,
                    kappa_interp_hot=None, kappa_interp_cold=None,
                    dust_type='thin', flat_prior=True, priors=None):
    """
    Calculate the log-probability of the model given the observed data.

    Parameters
    ----------
    theta : array
        For one component: (log_dust_mass_cold, temp_cold)
        For two components: (log_dust_mass_cold, temp_cold, log_dust_mass_hot, temperature_hot)
    obs_wave : array
        Observer-frame wavelength in microns
    obs_flux : array
        Observed flux density in Jy
    obs_flux_err : array
        Uncertainty in the observed flux density in Jy
    obs_limits : array
        Boolean array indicating upper limits
    kappa_interp : array
        Interpolated dust opacity data in cm^2/g
    redshift : float
        Redshift of the object
    distance : float
        Luminosity distance in cm
    radius : float, optional
        Radius of the dust shell in cm. Only used if dust_type is 'thick'.
    n_components : int
        Number of dust components (1 or 2)
    obs_wave_filters : list of arrays, optional
        Wavelengths of the filters used in the observations
    obs_trans_filters : list of arrays, optional
        Transmission of the filters used in the observations
    kappa_interp_hot : array, optional
        Pre-interpolated dust opacity data for the hot component in cm^2/g.
    kappa_interp_cold : array, optional
        Pre-interpolated dust opacity data for the cold component in cm^2/g.
    dust_type : str, default 'thin'
        Type of dust emission ('thin' for optically thin, 'thick' for optically thick).
    flat_prior : bool, default True
        If True, use a flat prior on the parameters.
    priors : dict, optional
        Dictionary containing prior bounds for the parameters.

    Returns
    -------
    log_prob : float
        Log-probability of the model given the observed data
    """

    # First check the prior
    lp = log_prior(theta, priors, n_components, flat_prior)
    if not np.isfinite(lp):
        return -np.inf

    # If prior is finite, add log-likelihood
    return lp + log_likelihood(theta, obs_wave, obs_flux, obs_flux_err, obs_limits,
                               kappa_interp, redshift, distance, radius, n_components,
                               obs_wave_filters, obs_trans_filters, kappa_interp_hot,
                               kappa_interp_cold, dust_type)


def no_warnings(sampler, pos, n_steps, emcee_progress=True):
    """
    Run MCMC sampler while ignoring specific warnings.

    Parameters:
    ----------
    sampler : emcee.EnsembleSampler
        The MCMC sampler to run
    pos : array
        Initial positions for the walkers
    n_steps : int
        Number of steps to take
    emcee_progress : bool, optional
        Whether to show progress bar (default: True)

    Returns:
    -------
    The result of sampler.run_mcmc
    """
    with warnings.catch_warnings():
        # Filter all the specific RuntimeWarnings
        warnings.filterwarnings("ignore", category=RuntimeWarning,
                                message="invalid value encountered in scalar subtract",
                                module="emcee.moves.red_blue")
        warnings.filterwarnings("ignore", category=RuntimeWarning,
                                message="overflow encountered in square")
        warnings.filterwarnings("ignore", category=RuntimeWarning,
                                message="overflow encountered in multiply")
        warnings.filterwarnings("ignore", category=RuntimeWarning,
                                message="overflow encountered in power")
        warnings.filterwarnings("ignore", category=RuntimeWarning,
                                message="divide by zero encountered in log")
        warnings.filterwarnings("ignore", category=RuntimeWarning,
                                message="invalid value encountered in double_scalars")

        # Run the sampler and return its result
        sampler.run_mcmc(pos, n_steps, progress=emcee_progress)

    return sampler


def mcmc_with_sigma_clipping(sampler, pos, n_steps, sigma_clip=2.0, repeats=3, emcee_progress=True):
    """
    Run MCMC with sigma clipping to help convergence.

    After each run, walkers outside the sigma_clip range in parameter space are replaced
    with new walkers drawn from within the sigma_clip range of the current distribution.

    Parameters
    ----------
    sampler : emcee.EnsembleSampler
        MCMC sampler object
    pos : numpy.ndarray
        Initial positions for walkers
    n_steps : int
        Number of steps for each MCMC run
    sigma_clip : float, default=3.0
        Number of sigma to clip walkers at between runs
    repeats : int, default=1
        Number of times to repeat the MCMC process with sigma clipping
    emcee_progress : bool, default=True
        Whether to display progress bar

    Returns
    -------
    sampler : emcee.EnsembleSampler
        Updated sampler with chain from final run
    """
    n_walkers, n_dim = pos.shape

    # First run
    sampler = no_warnings(sampler, pos, n_steps, emcee_progress=emcee_progress)

    # Repeat the process if requested
    for i in range(repeats - 1):
        print(f"Starting MCMC run {i + 2} of {repeats}...")

        # Get the last positions
        last_pos = sampler.chain[:, -1, :]

        # Check for invalid values in the last position
        invalid_mask = np.any(~np.isfinite(last_pos), axis=1)
        if np.any(invalid_mask):
            print(f"Found {np.sum(invalid_mask)} walkers with invalid positions. Replacing them...")
            # Replace invalid walkers with valid ones
            valid_indices = np.where(~invalid_mask)[0]

            # Replace invalid walkers with valid ones plus some noise
            for idx in np.where(invalid_mask)[0]:
                valid_idx = np.random.choice(valid_indices)
                last_pos[idx] = last_pos[valid_idx] + np.random.normal(0, 1e-4, n_dim)

        # Calculate the median and std for each parameter
        medians = np.median(last_pos, axis=0)
        stds = np.std(last_pos, axis=0)

        # Handle cases where std is zero or very small to avoid division by zero
        stds = np.maximum(stds, 1e-10)

        # Identify walkers outside the sigma_clip range
        valid_walkers = np.all(np.abs(last_pos - medians) < sigma_clip * stds, axis=1)
        valid_indices = np.where(valid_walkers)[0]

        # If there are walkers outside the range, replace them
        if len(valid_indices) < n_walkers:
            print(f"Found {n_walkers - len(valid_indices)} walkers outside {sigma_clip}-sigma range.")
            print("Replacing with new walkers drawn from within the clipped distribution...")

            # Create new positions for invalid walkers by sampling from valid ones
            new_pos = np.copy(last_pos)
            invalid_indices = np.where(~valid_walkers)[0]

            for idx in invalid_indices:
                # Sample a valid walker
                valid_idx = np.random.choice(valid_indices)
                # Copy position but add some noise within the acceptable range
                for dim in range(n_dim):
                    new_pos[idx, dim] = last_pos[valid_idx, dim] + \
                                      np.random.normal(0, stds[dim] / sigma_clip)

            # Run the next iteration with the new positions
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=RuntimeWarning,
                                        message="invalid value encountered in scalar subtract",
                                        module="emcee.moves.red_blue")
                sampler = no_warnings(sampler, new_pos, n_steps, emcee_progress=emcee_progress)
        else:
            print("All walkers are within the specified sigma range.")
            # Still run the next iteration with the last positions
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=RuntimeWarning,
                                        message="invalid value encountered in scalar subtract",
                                        module="emcee.moves.red_blue")
                sampler = no_warnings(sampler, last_pos, n_steps, emcee_progress=emcee_progress)

    return sampler


def plot_model(object_name, last_samples, results, obs_wave, obs_flux, obs_flux_err, obs_wave_samples,
               kappa_interp, redshift, distance, wave_kappa_hot, wave_kappa_cold, kappa, composition,
               n_components=1, fig_size=(8, 6), obs_wave_filters=None, obs_trans_filters=None, output_dir='.',
               kappa_interp_hot=None, kappa_interp_cold=None, kappa_hot=None, kappa_cold=None,
               composition_hot=None, composition_cold=None, dust_type='thin', radius=None):
    """
    Plot the MCMC results: corner plot and model fit, supporting 1 or 2 dust components.

    Parameters
    ----------
    object_name : str
        Name of the object being fitted
    last_samples : array
        Last samples from the MCMC chain
    results : dict
        Dictionary with the results of the MCMC fit
    obs_wave : array
        Observer-frame wavelength in microns
    obs_flux : array
        Observed flux density in Jy
    obs_flux_err : array
        Uncertainty in the observed flux density in Jy
    obs_wave_samples : array
        Wavelengths of the filters used in the observations
    obs_limits : array
        Boolean array indicating upper limits
    kappa_interp : array
        Pre-interpolated dust opacity data
    redshift :
        Redshift of the object
    distance : float
        Luminosity distance in cm
    wave_kappa : array
        Wavelengths of the original dust opacity data
    kappa : array
        Dust opacity data in cm^2/g
    n_components : int
        Number of dust components (1 or 2)
    fig_size : tuple
        Figure size
    obs_wave_filters : list of arrays
        Wavelengths of the filters used in the observations
    obs_trans_filters : list of arrays
        Transmission of the filters used in the observations
    output_dir : str
        Directory to save the plots
    kappa_interp_hot : array, optional
        Pre-interpolated dust opacity data for the hot component in cm^2/g.
    kappa_interp_cold : array, optional
        Pre-interpolated dust opacity data for the cold component in cm^2/g.
    kappa_hot : array
        Dust opacity data in cm^2/g for hot component
    kappa_cold : array
        Dust opacity data in cm^2/g for cold component
    composition_hot : str, optional
        Composition of the hot dust component
    composition_cold : str, optional
        Composition of the cold dust component
    dust_type : str, default 'thin'
        Type of dust emission ('thin' for optically thin, 'thick' for optically thick).
    """

    # Plot data and model fit
    print('Plotting Model...')
    plt.figure(figsize=fig_size)

    # Plot observed data points
    plt.errorbar(obs_wave.value, obs_flux.value, yerr=obs_flux_err.value,
                 fmt='o', color='k', label='Observed data', zorder=3, alpha=0.7)

    # Median model at observation points for comparison
    median_params = [results[key][0] for key in results.keys()][:-1]
    median_model = model_flux(median_params, obs_wave_samples, obs_flux, kappa_interp, redshift, distance,
                              radius, n_components, obs_wave_filters=obs_wave_filters,
                              obs_trans_filters=obs_trans_filters, kappa_interp_hot=kappa_interp_hot,
                              kappa_interp_cold=kappa_interp_cold, dust_type=dust_type)

    plt.scatter(obs_wave.value, median_model.value, color='green', marker='s',
                s=80, label='Filter Integrated Model', zorder=2, alpha=0.7)

    # Fine grid for smooth model curve
    wave_dense = np.logspace(np.log10(obs_wave.value.min()*0.5),
                             np.log10(obs_wave.value.max()*1.5), 200) * u.micron

    # If the hot and cold components have different compositions
    if kappa_cold is None:
        kappa_cold = kappa
    if kappa_hot is None:
        kappa_hot = kappa

    if composition_hot is None:
        composition_hot = composition
    if composition_cold is None:
        composition_cold = composition

    # Interpolate kappa to the fine grid's rest wavelengths
    wave_dense_rest = wave_dense / (1 + redshift)
    kappa_dense_hot = interpolate_kappa(wave_kappa_hot, kappa_hot, wave_dense_rest)
    kappa_dense_cold = interpolate_kappa(wave_kappa_cold, kappa_cold, wave_dense_rest)

    # Compute the median model on fine grid
    model_dense = model_flux(median_params, wave_dense, obs_flux, kappa_dense_hot, redshift, distance, radius, n_components,
                             kappa_interp_hot=kappa_dense_hot, kappa_interp_cold=kappa_dense_cold, dust_type=dust_type)
    total_M = results['total_dust_mass'][0]

    if n_components == 1:
        extra_label = ("\n" + r"$T = "
                       f"{results['temp_cold'][0]:.0f}"
                       r"^{+" + f"{results['temp_cold'][1]:.0f}" + r"}"
                       r"_{-" + f"{results['temp_cold'][2]:.0f}" + r"} \, {\rm K}$")
    else:
        extra_label = ''

    if total_M < 1e-2:
        total_label = (
            r"Best fit = $"
            f"{total_M:.1e}"
            r"^{+" + f"{results['total_dust_mass'][1]:.1e}" + r"}"
            r"_{-" + f"{results['total_dust_mass'][2]:.1e}" + r"} \, M_\odot$" +
            extra_label
        )
    else:
        total_label = (
            r"Best fit = $"
            f"{total_M:.2f}"
            r"^{+" + f"{results['total_dust_mass'][1]:.2f}" + r"}"
            r"_{-" + f"{results['total_dust_mass'][2]:.2f}" + r"} \, M_\odot$" +
            extra_label
        )
    plt.plot(wave_dense.value, model_dense.value, 'g', linestyle='-', alpha=0.9, linewidth=2, zorder=1,
             label=total_label)

    # Plot individual components for two-component model
    if n_components == 2:
        # Component 1
        comp1_params = (median_params[0], median_params[1])
        comp1_model = model_flux(comp1_params, wave_dense, obs_flux, kappa_dense_cold, redshift, distance,
                                 radius, n_components=1, dust_type=dust_type)

        # Extract values
        log_M_cold, log_M_cold_upper, log_M_cold_lower = results['log_dust_mass_cold']
        T_cold, T_cold_upper, T_cold_lower = results['temp_cold']
        log_M_hot, log_M_hot_upper, log_M_hot_lower = results['log_dust_mass_hot']
        T_hot, T_hot_upper, T_hot_lower = results['temperature_hot']

        # Convert log mass to actual mass in scientific notation
        M_cold = 10**log_M_cold
        M_cold_upper = M_cold * (10**log_M_cold_upper - 1)  # Convert log error to linear
        M_cold_lower = M_cold * (1 - 10**(-log_M_cold_lower))

        M_hot = 10**log_M_hot
        M_hot_upper = M_hot * (10**log_M_hot_upper - 1)
        M_hot_lower = M_hot * (1 - 10**(-log_M_hot_lower))

        if radius is None:
            radius_label = ''
        else:
            radius_label = f'\nR = {radius:.1e} cm'

        # Create formatted strings with appropriate precision
        if M_cold < 1e-2:
            cold_label = (
                f"Cold {dust_type} {composition_cold}\n"
                r"$M = "
                f"{M_cold:.1e}"
                r"^{+" + f"{M_cold_upper:.1e}" + r"}"
                r"_{-" + f"{M_cold_lower:.1e}" + r"} \, M_\odot$"
                "\n"
                r"$T = "
                f"{T_cold:.0f}"
                r"^{+" + f"{T_cold_upper:.0f}" + r"}"
                r"_{-" + f"{T_cold_lower:.0f}" + r"} \, {\rm K}$" +
                radius_label
            )
        else:
            cold_label = (
                f"Cold {dust_type} {composition_cold}\n"
                r"$M = "
                f"{M_cold:.2f}"
                r"^{+" + f"{M_cold_upper:.2f}" + r"}"
                r"_{-" + f"{M_cold_lower:.2f}" + r"} \, M_\odot$"
                "\n"
                r"$T = "
                f"{T_cold:.0f}"
                r"^{+" + f"{T_cold_upper:.0f}" + r"}"
                r"_{-" + f"{T_cold_lower:.0f}" + r"} \, {\rm K}$" +
                radius_label
            )
        if M_hot < 1e-2:
            hot_label = (
                f"Hot {dust_type} {composition_hot}\n"
                r"$M = "
                f"{M_hot:.1e}"
                r"^{+" + f"{M_hot_upper:.1e}" + r"}"
                r"_{-" + f"{M_hot_lower:.1e}" + r"} \, M_\odot$"
                "\n"
                r"$T = "
                f"{T_hot:.0f}"
                r"^{+" + f"{T_hot_upper:.0f}" + r"}"
                r"_{-" + f"{T_hot_lower:.0f}" + r"} \, {\rm K}$" +
                radius_label
            )
        else:
            hot_label = (
                f"Hot {dust_type} {composition_hot}\n"
                r"$M = "
                f"{M_hot:.2f}"
                r"^{+" + f"{M_hot_upper:.2f}" + r"}"
                r"_{-" + f"{M_hot_lower:.2f}" + r"} \, M_\odot$"
                "\n"
                r"$T = "
                f"{T_hot:.0f}"
                r"^{+" + f"{T_hot_upper:.0f}" + r"}"
                r"_{-" + f"{T_hot_lower:.0f}" + r"} \, {\rm K}$" +
                radius_label
            )
        plt.plot(wave_dense.value, comp1_model.value, color='b', linestyle='--', alpha=0.7, linewidth=1.5,
                 label=cold_label)
        # Component 2
        comp2_params = (median_params[2], median_params[3])
        comp2_model = model_flux(comp2_params, wave_dense, obs_flux, kappa_dense_hot, redshift, distance,
                                 radius, n_components=1, dust_type=dust_type)
        plt.plot(wave_dense.value, comp2_model.value, color='r', linestyle='--', alpha=0.7, linewidth=1.5,
                 label=hot_label)

    # Visualize uncertainty
    # Calculate model for a range of samples to create confidence interval
    model_values = np.zeros((len(last_samples), len(wave_dense)))

    for i, idx in enumerate(last_samples):
        sample_model = model_flux(idx, wave_dense, obs_flux, kappa_dense_hot, redshift, distance, radius, n_components,
                                  kappa_interp_hot=kappa_dense_hot, kappa_interp_cold=kappa_dense_cold, dust_type=dust_type)
        model_values[i] = sample_model.value
        plt.plot(wave_dense.value, sample_model.value, linewidth=0.1, color='green', alpha=0.1)

    # Calculate percentiles for lower and upper bounds
    lower_bound, upper_bound = np.percentile(model_values, [15.87, 84.13], axis=0)

    # Determine the optimal y-axis limits
    ymin = 10 ** (np.log10(np.nanmin(obs_flux.value)) - 0.2)
    ymax = 10 ** (np.log10(np.nanmax(obs_flux.value)) + 0.2)

    # Determine the optimal x-axis limits
    xmin = np.min([5, np.min(obs_wave.value) * 0.9])
    xmax = np.max([30, np.max(obs_wave.value) * 1.1])

    # Plot confidence region
    plt.fill_between(wave_dense.value, lower_bound, upper_bound, linewidth=0,
                     color='green', alpha=0.2)

    # Set labels and scales
    plt.xlabel(r'Observed Wavelength ($\mu$m)')
    plt.ylabel('Flux Density (Jy)')
    plt.yscale('log')
    plt.legend(loc='lower right', fontsize=10)
    plt.ylim(ymin, ymax)
    plt.xlim(xmin, xmax)
    plt.tight_layout()
    output_filename = os.path.join(output_dir, f"{object_name}_{n_components}_model_fit.pdf")
    plt.savefig(output_filename, bbox_inches='tight')
    plt.clf()
    plt.close('all')


def fit_dust_model(obs_wave, obs_flux, obs_flux_err, obs_limits, redshift, object_name,
                   composition='carbon', grain_size=0.1, n_components=1, n_walkers=32,
                   n_steps=1000, burn_in=0.75, n_cores=1, sigma_clip=2, repeats=3, emcee_progress=True,
                   obs_wave_filters=None, obs_trans_filters=None, n_filter_samples=1000, plot=False, plots_model=False,
                   plots_corner=False, plots_trace=False, output_dir='.', initial_pos=None, composition_hot=None,
                   composition_cold=None, dust_type='thin', radius=None, flat_prior=True, priors=None):
    """
    Run MCMC sampler without multiprocessing.

    Parameters
    ----------
    obs_wave : array
        Observer-frame wavelength in microns
    obs_flux : array
        Observed flux density in Jy
    obs_flux_err : array
        Uncertainty in the observed flux density in Jy
    obs_limits : array
        Boolean array indicating upper limits
    redshift : float
        Redshift of the object
    composition : str, default 'carbon'
        Composition of the dust (e.g., 'carbon', 'silicate')
    grain_size : float, default 0.1
        Grain size in microns
    n_components : int
        Number of dust components to model (1 or 2)
    n_walkers : int
        Number of walkers
    n_steps : int
        Number of steps for each walker
    burn_in : float
        Fraction of steps to discard as burn-in
    n_cores : int
        Number of threads
    sigma_clip : float
        Number of sigma to clip walkers at between runs
    repeats : int
        Number of times to repeat the MCMC process with sigma clipping
    emcee_progress : bool
        Whether to display progress bar
    obs_wave_filters : list of arrays
        Wavelengths of the filters used in the observations
    obs_trans_filters : list of arrays
        Transmission of the filters used in the observations
    n_filter_samples : int
        Number of samples to draw for the filter wavelengths
    plot : bool, default False
        Plot the model, corner plot, and traces
    plots_model : bool
        Plot model
    plots_corner : bool
        Plot corner plot
    plots_trace : bool
        Plot traces for parmeters
    output_dir : str
        Directory to save the plots and parameter files
    initial_pos : dict
        Initial positions for the walkers. If None, priors will be used to generate them.
    composition_hot : str, optional
        Composition of the hot dust component (if n_components=2)
    composition_cold : str, optional
        Composition of the cold dust component (if n_components=2)
    dust_type : str, default 'thin'
        Type of dust emission ('thin' for optically thin, 'thick' for optically thick).
    radius : float, optional
        Radius of the dust grains in microns (if applicable).
    flat_prior : bool, default True
        If True, use a flat prior on the parameters.
    priors : dict, optional
        Dictionary with prior bounds for the parameters. If None, default priors will be used.

    Returns
    -------
    sampler : emcee.EnsembleSampler
        The MCMC sampler object
    samples : array
        MCMC samples
    """
    # Number of parameters
    np.random.seed(42)

    # Calculate the distance in cm
    distance = calc_distance(redshift)

    # If composition_hot or composition_cold are not provided, use the same as composition
    if composition_cold is None:
        composition_cold = composition
    if composition_hot is None:
        composition_hot = composition

    # Import dust data
    wave_kappa, kappa = import_coefficients(grain_size=grain_size, composition=composition,
                                            interp_grain=False)

    # For hot component
    wave_kappa_hot, kappa_hot = import_coefficients(grain_size=grain_size, composition=composition_hot,
                                                    interp_grain=False)

    # For cold component
    wave_kappa_cold, kappa_cold = import_coefficients(grain_size=grain_size, composition=composition_cold,
                                                      interp_grain=False)

    if obs_wave_filters is not None:
        # Pre-calculate interpolated opacities for the filter wavelengths
        mins, maxs = np.array([(np.min(i), np.max(i)) for i in obs_wave_filters]).T
        min_wave = np.min(mins)
        max_wave = np.max(maxs)
        sampled_waves = np.linspace(min_wave, max_wave, n_filter_samples)
        rest_wave = sampled_waves / (1 + redshift) * u.micron
        kappa_interp = interpolate_kappa(wave_kappa, kappa, rest_wave)
        kappa_interp_hot = interpolate_kappa(wave_kappa_hot, kappa_hot, rest_wave)
        kappa_interp_cold = interpolate_kappa(wave_kappa_cold, kappa_cold, rest_wave)
        obs_wave_samples = sampled_waves * u.micron
    else:
        # Pre-calculate interpolated opacities for the reference wavelengths
        rest_wave = obs_wave / (1 + redshift)
        kappa_interp = interpolate_kappa(wave_kappa, kappa, rest_wave)
        kappa_interp_hot = interpolate_kappa(wave_kappa_hot, kappa_hot, rest_wave)
        kappa_interp_cold = interpolate_kappa(wave_kappa_cold, kappa_cold, rest_wave)
        obs_wave_samples = obs_wave

    # Set up initial positions for walkers
    def create_prior(n_walkers, initial_pos):
        if n_components == 1:
            log_dust_mass_cold = np.random.uniform(initial_pos['log_dust_mass_cold'][0],
                                                   initial_pos['log_dust_mass_cold'][1], n_walkers)
            temp_cold = np.random.uniform(initial_pos['temp_cold'][0],
                                          initial_pos['temp_cold'][1], n_walkers)
            pos = np.array([log_dust_mass_cold, temp_cold]).T
            return pos
        elif n_components == 2:
            log_dust_mass_cold = np.random.uniform(initial_pos['log_dust_mass_cold'][0],
                                                   initial_pos['log_dust_mass_cold'][1], n_walkers)
            temp_cold = np.random.uniform(initial_pos['temp_cold'][0],
                                          initial_pos['temp_cold'][1], n_walkers)
            log_dust_mass_hot = np.random.uniform(initial_pos['log_dust_mass_hot'][0],
                                                  initial_pos['log_dust_mass_hot'][1], n_walkers)
            temperature_hot = np.random.uniform(initial_pos['temperature_hot'][0],
                                                initial_pos['temperature_hot'][1], n_walkers)
            pos = np.array([log_dust_mass_cold, temp_cold,
                           log_dust_mass_hot, temperature_hot]).T
            return pos
        else:
            raise ValueError("n_components must be 1 or 2")

    if priors is None:
        priors = ref_priors
    if initial_pos is None:
        initial_pos = priors

    # Create initial positions
    pos_in = create_prior(n_walkers, initial_pos)
    pos_out = pos_in[0:1]
    while len(pos_out) < n_walkers:
        pos = pos_in[[np.isfinite(log_prior(i, priors, n_components=n_components)) for i in pos_in]]
        pos_out = np.append(pos_out, pos, axis=0)

    # Crop to correct length
    if len(pos_out) != n_walkers:
        pos = pos_out[1:n_walkers+1]
    else:
        pos = pos_out

    # Set up arguments for the log-probability function
    args = (obs_wave_samples, obs_flux, obs_flux_err, obs_limits, kappa_interp, redshift, distance, radius,
            n_components, obs_wave_filters, obs_trans_filters, kappa_interp_hot, kappa_interp_cold, dust_type,
            flat_prior, priors)

    # Run MCMC without multiprocessing
    n_dim = pos.shape[1]
    # Create a multiprocessing pool if n_cores > 1
    if n_cores > 1:
        with Pool(processes=n_cores) as pool:
            sampler = emcee.EnsembleSampler(n_walkers, n_dim, log_probability, args=args, pool=pool)

            # Show progress bar during sampling
            print("Running MCMC with parallel processing using", n_cores, "cores...")
            sampler = mcmc_with_sigma_clipping(sampler, pos, n_steps,
                                               sigma_clip=sigma_clip,
                                               repeats=repeats,
                                               emcee_progress=emcee_progress)
    else:
        # Use the existing non-parallel approach if n_cores=1
        sampler = emcee.EnsembleSampler(n_walkers, n_dim, log_probability, args=args)

        # Show progress bar during sampling
        print("Running MCMC without parallel processing...")
        sampler = mcmc_with_sigma_clipping(sampler, pos, n_steps,
                                           sigma_clip=sigma_clip,
                                           repeats=repeats,
                                           emcee_progress=emcee_progress)

    # Only consider the last bit of the chain for parameter estimation
    if repeats > 1:
        samples_crop = sampler.chain[:, -n_steps:].reshape((-1, n_dim))
    else:
        samples_crop = sampler.chain[:, -int(n_steps*(1-burn_in)):, :].reshape((-1, n_dim))
    last_samples = sampler.chain[:, -1, :]

    # Obtain the parametrs of the best fit
    if n_components == 1:
        log_dust_mass_cold_mcmc, temp_cold_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                                                      zip(*np.percentile(samples_crop, [15.87, 50, 84.13], axis=0)))
        # Calculate total dust mass
        total_dust_mass_samples = 10 ** samples_crop[:, 0]
        total_dust_mass_percentiles = np.percentile(total_dust_mass_samples, [15.87, 50, 84.13], axis=0)
        total_dust_mass_mcmc = (total_dust_mass_percentiles[1],
                                total_dust_mass_percentiles[2] - total_dust_mass_percentiles[1],
                                total_dust_mass_percentiles[1] - total_dust_mass_percentiles[0])
        results = {
            'log_dust_mass_cold': log_dust_mass_cold_mcmc,
            'temp_cold': temp_cold_mcmc,
            'total_dust_mass': total_dust_mass_mcmc
        }
    elif n_components == 2:
        log_dust_mass_cold_mcmc, temp_cold_mcmc, log_dust_mass_hot_mcmc, temperature_hot_mcmc = \
            map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                zip(*np.percentile(samples_crop, [15.87, 50, 84.13], axis=0)))
        # Calculate total dust mass
        total_dust_mass_samples = 10 ** samples_crop[:, 0] + 10 ** samples_crop[:, 2]
        total_dust_mass_percentiles = np.percentile(total_dust_mass_samples, [15.87, 50, 84.13], axis=0)
        total_dust_mass_mcmc = (total_dust_mass_percentiles[1],
                                total_dust_mass_percentiles[2] - total_dust_mass_percentiles[1],
                                total_dust_mass_percentiles[1] - total_dust_mass_percentiles[0])
        results = {
            'log_dust_mass_cold': log_dust_mass_cold_mcmc,
            'temp_cold': temp_cold_mcmc,
            'log_dust_mass_hot': log_dust_mass_hot_mcmc,
            'temperature_hot': temperature_hot_mcmc,
            'total_dust_mass': total_dust_mass_mcmc
        }

    # Make sure the directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Master plot command
    if plot:
        plots_corner = True
        plots_model = True
        plots_trace = True

    # Create requested plots
    if plots_model:
        plot_model(object_name, last_samples, results, obs_wave, obs_flux, obs_flux_err, obs_wave_samples,
                   kappa_interp, redshift, distance, wave_kappa_hot, wave_kappa_cold, kappa, composition,
                   n_components=n_components, obs_wave_filters=obs_wave_filters, obs_trans_filters=obs_trans_filters,
                   output_dir=output_dir, kappa_interp_hot=kappa_interp_hot, kappa_interp_cold=kappa_interp_cold,
                   kappa_hot=kappa_hot, kappa_cold=kappa_cold, composition_cold=composition_cold,
                   composition_hot=composition_hot, dust_type=dust_type, radius=radius)

    if plots_corner:
        plot_corner(object_name, samples_crop, results, n_components, output_dir=output_dir)

    if plots_trace:
        # Plot trace for each parameter
        params = list(results.keys())[:-1]
        for i, param in enumerate(params):
            is_log = False
            plot_trace(sampler.chain[:, :, i],
                       results[param],
                       results[param],
                       priors[param][0],
                       priors[param][1],
                       f"{param}",
                       f"{param}",
                       is_log,
                       n_steps,
                       burn_in,
                       repeats,
                       object_name,
                       n_components,
                       output_dir=output_dir)

    # Create lists for the table columns
    parameters = []
    median_values = []
    upper_values = []
    lower_values = []

    # Extract data from the dictionary
    for parameter, values in results.items():
        parameters.append(parameter)
        median_values.append(values[0])
        upper_values.append(values[1])
        lower_values.append(values[2])

    # Create Astropy Table
    output_table = table.Table([parameters, median_values, upper_values, lower_values],
                               names=('parameter', 'median', 'upper', 'lower'))
    output_filename = os.path.join(output_dir, f'parameters_{object_name}_{n_components}.txt')
    output_table.write(output_filename, format='ascii', overwrite=True)

    return results, last_samples


def compare_models(obs_wave, obs_flux, obs_flux_err, obs_limits, redshift, object_name,
                   results_1, results_2, composition1='carbon', composition2='carbon', grain_size=0.1,
                   obs_wave_filters=None, obs_trans_filters=None, fig_size=(8, 6),
                   output_dir='.', plot_comparison=True, composition_hot=None, composition_cold=None,
                   dust_type='thin', radius=None):
    """
    Compare 1-component vs 2-component dust models.

    Parameters
    ----------
    obs_wave : array
        Observer-frame wavelength in microns
    obs_flux : array
        Observed flux density in Jy
    obs_flux_err : array
        Uncertainty in the observed flux density in Jy
    obs_limits : array
        Boolean array indicating upper limits
    redshift : float
        Redshift of the object
    object_name : str
        Name of the object being fitted
    results_1 : dict
        Dictionary with the results of the 1-component MCMC fit
    results_2 : dict
        Dictionary with the results of the 2-component MCMC fit
    composition1 : str, default 'carbon'
        Composition of the dust for the 1-component model
    composition2 : str, default 'carbon'
        Composition of the dust for the 2-component model
    grain_size : float, default 0.1
        Grain size in microns
    obs_wave_filters : list of arrays
        Wavelengths of the filters used in the observations
    obs_trans_filters : list of arrays
        Transmission of the filters used in the observations
    fig_size : tuple
        Figure size
    output_dir : str
        Directory to save the plots
    plot_comparison : bool
        Whether to plot and save the comparison results
    composition_hot : str, optional
        Composition of the hot dust component (if n_components=2)
    composition_cold : str, optional
        Composition of the cold dust component (if n_components=2)
    dust_type : str, default 'thin'
        Type of dust emission ('thin' for optically thin, 'thick' for optically thick).

    Returns
    -------
    results_dict : dict
        Dictionary containing results for both models and comparison metrics
    """

    # Calculate the distance in cm
    distance = calc_distance(redshift)

    # Calculate likelihood for the first model
    wave_kappa_1, kappa_1 = import_coefficients(grain_size=grain_size, composition=composition1,
                                                interp_grain=False)
    obs_wave_rest_1 = obs_wave / (1 + redshift)
    kappa_interp_1 = interpolate_kappa(wave_kappa_1, kappa_1, obs_wave_rest_1)

    # Get the best parameters for the first model
    best_params_1 = [results_1[key][0] for key in results_1.keys()][:-1]
    log_like_1 = log_likelihood(best_params_1, obs_wave, obs_flux, obs_flux_err, obs_limits,
                                kappa_interp_1, redshift, distance, radius, n_components=1,
                                obs_wave_filters=obs_wave_filters, obs_trans_filters=obs_trans_filters,
                                kappa_interp_hot=None, kappa_interp_cold=None, dust_type=dust_type)

    # Calculate likelihood for the second model
    if composition_cold is None:
        composition_cold = composition2
    if composition_hot is None:
        composition_hot = composition2

    # Interpolate cold component
    wave_kappa_2_cold, kappa_2_cold = import_coefficients(grain_size=grain_size, composition=composition_cold,
                                                          interp_grain=False)
    obs_wave_rest_2_cold = obs_wave / (1 + redshift)
    kappa_interp_cold = interpolate_kappa(wave_kappa_2_cold, kappa_2_cold, obs_wave_rest_2_cold)

    # Interpolate hot component
    wave_kappa_2_hot, kappa_2_hot = import_coefficients(grain_size=grain_size, composition=composition_hot,
                                                        interp_grain=False)
    obs_wave_rest_2_hot = obs_wave / (1 + redshift)
    kappa_interp_hot = interpolate_kappa(wave_kappa_2_hot, kappa_2_hot, obs_wave_rest_2_hot)

    # Get the best parameters for the second model
    best_params_2 = [results_2[key][0] for key in results_2.keys()][:-1]
    log_like_2 = log_likelihood(best_params_2, obs_wave, obs_flux, obs_flux_err, obs_limits,
                                kappa_interp_hot, redshift, distance, radius, n_components=2,
                                obs_wave_filters=obs_wave_filters, obs_trans_filters=obs_trans_filters,
                                kappa_interp_hot=kappa_interp_hot, kappa_interp_cold=kappa_interp_cold,
                                dust_type=dust_type)

    # Calculate AIC and BIC
    n_data = len(obs_flux)
    n_params_1 = len(best_params_1)  # one component with log_mass and temperature
    n_params_2 = len(best_params_2)  # two components with log_mass and temperature

    aic_1comp = 2 * n_params_1 - 2 * log_like_1
    aic_2comp = 2 * n_params_2 - 2 * log_like_2
    delta_aic = aic_1comp - aic_2comp

    bic_1comp = np.log(n_data) * n_params_1 - 2 * log_like_1
    bic_2comp = np.log(n_data) * n_params_2 - 2 * log_like_2
    delta_bic = bic_1comp - bic_2comp

    # Print comparison
    print("\nModel Comparison:")
    print(f"1-component model: Log-likelihood = {log_like_1:.2f}, AIC = {aic_1comp:.2f}, BIC = {bic_1comp:.2f}")
    print(f"2-component model: Log-likelihood = {log_like_2:.2f}, AIC = {aic_2comp:.2f}, BIC = {bic_2comp:.2f}")
    print(f"Delta AIC (1comp - 2comp) = {delta_aic:.2f} (positive favors 2-component model)")
    print(f"Delta BIC (1comp - 2comp) = {delta_bic:.2f} (positive favors 2-component model)")

    if delta_aic > 0 and delta_bic > 0:
        print("Both AIC and BIC favor the 2-component model.")
    elif delta_aic > 0 and delta_bic <= 0:
        print("AIC favors the 2-component model, but BIC (which penalizes complexity more) favors the 1-component model.")
    elif delta_aic <= 0 and delta_bic <= 0:
        print("Both AIC and BIC favor the 1-component model.")

    # Create a figure that shows both models
    if plot_comparison:
        plt.figure(figsize=fig_size)

        # Plot observed data
        plt.errorbar(obs_wave.value, obs_flux.value, yerr=obs_flux_err.value,
                     fmt='o', color='black', label='Observed data', zorder=3)

        # Fine grid for smooth model curves
        wave_dense = np.logspace(np.log10(obs_wave.value.min()*0.5),
                                 np.log10(obs_wave.value.max()*1.5), 200) * u.micron

        # Interpolate kappa to the fine grid's rest wavelengths
        wave_dense_rest = wave_dense / (1 + redshift)
        kappa_dense_1 = interpolate_kappa(wave_kappa_1, kappa_1, wave_dense_rest)
        kappa_dense_hot = interpolate_kappa(wave_kappa_2_cold, kappa_2_cold, wave_dense_rest)
        kappa_dense_cold = interpolate_kappa(wave_kappa_2_hot, kappa_2_hot, wave_dense_rest)

        # Plot 1-component model
        model_1comp = model_flux(best_params_1, wave_dense, obs_flux, kappa_dense_1,
                                 redshift, distance, radius, n_components=1, dust_type=dust_type)
        plt.plot(wave_dense.value, model_1comp.value, 'b-', linewidth=2,
                 label='1-component', alpha=0.7)

        # Plot 2-component model
        model_2comp = model_flux(best_params_2, wave_dense, obs_flux, kappa_dense_hot,
                                 redshift, distance, radius, n_components=2,
                                 kappa_interp_hot=kappa_dense_hot, kappa_interp_cold=kappa_dense_cold,
                                 dust_type=dust_type)
        plt.plot(wave_dense.value, model_2comp.value, 'r-', linewidth=2,
                 label='2-component', alpha=0.7)

        # Plot individual components of the 2-component model
        # Cold component
        comp1_params = (best_params_2[0], best_params_2[1])
        comp1_model = model_flux(comp1_params, wave_dense, obs_flux, kappa_dense_cold,
                                 redshift, distance, radius, n_components=1, dust_type=dust_type)
        plt.plot(wave_dense.value, comp1_model.value, 'r--', alpha=0.5, linewidth=1.5,
                 label='Cold')

        # Hot component
        comp2_params = (best_params_2[2], best_params_2[3])
        comp2_model = model_flux(comp2_params, wave_dense, obs_flux, kappa_dense_hot,
                                 redshift, distance, radius, n_components=1, dust_type=dust_type)

        # Determine the optimal y-axis limits
        ymin = 10 ** (np.log10(np.nanmin(obs_flux.value)) - 0.2)
        ymax = 10 ** (np.log10(np.nanmax(obs_flux.value)) + 0.2)

        # Determine the optimal x-axis limits
        xmin = np.min([5, np.min(obs_wave.value) * 0.9])
        xmax = np.max([30, np.max(obs_wave.value) * 1.1])

        plt.plot(wave_dense.value, comp2_model.value, 'r:', alpha=0.5, linewidth=1.5,
                 label='Hot')
        plt.title(rf'AIC = {delta_aic:.2f} $-$ BIC = {delta_bic:.2f}')
        plt.xlabel(r'Observed Wavelength ($\mu$m)')
        plt.ylabel('Flux Density (Jy)')
        plt.yscale('log')
        plt.legend(loc='lower right')
        plt.ylim(ymin, ymax)
        plt.xlim(xmin, xmax)

        plt.tight_layout()
        output_filename = os.path.join(output_dir, f"comparison_{object_name}.pdf")
        plt.savefig(output_filename, bbox_inches='tight')
        plt.clf()
        plt.close('all')

    # Return combined results
    return {
        '1comp': results_1,
        '2comp': results_2,
        'comparison': {
            'log_like_1': log_like_1,
            'log_like_2': log_like_2,
            'aic_1comp': aic_1comp,
            'aic_2comp': aic_2comp,
            'delta_aic': delta_aic,
            'bic_1comp': bic_1comp,
            'bic_2comp': bic_2comp,
            'delta_bic': delta_bic
        }
    }


def full_model(filename, object_name, redshift, n_steps, n_walkers, composition='carbon',
               grain_size=0.1, n_cores=1, sigma_clip=2, repeats=3, emcee_progress=True,
               plot=False, plots_model=False, plots_corner=False, plots_trace=False, output_dir='.',
               burn_in=0.75, n_filter_samples=1000, initial_pos=None, priors=None, composition_hot=None,
               composition_cold=None, dust_type='thin', radius=None, flat_prior=True):
    """
    Run the full model fitting process.

    Parameters
    ----------
    filename : str
        Path to the input data file
    object_name : str
        Name of the object being fitted
    redshift : float
        Redshift of the object
    n_steps : int
        Number of steps for each MCMC run
    n_walkers : int
        Number of walkers
    composition : str, default 'carbon'
        Composition of the dust (e.g., 'carbon', 'silicate')
    grain_size : float, default 0.1
        Grain size in microns
    n_cores : int, default 1
        Number of threads for multiprocessing
    sigma_clip : float, default 2.0
        Number of sigma to clip walkers at between runs
    repeats : int, default 3
        Number of times to repeat the MCMC process with sigma clipping
    emcee_progress : bool, default True
        Whether to display progress bar during MCMC sampling
    plot : bool, default True
        Whether to plot and save all the diagnosis plots
    plots_model : bool
        Plot model
    plots_corner : bool
        Plot corner plot
    plots_trace : bool
        Plot traces for parmeters
    output_dir : str, default '.'
        Directory to save the plots and parameter files
    burn_in : float, default 0.75
        Fraction of steps to discard as burn-in
    n_filter_samples : int, default 1000
        Number of samples to draw for the filter wavelengths
    initial_pos : dict, default None
        Initial positions for the walkers. If None, priors will be used to generate them.
    priors : dict, default None
        Dictionary containing prior ranges for the parameters.
    composition_hot : str, optional
        Composition of the hot dust component (if n_components=2)
    composition_cold : str, optional
        Composition of the cold dust component (if n_components=2)
    dust_type : str, default 'thin'
        Type of dust emission ('thin' for optically thin, 'thick' for optically thick).
    radius : float, optional
        Radius of the dust grains in microns (if applicable).
    flat_prior : bool, default True
        If True, use a flat prior on the parameters.

    Returns
    -------
    results : dict
        Dictionary containing results for both models and comparison metrics if applicable.
    """

    # Import data from file
    obs_wave, obs_flux, obs_flux_err, obs_limits, obs_filters, obs_wave_filters, obs_trans_filters = import_data(filename)

    # Fit both models
    results_2, last_samples_2 = fit_dust_model(obs_wave, obs_flux, obs_flux_err, obs_limits, redshift, object_name,
                                               composition=composition, grain_size=grain_size, n_components=2, n_walkers=n_walkers,
                                               n_steps=n_steps, burn_in=burn_in, n_cores=n_cores, sigma_clip=sigma_clip, repeats=repeats,
                                               emcee_progress=emcee_progress, obs_wave_filters=obs_wave_filters, obs_trans_filters=obs_trans_filters,
                                               n_filter_samples=n_filter_samples, plot=plot, plots_model=plots_model, plots_corner=plots_corner,
                                               plots_trace=plots_trace, output_dir=output_dir, initial_pos=initial_pos, composition_hot=composition_hot,
                                               composition_cold=composition_cold, dust_type=dust_type, radius=radius, flat_prior=flat_prior, priors=priors)

    results_1, last_samples_1 = fit_dust_model(obs_wave, obs_flux, obs_flux_err, obs_limits, redshift, object_name,
                                               composition=composition, grain_size=grain_size, n_components=1, n_walkers=n_walkers, n_steps=n_steps,
                                               burn_in=burn_in, n_cores=n_cores, sigma_clip=sigma_clip, repeats=repeats, emcee_progress=emcee_progress,
                                               obs_wave_filters=obs_wave_filters, obs_trans_filters=obs_trans_filters, n_filter_samples=n_filter_samples,
                                               plot=plot, plots_model=plots_model, plots_corner=plots_corner, plots_trace=plots_trace, output_dir=output_dir,
                                               initial_pos=initial_pos, dust_type=dust_type, radius=radius, flat_prior=flat_prior, priors=priors)

    results = compare_models(obs_wave, obs_flux, obs_flux_err, obs_limits, redshift, object_name,
                             results_1, results_2, composition1=composition, composition2=composition, grain_size=grain_size,
                             obs_wave_filters=obs_wave_filters, obs_trans_filters=obs_trans_filters, output_dir=output_dir,
                             plot_comparison=plot, composition_hot=composition_hot, composition_cold=composition_cold,
                             dust_type=dust_type, radius=radius)

    return results
