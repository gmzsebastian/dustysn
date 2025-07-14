.. _reference:

Reference
=========

Dust
----

The functions in ``dustysn`` are based off the mass absorption coefficients
from `Sarangi et al. 2022 <https://ui.adsabs.harvard.edu/abs/2022A%26A...668A..57S/abstract>`__.
Currently, the supported compositions are carbon and silicate dust. The carbon dust can have
grain sizes of 0.01, 0.1, or 1.0 microns, while the silicate dust can have grain sizes of 0.1 microns.

Models
------

The flux calculation in this model is based on the intrinsic luminosity of the dust emission and the observed flux density. First, the intrinsic luminosity of the dust is calculated from the blackbody emission model, given the temperature and dust mass.

.. math::

   L(\nu) = M_{\rm dust} \kappa(\nu) B_{\nu}(T) 4 \pi

where:

* :math:`L(\nu)` is the intrinsic luminosity in erg/s/Hz or erg/s/Å
* :math:`M_{\rm dust}` is the dust mass in solar masses
* :math:`\kappa(\nu)` is the mass absorption coefficient in cm²/g, which is interpolated to the rest-frame wavelength
* :math:`B_{\nu}(T)` is the Planck function for blackbody radiation at a given temperature, which is calculated from the temperature of the dust (in Kelvin)

The Planck function :math:`B_{\nu}(T)` is computed using the ``BlackBody`` Astropy model.

Once the luminosity is computed, the flux density :math:`f_{\rm obs}(\nu)` is calculated from the luminosity and the luminosity distance :math:`d` to the object, with an adjustment for redshift. The flux is given by:

.. math::

   f_{\rm obs}(\nu) = \frac{L(\nu)}{4 \pi d^2} (1 + z)

where:

* :math:`f_{\rm obs}(\nu)` is the observed flux in Jy (Jansky), erg/s/cm²/Hz, or erg/s/cm²/Å
* :math:`L(\nu)` is the intrinsic luminosity at the observed wavelength
* :math:`d` is the luminosity distance in cm
* :math:`z` is the redshift of the object

Filters
-------

Given that this package was designed with the main goal of fitting JWST photometry, we include information for all MIRI and NIRCam filters.
The filter passbands were obtained from CRDS based on the data from 2025 Jan 15.
