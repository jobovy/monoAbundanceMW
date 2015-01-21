monoAbundanceMW
-----------------

**Results from the Bovy et al. (2012) mono-abundance papers**

DESCRIPTION
============

This module provides access to the fit results from `Bovy et
al. (2012a) <http://adsabs.harvard.edu/abs/2012ApJ...753..148B>`__,
`Bovy, Rix, & Hogg (2012)
<http://adsabs.harvard.edu/abs/2012ApJ...751..131B>`__, and `Bovy et
al. (2012b)
<http://adsabs.harvard.edu/abs/2012ApJ...755..115B>`__. They can
either be used through a python interface or are directly accessible
through a FITS file (see below).

AUTHOR
======

Jo Bovy - bovy at ias dot edu

INSTALLATION
============

Standard python setup.py build/install

Either

``sudo python setup.py install``

or 

``python setup.py install --prefix=/some/directory/``

USE
===

Python interface
+++++++++++++++++

A Python interface to the results is available as::

  import monoAbundanceMW as map
  map.hz(-0.5,0.25)
  443.2816423047887 #pc

for example. See the ``__init__.py`` file in ``monoAbundanceMW/`` for
more info on various functions.


FITS file
++++++++++

The results are also available in `this FITS file
<https://github.com/jobovy/monoAbundanceMW/blob/master/monoAbundanceMW/data/monoAbundanceResults.fits>`__,
which you can download by clicking on ``raw``. This file has the fit
results for each bin, and contains the [Fe/H], [a/Fe] that the bin is
centered on, and all of the parameters describing the density and
velocity dispersion. The parameter names are pretty self-explanatory,
but here's a quick overview::

		  feh: metallicity
		  afe: alpha-enhancement
		  hz: scale height in pc
		  hr: scale length in kpc
		  hz_err, hr_err: errors on the two above
		  mass: total surface density in Msolar/pc^2
		  zmedian: median height of the data sample in pc (this is the pivot for
		  the kinematics fits)
		  sz: vertical dispersion in km/s at zmedian
		  p1: slope of the vertical dispersion in km/s/kpc
		  p2: second derivative of the vertical dispersion in km/s/kpc/kpc
		  hsz: radial scale length of the vertical velocity dispersion in kpc
		  sz_err, p1_err, p2_err, hsz_err: errors on the previous four
		  szp1_corr, szp2_corr, ...: correlation coefficients for the four
		  kinematics parameters
		  sr: radial dispersion in km/s at zmedian
		  hsr: radial scale length of the radial velocity dispersion in kpc
		  sr_err: error on sr

[For the radial velocity dispersion only the results from a fit with constant dispersion(height) are included, so the p1 and p2 parameters for sr are set to zero].

Please contact the author for more information.