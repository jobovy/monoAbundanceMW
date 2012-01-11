import os, os.path
import math
import numpy
import fitsio
_DATADIR= 'data'
_DATANAME= os.path.join(os.path.dirname(__file__),
                        _DATADIR,'monoAbundanceResults.fits')
#Binsize
_DFEH=0.1
_DAFE=0.05
#Load fits
results= fitsio.read(_DATANAME)
def abundanceDist(feh,afe,z=None):
    """
    NAME:

       abundanceDist

    PURPOSE:

       return the mass-weighted abundance distribution as a function of feh,afe

    INPUT:

       feh - metallicity

       afe - alpha-enhancement

       z= default=None: integrated over height; if set to number, evaluate at this height, if set to range, integrated between these numbers (can have None to indicate zero or infinity) [pc]

    OUTPUT:

       surface-mass density (or mass-density if z=number)

    HISTORY:

       2012-01-10 - Written - Bovy (IAS/@Tucson)

    """
    #First determine whether this point lies within the fit range
    if numpy.sum((numpy.fabs(results['feh']-feh) < _DFEH/2.)\
                     *(numpy.fabs(results['afe']-afe) < _DAFE/2.)) == 0.:
        return 0.
    #Then find the relevant bin
    indx= (numpy.fabs(results['feh']-feh) < _DFEH/2.)\
        *(numpy.fabs(results['afe']-afe) < _DAFE/2.)
    if z is None:
        return results['mass'][indx]
    elif isinstance(z,(list,numpy.ndarray)):
        hz= results['hz'][indx]
        if z[0] is None and not z[1] is None:
            return results['mass'][indx]*(1.-numpy.exp(-z[1]/hz))
        elif z[1] is None and not z[0] is None:
            return results['mass'][indx]*numpy.exp(-z[0]/hz)
        elif z[0] is None and z[1] is None:
            return results['mass'][indx]
        else:
            return results['mass'][indx]*(numpy.exp(-z[0]/hz)-numpy.exp(-z[1]/hz))
    else:
        hz= results['hz'][indx]
        return results['mass'][indx]/2./hz*math.exp(-numpy.fabs(z)/hz)
