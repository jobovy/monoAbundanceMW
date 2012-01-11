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
    if numpy.sum((numpy.fabs(results['feh']-feh) <= _DFEH/2.)\
                     *(numpy.fabs(results['afe']-afe) <= _DAFE/2.)) == 0.:
        return 0.
    #Then find the relevant bin
    indx= (numpy.fabs(results['feh']-feh) <= _DFEH/2.)\
        *(numpy.fabs(results['afe']-afe) <= _DAFE/2.)
    if z is None:
        return results['mass'][indx][0]
    elif isinstance(z,(list,numpy.ndarray)):
        hz= results['hz'][indx][0]
        if z[0] is None and not z[1] is None:
            return results['mass'][indx][0]*(1.-numpy.exp(-z[1]/hz))
        elif z[1] is None and not z[0] is None:
            return results['mass'][indx][0]*numpy.exp(-z[0]/hz)
        elif z[0] is None and z[1] is None:
            return results['mass'][indx][0]
        else:
            return results['mass'][indx][0]*(numpy.exp(-z[0]/hz)-numpy.exp(-z[1]/hz))
    else:
        hz= results['hz'][indx][0]
        return results['mass'][indx][0]/2./hz*math.exp(-numpy.fabs(z)/hz)

def hz(feh,afe,err=False):
    """
    NAME:

       hz

    PURPOSE:

       return the vertical scale height as a function of feh,afe

    INPUT:

       feh - metallicity

       afe - alpha-enhancement

       err= (default: False) if True, also return error

    OUTPUT:

       vertical scale height in pc (, +uncertainty)

    HISTORY:

       2012-01-10 - Written - Bovy (IAS/@Tucson)

    """
    #First determine whether this point lies within the fit range
    if numpy.sum((numpy.fabs(results['feh']-feh) < _DFEH/2.)\
                     *(numpy.fabs(results['afe']-afe) < _DAFE/2.)) == 0.:
        return 0.
    #Then find the relevant bin
    indx= (numpy.fabs(results['feh']-feh) <= _DFEH/2.)\
        *(numpy.fabs(results['afe']-afe) <= _DAFE/2.)
    if err:
        return (results['hz'][indx][0],results['hz_err'][indx][0])
    else:
        return results['hz'][indx][0]

def hr(feh,afe,err=False):
    """
    NAME:

       hr

    PURPOSE:

       return the radial scale length as a function of feh,afe

    INPUT:

       feh - metallicity

       afe - alpha-enhancement

       err= (default: False) if True, also return error

    OUTPUT:
    
       radial scale length in kpc (, +uncertainty)

    HISTORY:

       2012-01-10 - Written - Bovy (IAS/@Tucson)

    """
    #First determine whether this point lies within the fit range
    if numpy.sum((numpy.fabs(results['feh']-feh) <= _DFEH/2.)\
                     *(numpy.fabs(results['afe']-afe) <= _DAFE/2.)) == 0.:
        return 0.
    #Then find the relevant bin
    indx= (numpy.fabs(results['feh']-feh) <= _DFEH/2.)\
        *(numpy.fabs(results['afe']-afe) <= _DAFE/2.)
    if err:
        return (results['hr'][indx][0],results['hr_err'][indx][0])
    else:
        return results['hr'][indx][0]

def sigmaz(feh,afe,z=1000.,R=8.,err=False):
    """
    NAME:

       sigmaz

    PURPOSE:

       return the vertical velocity dispersion at height Z as a function of feh,afe

    INPUT:

       feh - metallicity

       afe - alpha-enhancement

       z= (default: 1000) height [pc]

       R= (default 8) Galactocentric radius [kpc]

       err= (default: False) if True, also return error

    OUTPUT:
    
       vertical velocity dispersion at height Z

    HISTORY:

       2012-01-10 - Written - Bovy (IAS/@Tucson)

    """
    #First determine whether this point lies within the fit range
    if numpy.sum((numpy.fabs(results['feh']-feh) <= _DFEH/2.)\
                     *(numpy.fabs(results['afe']-afe) <= _DAFE/2.)) == 0.:
        return 0.
    #Then find the relevant bin
    indx= (numpy.fabs(results['feh']-feh) <= _DFEH/2.)\
        *(numpy.fabs(results['afe']-afe) <= _DAFE/2.)
    if err:
        raise NotImplementedError("Err for sigmaz not implemented yet")
        return (results['hr'][indx][0],results['hr_err'][indx][0])
    else:
        d= (z-results['zmedian'][indx])/1000.
        return (results['sz'][indx][0]+d*results['p1'][indx][0]\
            +d**2.*results['p2'][indx][0])*math.exp(-(R-8.)/results['hsz'][indx][0])

def sigmazSlope(feh,afe,err=False):
    """
    NAME:

       sigmazSlope

    PURPOSE:

       return the slope of the vertical velocity dispersion as a function of feh,afe

    INPUT:

       feh - metallicity

       afe - alpha-enhancement

       err= (default: False) if True, also return error

    OUTPUT:
    
       slope of the vertical velocity dispersion [km s^-1 kpc^-1]

    HISTORY:

       2012-01-10 - Written - Bovy (IAS/@Tucson)

    """
    #First determine whether this point lies within the fit range
    if numpy.sum((numpy.fabs(results['feh']-feh) <= _DFEH/2.)\
                     *(numpy.fabs(results['afe']-afe) <= _DAFE/2.)) == 0.:
        return 0.
    #Then find the relevant bin
    indx= (numpy.fabs(results['feh']-feh) <= _DFEH/2.)\
        *(numpy.fabs(results['afe']-afe) <= _DAFE/2.)
    if err:
        return (results['p1'][indx][0],results['p1_err'][indx][0])
    else:
        return results['p1'][indx][0]

