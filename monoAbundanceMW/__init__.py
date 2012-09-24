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
def afes():
    """
    NAME:

       afes

    PURPOSE:

       return the [alpha/Fe] for which we have measurements

    INPUT:

       (none)

    OUTPUT:

       [alpha/Fe] for which we have measurements

    HISTORY:

       2012-01-11 - Written - Bovy (IAS/@Tucson)

    """
    return results['afe']

def abundanceDist(feh,afe,z=None,r=None,number=False):
    """
    NAME:

       abundanceDist

    PURPOSE:

       return the mass-weighted abundance distribution as a function of feh,afe

    INPUT:

       feh - metallicity

       afe - alpha-enhancement

       z= default=None: integrated over height; if set to number, evaluate at this height, if set to range, integrated between these numbers (can have None to indicate zero or infinity) [pc]

       r= default=None: Galactocentric radius in kpc

       number= (default: False) if True, return the number density of G-type dwarfs

    OUTPUT:

       surface-mass density (or mass-density if z=number)

    HISTORY:

       2012-01-10 - Written - Bovy (IAS/@Tucson)

    """
    #First determine whether this point lies within the fit range
    if numpy.sum((numpy.fabs(results['feh']-feh) <= _DFEH/2.)\
                     *(numpy.fabs(results['afe']-afe) <= _DAFE/2.)) == 0.:
        return numpy.nan
    #Then find the relevant bin
    indx= (numpy.fabs(results['feh']-feh) <= _DFEH/2.)\
        *(numpy.fabs(results['afe']-afe) <= _DAFE/2.)
    if not r is None:
        rFactor= numpy.exp(-(r-8.)/hr(results['feh'][indx],results['afe'][indx]))
    else: rFactor= 1.
    if number:
        numberFactor= _omegafeh(feh)/_mgfeh(feh)
    else:
        numberFactor= 1.
    if z is None:
        return results['mass'][indx][0]*numberFactor*rFactor
    elif isinstance(z,(list,numpy.ndarray)):
        hz= results['hz'][indx][0]
        if z[0] is None and not z[1] is None:
            return results['mass'][indx][0]*(1.-numpy.exp(-z[1]/hz))*numberFactor*rFactor
        elif z[1] is None and not z[0] is None:
            return results['mass'][indx][0]*numpy.exp(-z[0]/hz)*numberFactor*rFactor
        elif z[0] is None and z[1] is None:
            return results['mass'][indx][0]*numberFactor*rFactor
        else:
            return results['mass'][indx][0]*(numpy.exp(-z[0]/hz)-numpy.exp(-z[1]/hz))*numberFactor*rFactor
    else:
        hz= results['hz'][indx][0]
        return results['mass'][indx][0]/2./hz*math.exp(-numpy.fabs(z)/hz)*numberFactor*rFactor

def dfehdr(z=1000.,r=None,feh=None,afe=None):
    """
    NAME:

       dfehdr

    PURPOSE:

       return the radial [Fe/H] gradient at height z

    INPUT:

       z= height in pc

       r= Galactocentric distance (kpc)

       feh= (default: None) range in feh to consider

       afe= (default: None) range in afe to consider

    OUTPUT:

       d<[Fe/H]>/dR(z)

    HISTORY:

       2012-07-19 - Written - Bovy (IAS/@MPIA)

    """
    #First get the weights and hrs
    w= numpy.zeros(len(results['afe']))
    hrs= numpy.zeros_like(w)
    for ii in range(len(results['afe'])):
        w[ii]= abundanceDist(results['feh'][ii],results['afe'][ii],
                             z=z,r=r,number=False)
        hrs[ii]= hr(results['feh'][ii],results['afe'][ii])
    if False:
        #Cut out pops with undetermined scale lengths (very little mass)
        indx= (hrs < 4.5)
        #indx= (results['afe'] >= 0.05)*(results['afe'] < 0.1)
    else:
        indx= numpy.zeros(len(w),dtype='bool')+True
    if not feh is None:
        indx= indx*(results['feh'] > feh[0])*(results['feh'] < feh[1])
    if not afe is None:
        indx= indx*(results['afe'] > afe[0])*(results['afe'] < afe[1])
    #Then return
    return -numpy.sum(w[indx]/hrs[indx]*results['feh'][indx])\
        /numpy.sum(w[indx])\
        +numpy.sum(w[indx]*results['feh'][indx])\
        /numpy.sum(w[indx])**2.\
        *numpy.sum(w[indx]/hrs[indx])

def fehs():
    """
    NAME:

       fehs

    PURPOSE:

       return the [Fe/H] for which we have measurements

    INPUT:

       (none)

    OUTPUT:

       [Fe/H] for which we have measurements

    HISTORY:

       2012-01-11 - Written - Bovy (IAS/@Tucson)

    """
    return results['feh']

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
        return numpy.nan
    #Then find the relevant bin
    indx= (numpy.fabs(results['feh']-feh) <= _DFEH/2.)\
        *(numpy.fabs(results['afe']-afe) <= _DAFE/2.)
    if err:
        return (results['hr'][indx][0],results['hr_err'][indx][0])
    else:
        return results['hr'][indx][0]

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
    if numpy.sum((numpy.fabs(results['feh']-feh) <= _DFEH/2.)\
                     *(numpy.fabs(results['afe']-afe) <= _DAFE/2.)) == 0.:
        return numpy.nan
    #Then find the relevant bin
    indx= (numpy.fabs(results['feh']-feh) <= _DFEH/2.)\
        *(numpy.fabs(results['afe']-afe) <= _DAFE/2.)
    if err:
        return (results['hz'][indx][0],results['hz_err'][indx][0])
    else:
        return results['hz'][indx][0]

def meanfeh(z=1000.,r=None):
    """
    NAME:

       meanfeh

    PURPOSE:

       return the mean [Fe/H] as a function of height z

    INPUT:

       z= height from the plane (pc)

       r= Galactocentric distance (kpc)

    OUTPUT:
    
       mean [Fe/H]at height z

    HISTORY:

       2012-07-19 - Written - Bovy (IAS@MPIA)

    """
    #First get the weights and hrs
    w= numpy.zeros(len(results['afe']))
    for ii in range(len(results['afe'])):
        w[ii]= abundanceDist(results['feh'][ii],results['afe'][ii],
                             z=z,r=r,number=False)
    if False:
        #Cut out pops with undetermined scale lengths (very little mass)
        #indx= (hrs < 4.5)
        indx= (results['afe'] > 0.2)#hrs < 4.5)
    else:
        indx= numpy.zeros(len(w),dtype='bool')+True
    #Then return
    return numpy.sum(w[indx]*results['feh'][indx])/numpy.sum(w[indx])

def medianfeh(z=1000.,r=None):
    """
    NAME:

       medianfeh

    PURPOSE:

       return the median [Fe/H] as a function of height z

    INPUT:

       z= height from the plane (pc)
       
       r= Galactocentric distance (kpc)

    OUTPUT:
    
       median [Fe/H]at height z and radius r

    HISTORY:

       2012-07-19 - Written - Bovy (IAS@MPIA)

    """
    #First get the weights and hrs
    w= numpy.zeros(len(results['afe']))
    for ii in range(len(results['afe'])):
        w[ii]= abundanceDist(results['feh'][ii],results['afe'][ii],
                             z=z,r=r,number=False)
    if False:
        #Cut out pops with undetermined scale lengths (very little mass)
        indx= (hrs < 4.5)
    else:
        indx= numpy.zeros(len(w),dtype='bool')+True
    w= w[indx]
    feh= results['feh'][indx]
    sortindx= numpy.argsort(feh)
    feh= feh[sortindx]
    w= w[sortindx]
    w/= numpy.sum(w)
    tot, ii= 0., 0
    while tot < 0.5:
        tot+= w[ii]
        ii+= 1
    ii-= 1
    #Then return
    print w, feh
    if feh[ii] == feh[ii+1]:
        return feh[ii]
    else:
        nthisfeh= numpy.sum((feh == feh[ii]))
        return feh[ii]-(tot-0.5)/w[ii]*_DFEH/nthisfeh+_DFEH/2./nthisfeh #NOT SURE THAT THIS WORKS

def meanhr(z=1000.):
    """
    NAME:

       meanhr

    PURPOSE:

       return the mean radial scale length as a function of height z

    INPUT:

       z= height from the plane (pc)

    OUTPUT:
    
       mean radial scale length in kpc at height z

    HISTORY:

       2012-05-25 - Written - Bovy (IAS)

    """
    #First get the weights and hrs
    w= numpy.zeros(len(results['afe']))
    hrs= numpy.zeros_like(w)
    for ii in range(len(results['afe'])):
        w[ii]= abundanceDist(results['feh'][ii],results['afe'][ii],
                             z=z,number=False)
        hrs[ii]= hr(results['feh'][ii],results['afe'][ii])
    #Cut out pops with undetermined scale lengths (very little mass)
    indx= (hrs < 4.5)
    #Then return
    return numpy.sum(w[indx]*hrs[indx])/numpy.sum(w[indx])

def meanhz(z=1000.):
    """
    NAME:

       meanhz

    PURPOSE:

       return the mean vertical scale height as a function of height z

    INPUT:

       z= height from the plane (pc)

    OUTPUT:
    
       mean vertical scale height in pc at height z

    HISTORY:

       2012-05-25 - Written - Bovy (IAS)

    """
    #First get the weights and hrs
    w= numpy.zeros(len(results['afe']))
    hzs= numpy.zeros_like(w)
    for ii in range(len(results['afe'])):
        w[ii]= abundanceDist(results['feh'][ii],results['afe'][ii],
                             z=z,number=False)
        hzs[ii]= hz(results['feh'][ii],results['afe'][ii])
    #Then return
    return numpy.sum(w*hzs)/numpy.sum(w)

def meansigmaz(z=1000.,r=None):
    """
    NAME:

       meansigmaz

    PURPOSE:

       return the mean vertical velocity dispersion at height Z

    INPUT:

       z= (default: 1000) height [pc]

       r= Galactocentric distance (kpc)

    OUTPUT:
    
       mean vertical velocity dispersion at height Z

    HISTORY:

       2012-05-25 - Written - Bovy (IAS/@Tucson)

    """
    #First get the weights and sigmazs
    w= numpy.zeros(len(results['afe']))
    sz= numpy.zeros_like(w)
    for ii in range(len(results['afe'])):
        w[ii]= abundanceDist(results['feh'][ii],results['afe'][ii],
                             z=z,r=r,number=False)
        sz[ii]= sigmaz(results['feh'][ii],results['afe'][ii],
                       z=z,r=r)
    #Then return
    return numpy.sqrt(numpy.sum(w*sz**2.)/numpy.sum(w))

def sigmaz(feh,afe,z=None,r=None,err=False):
    """
    NAME:

       sigmaz

    PURPOSE:

       return the vertical velocity dispersion at height Z as a function of feh,afe

    INPUT:

       feh - metallicity

       afe - alpha-enhancement

       z= (default: median height of data sample) height [pc]

       r= (default 8) Galactocentric radius [kpc]

       err= (default: False) if True, also return error

    OUTPUT:
    
       vertical velocity dispersion at height Z

    HISTORY:

       2012-01-10 - Written - Bovy (IAS/@Tucson)

       2012-05-03 - change default z - Bovy (IAS)

    """
    if r is None:
        r= 8.
    #First determine whether this point lies within the fit range
    if numpy.sum((numpy.fabs(results['feh']-feh) <= _DFEH/2.)\
                     *(numpy.fabs(results['afe']-afe) <= _DAFE/2.)) == 0.:
        return numpy.nan
    #Then find the relevant bin
    indx= (numpy.fabs(results['feh']-feh) <= _DFEH/2.)\
        *(numpy.fabs(results['afe']-afe) <= _DAFE/2.)
    if z is None:
        z= results['zmedian'][indx][0]
    if err:
        if z != results['zmedian'][indx][0]:
            raise NotImplementedError("Err for sigmaz not implemented for z =/= zmedian")
        return (results['sz'][indx][0],results['sz_err'][indx][0])
    else:
        d= (z-results['zmedian'][indx][0])/1000.
        return (results['sz'][indx][0]+d*results['p1'][indx][0]\
            +d**2.*results['p2'][indx][0])*math.exp(-(r-8.)/results['hsz'][indx][0])

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
        return numpy.nan
    #Then find the relevant bin
    indx= (numpy.fabs(results['feh']-feh) <= _DFEH/2.)\
        *(numpy.fabs(results['afe']-afe) <= _DAFE/2.)
    if err:
        return (results['p1'][indx][0],results['p1_err'][indx][0])
    else:
        return results['p1'][indx][0]

def sigmazCurv(feh,afe,err=False):
    """
    NAME:

       sigmazCurv

    PURPOSE:

       return the second Z-derivative of the vertical velocity dispersion as a function of feh,afe

    INPUT:

       feh - metallicity

       afe - alpha-enhancement

       err= (default: False) if True, also return error

    OUTPUT:
    
       second derivative of the vertical velocity dispersion [km s^-1 kpc^-2]

    HISTORY:

       2012-05-03 - Written - Bovy (IAS)

    """
    #First determine whether this point lies within the fit range
    if numpy.sum((numpy.fabs(results['feh']-feh) <= _DFEH/2.)\
                     *(numpy.fabs(results['afe']-afe) <= _DAFE/2.)) == 0.:
        return numpy.nan
    #Then find the relevant bin
    indx= (numpy.fabs(results['feh']-feh) <= _DFEH/2.)\
        *(numpy.fabs(results['afe']-afe) <= _DAFE/2.)
    if err:
        return (results['p2'][indx][0],results['p2_err'][indx][0])
    else:
        return results['p2'][indx][0]

def _mgfeh(feh):
    return 0.956+0.205*feh+0.051*feh**2.
def _omegafeh(feh):
    return 0.0425+0.0198*feh+0.0057*feh**2.
