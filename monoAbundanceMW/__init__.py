import os, os.path
import math
import numpy
import fitsio
try:
    from galpy.util import bovy_plot
except ImportError:
    _GALPYLOADED= False
else:
    _GALPYLOADED= True
_DATADIR= 'data'
_DATANAME= os.path.join(os.path.dirname(__file__),
                        _DATADIR,'monoAbundanceResults.fits')
_DATANAME_K= os.path.join(os.path.dirname(__file__),
                          _DATADIR,'monoAbundanceResults_k.fits')
#Binsize
_DFEH=0.1
_DAFE=0.05
#Load fits
results= fitsio.read(_DATANAME)
results_k= fitsio.read(_DATANAME_K)
def afes(k=False):
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
    if k:
        return results_k['afe']
    else:
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

       d<[Fe/H]>/dR(z) in dex/kpc

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

def dfehdz(z=1000.,r=None,feh=None,afe=None):
    """
    NAME:

       dfehdz

    PURPOSE:

       return the vertical [Fe/H] gradient at (R,z)

    INPUT:

       z= height in pc

       r= Galactocentric distance (kpc)

       feh= (default: None) range in feh to consider

       afe= (default: None) range in afe to consider

    OUTPUT:

       d<[Fe/H]>/dZ(R,z) in dex/kpc

    HISTORY:

       2013-08-02 - Written - Bovy (IAS)

    """
    #First get the weights and hzs
    w= numpy.zeros(len(results['afe']))
    hzs= numpy.zeros_like(w)
    for ii in range(len(results['afe'])):
        w[ii]= abundanceDist(results['feh'][ii],results['afe'][ii],
                             z=z,r=r,number=False)
        hzs[ii]= hz(results['feh'][ii],results['afe'][ii])/1000.
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
    return -numpy.sum(w[indx]/hzs[indx]*results['feh'][indx])\
        /numpy.sum(w[indx])\
        +numpy.sum(w[indx]*results['feh'][indx])\
        /numpy.sum(w[indx])**2.\
        *numpy.sum(w[indx]/hzs[indx])

def fehs(k=False):
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
    if k:
        return results_k['feh']
    else:
        return results['feh']

def hr(feh,afe,err=False,k=False):
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
    if k:
        tresults= results_k
    else:
        tresults= results
    #First determine whether this point lies within the fit range
    if numpy.sum((numpy.fabs(tresults['feh']-feh) <= _DFEH/2.)\
                     *(numpy.fabs(tresults['afe']-afe) <= _DAFE/2.)) == 0.:
        return numpy.nan
    #Then find the relevant bin
    indx= (numpy.fabs(tresults['feh']-feh) <= _DFEH/2.)\
        *(numpy.fabs(tresults['afe']-afe) <= _DAFE/2.)
    if err:
        return (tresults['hr'][indx][0],tresults['hr_err'][indx][0])
    else:
        return tresults['hr'][indx][0]

def hz(feh,afe,err=False,k=False):
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
    if k:
        tresults= results_k
    else:
        tresults= results
    #First determine whether this point lies within the fit range
    if numpy.sum((numpy.fabs(tresults['feh']-feh) <= _DFEH/2.)\
                     *(numpy.fabs(tresults['afe']-afe) <= _DAFE/2.)) == 0.:
        return numpy.nan
    #Then find the relevant bin
    indx= (numpy.fabs(tresults['feh']-feh) <= _DFEH/2.)\
        *(numpy.fabs(tresults['afe']-afe) <= _DAFE/2.)
    if err:
        return (tresults['hz'][indx][0],tresults['hz_err'][indx][0])
    else:
        return tresults['hz'][indx][0]

def meanfeh(z=1000.,r=None,indx=None):
    """
    NAME:

       meanfeh

    PURPOSE:

       return the mean [Fe/H] as a function of height z

    INPUT:

       z= height from the plane (pc)

       r= Galactocentric distance (kpc)

       indx= if given, use only these MAPs

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
    if not indx is None:
        w[True-indx]= 0.
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

def sigmaz(feh,afe,z=None,r=None,err=False,smooth=False,smoothfunc=numpy.mean,
           k=False):
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
       
       smooth= (default: False) if True, return a smoothed estimate (using the nine bins around the requested bin, including the input FeH and aFe) (does not work with err=True)

       smoothfunc= function to use for smoothing (default: numpy.mean)

    OUTPUT:
    
       vertical velocity dispersion at height Z

    HISTORY:

       2012-01-10 - Written - Bovy (IAS/@Tucson)

       2012-05-03 - change default z - Bovy (IAS)

    """
    if smooth and not err:
        #Smooth sz                                                      
        up= sigmaz(feh+0.1,afe,z=z,r=r,smooth=False)
        down= sigmaz(feh-0.1,afe,z=z,r=r,smooth=False)
        left= sigmaz(feh,afe-0.05,z=z,r=r,smooth=False)
        right= sigmaz(feh,afe+0.05,z=z,r=r,smooth=False)
        here= sigmaz(feh,afe,z=z,r=r,smooth=False)
        upright= sigmaz(feh+0.1,afe+0.05,z=z,r=r,smooth=False)
        upleft= sigmaz(feh+0.1,afe-0.05,z=z,r=r,smooth=False)
        downright= sigmaz(feh-0.1,afe+0.05,z=z,r=r,smooth=False)
        downleft= sigmaz(feh-0.1,afe-0.05,z=z,r=r,smooth=False)
        allsz= numpy.array([here,up,down,left,right,upright,upleft,downright,
                            downleft])
        indx= True-numpy.isnan(allsz)
        return smoothfunc(allsz[True-numpy.isnan(allsz)])
    elif smooth and err:
        raise NotImplementedError("sigmaz with smooth=True and err=True not implemented yet")
    if k:
        tresults= results_k
    else:
        tresults= results
    if r is None:
        r= 8.
    #First determine whether this point lies within the fit range
    if numpy.sum((numpy.fabs(tresults['feh']-feh) <= _DFEH/2.)\
                     *(numpy.fabs(tresults['afe']-afe) <= _DAFE/2.)) == 0.:
        return numpy.nan
    #Then find the relevant bin
    indx= (numpy.fabs(tresults['feh']-feh) <= _DFEH/2.)\
        *(numpy.fabs(tresults['afe']-afe) <= _DAFE/2.)
    if z is None:
        z= tresults['zmedian'][indx][0]
    if err:
        if z != tresults['zmedian'][indx][0]:
            raise NotImplementedError("Err for sigmaz not implemented for z =/= zmedian")
        return (tresults['sz'][indx][0],tresults['sz_err'][indx][0])
    else:
        d= (z-tresults['zmedian'][indx][0])/1000.
        return (tresults['sz'][indx][0]+d*tresults['p1'][indx][0]\
            +d**2.*tresults['p2'][indx][0])*math.exp(-(r-8.)/tresults['hsz'][indx][0])

def sigmar(feh,afe,z=None,r=None,err=False,smooth=False,smoothfunc=numpy.mean,
           k=False):
    """
    NAME:

       sigmar

    PURPOSE:

       return the radial velocity dispersion at height Z as a function of feh,afe

    INPUT:

       feh - metallicity

       afe - alpha-enhancement

       z= (default: median height of data sample) height [pc]

       r= (default 8) Galactocentric radius [kpc]

       err= (default: False) if True, also return error

       smooth= (default: False) if True, return a smoothed estimate (using the nine bins around the requested bin, including the input FeH and aFe) (does not work with err=True)

       smoothfunc= function to use for smoothing (default: numpy.mean)

    OUTPUT:
    
       radial velocity dispersion at height Z

    BUGS:

       needs to be updated to use the slope and quadratic coefficient measurements

    HISTORY:

       2013-03-08 - Written - Bovy (IAS)

    """
    if smooth and not err:
        #Smooth sr
        up= sigmar(feh+0.1,afe,z=z,r=r,smooth=False)
        down= sigmar(feh-0.1,afe,z=z,r=r,smooth=False)
        left= sigmar(feh,afe-0.05,z=z,r=r,smooth=False)
        right= sigmar(feh,afe+0.05,z=z,r=r,smooth=False)
        here= sigmar(feh,afe,z=z,r=r,smooth=False)
        upright= sigmar(feh+0.1,afe+0.05,z=z,r=r,smooth=False)
        upleft= sigmar(feh+0.1,afe-0.05,z=z,r=r,smooth=False)
        downright= sigmar(feh-0.1,afe+0.05,z=z,r=r,smooth=False)
        downleft= sigmar(feh-0.1,afe-0.05,z=z,r=r,smooth=False)
        allsr= numpy.array([here,up,down,left,right,upright,upleft,downright,
                            downleft])
        indx= True-numpy.isnan(allsr)
        return smoothfunc(allsr[True-numpy.isnan(allsr)])
    elif smooth and err:
        raise NotImplementedError("sigmaz with smooth=True and err=True not implemented yet")
    if r is None:
        r= 8.
    if k:
        tresults= results_k
    else:
        tresults= results
    #First determine whether this point lies within the fit range
    if numpy.sum((numpy.fabs(tresults['feh']-feh) <= _DFEH/2.)\
                     *(numpy.fabs(tresults['afe']-afe) <= _DAFE/2.)) == 0.:
        return numpy.nan
    #Then find the relevant bin
    indx= (numpy.fabs(tresults['feh']-feh) <= _DFEH/2.)\
        *(numpy.fabs(tresults['afe']-afe) <= _DAFE/2.)
    if z is None:
        z= tresults['zmedian'][indx][0]
    if err:
        if z != tresults['zmedian'][indx][0]:
            raise NotImplementedError("Err for sigmaz not implemented for z =/= zmedian")
        return (tresults['sr'][indx][0],tresults['sr_err'][indx][0])
    else:
        d= (z-tresults['zmedian'][indx][0])/1000.
        return (tresults['sr'][indx][0]+0.*d*tresults['p1'][indx][0]\
                    +0.*d**2.*tresults['p2'][indx][0])*math.exp(-(r-8.)/tresults['hsr'][indx][0])

def sigmazSlope(feh,afe,err=False,k=False):
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
    if k:
        tresults= results_k
    else:
        tresults= results
    #First determine whether this point lies within the fit range
    if numpy.sum((numpy.fabs(tresults['feh']-feh) <= _DFEH/2.)\
                     *(numpy.fabs(tresults['afe']-afe) <= _DAFE/2.)) == 0.:
        return numpy.nan
    #Then find the relevant bin
    indx= (numpy.fabs(tresults['feh']-feh) <= _DFEH/2.)\
        *(numpy.fabs(tresults['afe']-afe) <= _DAFE/2.)
    if err:
        return (tresults['p1'][indx][0],tresults['p1_err'][indx][0])
    else:
        return tresults['p1'][indx][0]

def sigmazCurv(feh,afe,err=False,k=False):
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
    if k:
        tresults= results_k
    else:
        tresults= results
    #First determine whether this point lies within the fit range
    if numpy.sum((numpy.fabs(tresults['feh']-feh) <= _DFEH/2.)\
                     *(numpy.fabs(tresults['afe']-afe) <= _DAFE/2.)) == 0.:
        return numpy.nan
    #Then find the relevant bin
    indx= (numpy.fabs(tresults['feh']-feh) <= _DFEH/2.)\
        *(numpy.fabs(tresults['afe']-afe) <= _DAFE/2.)
    if err:
        return (tresults['p2'][indx][0],tresults['p2_err'][indx][0])
    else:
        return tresults['p2'][indx][0]

def plotPixelFunc(feh,afe,z,
                  fehmin=-1.6,fehmax=0.5,afemin=-0.05,afemax=0.55,
                  **kwargs):
    """
    NAME:

       plotPixelFunc

    PURPOSE:

       Plot a function as a function of [Fe/H] and [a/Fe], a la Bovy et 
       al. (2012)

    INPUT:

       feh - feh values of z

       afe - afe values of z

       z - function to plot (array with the same length as feh and afe)

       bovy_plot.bovy_dens2d kwargs

    OUTPUT:

       plot to output device

    HISTORY:
       2013-03-19 - Written - Bovy (IAS)
    """
    if not _GALPYLOADED:
        raise ImportError("Use of plotPixelFunc requires galpy's bovy_plot to be installed; install galpy starting at https://github.com/jobovy/galpy")
    #First create 2D
    gfeh= numpy.linspace(fehmin+0.05,fehmax-0.05,int((fehmax-fehmin)/0.1))
    gafe= numpy.linspace(afemin+0.025,afemax-0.025,int((afemax-afemin)/0.05))
    z2d= numpy.empty((int((fehmax-fehmin)/0.1),int((afemax-afemin)/0.05)))
    for ii in range(z2d.shape[0]):
        for jj in range(z2d.shape[1]):
            if numpy.amin((gfeh[ii]-feh)**2./0.01**2.
                          +(gafe[jj]-afe)**2./0.05**2.) < 10.**-4.:
                indx= numpy.argmin((gfeh[ii]-feh)**2./0.01**2.
                                   +(gafe[jj]-afe)**2./0.05**2.)
                z2d[ii,jj]= z[indx]
            else:
                z2d[ii,jj]= numpy.nan
    #Now plot
    xrange=[-1.6,0.5]
    yrange=[-0.05,0.55]
    if not kwargs.has_key('colorbar'):
        kwargs['colorbar']= True
    if not kwargs.has_key('shrink'):
        kwargs['shrink']= 0.78
    if not kwargs.has_key('vmin'):
        kwargs['vmin']= numpy.nanmin(z2d)
    if not kwargs.has_key('vmax'):
        kwargs['vmax']= numpy.nanmax(z2d)
    return bovy_plot.bovy_dens2d(z2d.T,origin='lower',cmap='jet',
                                 interpolation='nearest',
                                 xlabel=r'$[\mathrm{Fe/H}]$',
                                 ylabel=r'$[\alpha/\mathrm{Fe}]$',
                                 xrange=xrange,yrange=yrange,
                                 contours=False,
                                 **kwargs)

def _mgfeh(feh):
    return 0.956+0.205*feh+0.051*feh**2.
def _omegafeh(feh):
    return 0.0425+0.0198*feh+0.0057*feh**2.
