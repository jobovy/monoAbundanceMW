import os, os.path
import numpy
_DATADIR= 'data'
_DENSNAME= os.path.join(_DATADIR,'pixelFitG_DblExp_BigPix0.1.sav')
_VELNAME= os.path.join(_DATADIR,'pixelFitG_Vel_HWR_BigPix0.1.sav')
_MASSNAME= os.path.join(_DATADIR,'pixelFitG_Mass_DblExp_BigPix0.1_simpleage.sav')
def _load_dens_fits():
    """load the density fits"""
    
