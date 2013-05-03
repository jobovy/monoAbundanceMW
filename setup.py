from setuptools import setup #, Extension
import os, os.path
import re

longDescription= ""


setup(name='monoAbundanceMW',
      version='1.',
      description='mono-abundance MW',
      author='Jo Bovy',
      author_email='bovy@ias.edu',
      license='New BSD',
      long_description=longDescription,
      url='https://github.com/jobovy/monoAbundanceMW',
      package_dir = {'monoAbundanceMW/': ''},
      packages=['monoAbundanceMW'],
      package_data={'monoAbundanceMW': ['data/monoAbundanceResults.fits',
                                        'data/monoAbundanceResults_k.fits',]},
      dependency_links = ['https://github.com/esheldon/fitsio/tarball/master#egg=fitsio'],
      install_requires=['numpy','fitsio']
      )
