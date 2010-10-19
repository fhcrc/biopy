from distutils.core import setup, Extension
import numpy

module1 = Extension('cchelp',
                    include_dirs = [numpy.get_include()],
                    sources = ['biopy/cchelp.cc'])

module2 = Extension('cnexus',
                    sources = ['biopy/cnexus.c'])

setup (name = 'biopy',
       version = '1.0',
       description = 'biopy',
       packages = ['biopy'],
       package_dir={'biopy': 'biopy'},
       ext_modules = [module1,module2])
