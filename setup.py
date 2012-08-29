from distutils.core import setup, Extension
import numpy, glob

module1 = Extension('cchelp',
                    include_dirs = [numpy.get_include()],
                    sources = ['biopy/cchelp.cc'])

module2 = Extension('cnexus',
                    sources = ['biopy/cnexus.c'])

module3 = Extension('ctree',
                    sources = ['biopy/ctree.cc'])

classifiers=[
  "Development Status :: 3 - Alpha",
  "Environment :: Console",
  "Intended Audience :: Developers",
  "Intended Audience :: Science/Research",
  "License :: OSI Approved :: GNU Affero General Public License v3",
  "Programming Language :: Python",
  "Programming Language :: C",
  "Topic :: Scientific/Engineering :: Bio-Informatics"
  ]

from biopy import __version__

setup (name = 'biopy',
       version = __version__,
       description = 'Bioinformatics Python utilities',
       long_description = 'Bioinformatics Python utilities',
       author = 'Joseph Heled',
       author_email = 'jheled@gmail.com',
       url = 'http://http://code.google.com/p/biopy/',
       license = 'LGPL (V3)',
       classifiers = classifiers,
       platforms = ["Linux", "Mac OS-X"],
       packages = ['biopy'],
       package_dir={'biopy': 'biopy'},
       ext_modules = [module1,module2,module3],
       scripts = glob.glob('scripts/*'),
       data_files=[('doc',glob.glob('html/*.*')),
                   ('doc/_images', glob.glob('html/_images/*.*')),
                   ('doc/_images/math', glob.glob('html/_images/math/*.*')),
                   ('doc/_static', glob.glob('html/_static/*.*'))]
       )
