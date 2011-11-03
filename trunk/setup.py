from distutils.core import setup, Extension
import numpy, glob

module1 = Extension('cchelp',
                    include_dirs = [numpy.get_include()],
                    sources = ['biopy/cchelp.cc'])

module2 = Extension('cnexus',
                    sources = ['biopy/cnexus.c'])

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

setup (name = 'biopy',
       version = '0.1.7',
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
       ext_modules = [module1,module2],
       scripts = glob.glob('scripts/*'),
#       data_files=[('doc',glob.glob('html/*.html') + glob.glob('html/*.js')
#                    + glob.glob('html/_images/*') + glob.glob('html/_static/*'))]
       )
