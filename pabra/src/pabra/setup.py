'''
Created on Jun 14, 2010

@author: nbecker

setup module to compile cython parts
'''

from distutils.core import setup
from distutils.extension import Extension
import numpy
from Cython.Distutils import build_ext

#import os
#home = os.environ['HOME']  # not working anymore after $HOME change!

lib_dirs = ['/usr/lib64/']

gcy_sourcefiles = ["gillespieCy.pyx" ] # , "cgillespieCy.pxd"]
gcy_libs = ["m", "gsl", "gslcblas"]
gcy_includes = [numpy.get_include(), '.', '/usr/include/', 
                '/opt/local/include/']


#numpy.get_include(),

ext_modules = [Extension("gillespieCy", gcy_sourcefiles, libraries=gcy_libs,
                         include_dirs=gcy_includes,
#                         runtime_library_dirs=lib_dirs,
#                         , extra_compile_args=["-O3"]  # tune optimization
                         ),
               ]

print 'numpy include dirs:\n', numpy.get_include() 

print 'numpy version:', numpy.version.version
#    include_dirs = [numpy.get_include(), '.']

#print 'sys.path', sys.path

#print 'environment', os.environ

setup(
    name = 'Faster Gillespie',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)
# ('/opt/local/Library/Frameworks/Python.framework/'
#    'Versions/2.6/lib/python2.6/site-packages/'
#    'numpy/core/include')

