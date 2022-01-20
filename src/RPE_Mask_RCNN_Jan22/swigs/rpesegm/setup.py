#! /usr/bin/env python

# System imports
import platform
import setuptools
from distutils.core import *
from distutils      import sysconfig

# Third-party modules - we depend on numpy for everything
import numpy

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

if platform.system().lower() == 'windows':
    _complile_args = [] # ["/EHsc"]
else:
    _complile_args = ["-std=c++11", "-Wno-unused-variable", "-Wno-unused-but-set-variable", "-Wno-maybe-uninitialized"]

# rpesegm extension module
_rpesegm = Extension("_rpesegm",
                   ["rpesegm.i", "rpesegm.cpp", "assembly3d.cpp", "geom.cpp", "raster.cpp", "csv.cpp"],
                   include_dirs = [numpy_include],
                   swig_opts = ["-c++"],
                   extra_compile_args = _complile_args,
                   )

# rpeutil extension module
_rpeutil = Extension("_rpeutil",
                   ["rpeutil.i", "rpeutil.cpp", "flatcont.cpp", "geom.cpp", "raster.cpp", "csv.cpp"],
                   include_dirs = [numpy_include],
                   swig_opts = ["-c++"],
                   extra_compile_args = _complile_args,
                   )

# rpesegm setup
setup(  name        = "RPE Segmentation",
        description = "Native support for RPE Map Segmentation",
        author      = "Andrei Volkov",
        version     = "1.0",
        ext_modules = [_rpesegm, _rpeutil],
        )
