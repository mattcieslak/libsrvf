import numpy as np
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from glob import glob

srvf_ext = Extension("libsrvf",
            sources = ["libsrvf.pyx"] + glob("src/*.cc"), 
            language="c++")

setup(
    name="srvf",
    version="0.0",
    ext_modules = cythonize(srvf_ext,language="c++"),
    include_dirs = [np.get_include(),"include"],
    
)


