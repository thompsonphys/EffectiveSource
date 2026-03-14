from setuptools import setup, Extension
import os

prefix = os.environ.get("CONDA_PREFIX", "/usr")

effsource_circular = Extension(
    '_effsource_circular',
    sources=[
        'effsource_circular_wrap.c',
        'kerr-circular.c',
    ],
    include_dirs=[prefix + '/include', '.'],
    library_dirs=[prefix + "/lib"],
    runtime_library_dirs=[prefix + "/lib"],
    libraries=['gsl', 'gslcblas', 'm'],
    extra_compile_args=['-std=gnu99', '-O3'],
)

effsource_equatorial = Extension(
    '_effsource_equatorial',
    sources=[
        'effsource_equatorial_wrap.c',
        'kerr-equatorial.c',
        'kerr-equatorial-coeffs.c',
        'kerr-equatorial-dtcoeffs.c',
        'kerr-equatorial-dttcoeffs.c',
    ],
    include_dirs=[prefix + '/include', '.'],
    library_dirs=[prefix + "/lib"],
    runtime_library_dirs=[prefix + "/lib"],
    libraries=['gsl', 'gslcblas', 'm'],
    extra_compile_args=['-std=gnu99', '-O0'],
)

setup(
    ext_modules=[effsource_circular, effsource_equatorial]
)
