from setuptools import setup, Extension
import os

prefix = os.environ.get("CONDA_PREFIX", "/usr")

effsource_circular = Extension(
    '_effsource_circular',
    sources=[
        'effsource_circular.i',
        # 'kerr-circular.c',
        'kerr-circular-ctx.c',
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
        'effsource_equatorial.i',
        'kerr-equatorial-ctx.c',
        'kerr-equatorial-coeffs-ctx.c',
        'kerr-equatorial-dtcoeffs-ctx.c',
        'kerr-equatorial-dttcoeffs-ctx.c',
    ],
    include_dirs=[prefix + '/include', '.'],
    library_dirs=[prefix + "/lib"],
    runtime_library_dirs=[prefix + "/lib"],
    libraries=['gsl', 'gslcblas', 'm'],
    extra_compile_args=['-std=gnu99', '-O0'],
)

effsource_legacy = Extension(
    '_effsource_legacy',
    sources=[
        'effsource_legacy.i',
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

effsource_legacy_elliptic = Extension(
    '_effsource_legacy_elliptic',
    sources=[
        'effsource_legacy_elliptic.i',
        'legacy_elliptic/kerr-equatorial.c',
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
    ext_modules=[effsource_circular, effsource_equatorial,
                 effsource_legacy, effsource_legacy_elliptic]
)
