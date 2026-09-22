from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import os

prefix = os.environ.get("CONDA_PREFIX", "/usr")

# The generated coefficient tables and the seven-channel kernel split carry
# single expressions up to 600 kB, which GCC cannot optimise in bounded memory:
# -O1 exceeds 5 GB and -O2 needs ~52 GB on kerr-equatorial-dttcoeffs-ctx.c.
# They are built at -O0; measured cost of doing so is 3.4% on set_particle.
# Everything else, notably kerr-equatorial-ctx.c which holds the pointwise
# evaluator, keeps -O2: at -O0 calc_offset costs 432 ns against 181 ns.
O0_SOURCES = {
    'kerr-equatorial-split-ctx.c',
    'kerr-equatorial-coeffs-ctx.c',
    'kerr-equatorial-dtcoeffs-ctx.c',
    'kerr-equatorial-dttcoeffs-ctx.c',
}


class build_ext_per_file_opt(build_ext):
    """Compile the listed sources at -O0, leaving every other file's flags."""

    def build_extension(self, ext):
        original = self.compiler._compile

        def patched(obj, src, ext_, cc_args, extra_postargs, pp_opts):
            args = list(extra_postargs)
            if os.path.basename(src) in O0_SOURCES:
                args = [a for a in args if not a.startswith('-O')] + ['-O0']
            return original(obj, src, ext_, cc_args, args, pp_opts)

        self.compiler._compile = patched
        try:
            super().build_extension(ext)
        finally:
            self.compiler._compile = original


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
        'kerr-equatorial-split-ctx.c',
        'kerr-equatorial-coeffs-ctx.c',
        'kerr-equatorial-dtcoeffs-ctx.c',
        'kerr-equatorial-dttcoeffs-ctx.c',
    ],
    include_dirs=[prefix + '/include', '.'],
    library_dirs=[prefix + "/lib"],
    runtime_library_dirs=[prefix + "/lib"],
    libraries=['gsl', 'gslcblas', 'm'],
    extra_compile_args=['-std=gnu99', '-O2'],
)

effsource_original = Extension(
    '_effsource_original',
    sources=[
        'effsource_original.i',
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

effsource_original_agm = Extension(
    '_effsource_original_agm',
    sources=[
        'effsource_original_agm.i',
        'original_agm/kerr-equatorial.c',
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
                 effsource_original, effsource_original_agm],
    cmdclass={'build_ext': build_ext_per_file_opt},
)
