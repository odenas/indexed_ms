
import glob
from setuptools import setup, Extension
from Cython.Build import cythonize

extensions = [
    Extension(
        'mstat.ms',
        ['src/ms.pyx'],
        include_dirs=["../sdsl-lite/build/include", "../fast_ms/src"],
        libraries=["sdsl", "divsufsort", "divsufsort64", "dl", "pthread"],
        library_dirs=["../sdsl-lite/build/lib"],
    )
]

setup(
        name='mstat',
        description='Python interface to matching statistics',
        author='Olgert Denas',
        author_email='gertidenas@gmail.com',
        version=0.1,
        packages=['mstat'],
        package_dir={'mstat': 'src'},
        scripts=glob.glob('bin/*py'),
        ext_modules=cythonize(extensions, language="c++"),
        #ext_modules=[Extension("dimer.genome.bed", ["lib/genome/bed.pyx"]),
        #              Extension("dimer.genome.peak", ["lib/genome/peak.pyx"])],
        install_requires=['numpy > 1.5']
)
