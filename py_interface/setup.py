
import glob
from setuptools import *
import Cython.Distutils
command_classes = {'build_ext': Cython.Distutils.build_ext}

setup(
        name='mstat',
        description='Python interface to matching statistics',
        author='Olgert Denas',
        author_email='gertidenas@gmail.com',
        version=0.1,
        packages=['mstat'],
        package_dir={'mstat': 'src'},
        scripts=glob.glob('bin/*py'),
        ext_modules=[],
        #ext_modules=[Extension("dimer.genome.bed", ["lib/genome/bed.pyx"]),
        #              Extension("dimer.genome.peak", ["lib/genome/peak.pyx"])],
        install_requires=['Cython >= 0.27', 'numpy > 1.5'],
        cmdclass=command_classes
        )



