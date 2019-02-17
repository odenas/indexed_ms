
import glob
from setuptools import *

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
        install_requires=['numpy > 1.5']
)
