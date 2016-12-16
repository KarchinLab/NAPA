from setuptools import setup;
from setuptools import find_packages;

setup(name='napa',
      version='0.87',
      description='Network Analysis of Protein Adaptation (intra-protein residue coevolution network construction and analysis)',
      url='http://karchinlab.org/napa',
      author='Violeta Beleva-Guthrie',
      author_email='vbeleva@gmail.com',
      license='',
      packages=find_packages(exclude=['tests']),
      install_requires=[
          'dill',
          'joblib',
          'numpy',
          'scipy',
          'networkx',
          'ete2'
      ],
      classifiers=[
      'Development Status :: 4 - Beta',
      'Intended Audience :: Developers',
      'Topic :: Scientific/Engineering :: Bio-Informatics'
      ],
      zip_safe=True)
