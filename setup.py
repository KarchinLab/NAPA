from setuptools import setup;
from setuptools import find_packages;

def readme():
    with open('README.rst') as f:
        return f.read()


setup(name='napa',
      version='1.2',
      description='Network Analysis of Protein Adaptation (intra-protein residue coevolution network construction and analysis)',
      long_description=readme(),
      url='http://karchinlab.org/apps/appNapa.html',
      author='Violeta Beleva Guthrie',
      author_email='vbeleva@gmail.com',
      license='',
      packages=find_packages(exclude=['tests']),
      install_requires=[
          'joblib',
          'numpy',
          'scipy',
          'networkx',
          'ete2',
          'pyyaml'
      ],
      classifiers=[
      'Development Status :: 5 - Production/Stable',
      'Intended Audience :: Developers',
      'Topic :: Scientific/Engineering :: Bio-Informatics'
      ],
      zip_safe=True)
