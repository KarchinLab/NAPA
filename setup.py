from setuptools import setup

setup(name='napa',
      version='0.81',
      description='Network Analysis of Protein Adaptation (intra-protein residue coevolution network construction and analysis)',
      url='http://karchinlab.org/napa',
      author='Violeta Beleva-Guthrie',
      author_email='vbeleva@gmail.com',
      license='',
      packages=['napa'],
      install_requires=[
          'numpy',
          'scipy',
          'networkx',
          'ete2'
      ],
      zip_safe=False)