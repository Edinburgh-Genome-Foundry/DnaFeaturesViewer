import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, find_packages

exec(open('dna_features_viewer/version.py').read()) # loads __version__

setup(name='dna_features_viewer',
      version=__version__,
      author='Zulko',
    description='Plot features from DNA sequences (e.g. Genbank) with Python',
    long_description=open('README.rst').read(),
    license='see LICENSE.txt',
    keywords="DNA Sequence Feature Genbank Biopython Matplotlib",
    packages= find_packages(exclude='docs'),
    install_requires=["matplotlib"])
