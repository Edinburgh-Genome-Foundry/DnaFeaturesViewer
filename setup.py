import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, find_packages

exec(open('dna_features_viewer/version.py').read()) # loads __version__

setup(name='dna_features_viewer',
      version=__version__,
      author='Zulko',
    description='Plot features from DNA sequences (e.g. Genbank) with Python',
    long_description=open('pypi-readme.rst').read(),
    url='https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer',
    license='MIT',
    keywords="DNA Sequence Feature Genbank Biopython Matplotlib",
    packages= find_packages(exclude='docs'),
    install_requires=["matplotlib>=3", "Biopython"])
