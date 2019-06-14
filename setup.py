import os
from setuptools import setup
import sys

# Ordinary dependencies
DEPENDENCIES = []
with open("requirements/requirements-all.txt", "r") as reqs_file:
    for line in reqs_file:
        if not line.strip():
            continue
        #DEPENDENCIES.append(line.split("=")[0].rstrip("<>"))
        DEPENDENCIES.append(line)

# Additional keyword arguments for setup()
extra = {"install_requires": DEPENDENCIES}

# 2to3
if sys.version_info >= (3, ):
    extra["use_2to3"] = True

with open("refgenie/_version.py", 'r') as versionfile:
    version = versionfile.readline().split()[-1].strip("\"'\n")

# Handle the pypi README formatting.
try:
    import pypandoc
    long_description = pypandoc.convert_file('README.md', 'rst')
except(IOError, ImportError, OSError):
    print("Warning -- couldn't convert README to rst")
    long_description = open('README.md').read()

setup(
    name='refgenie',
    packages=["refgenie"],
    version=version,
    description='Refgenie creates a standardized folder structure for reference genome files and indexes. '
                'You can download pre-built genomes or build your own for any fasta file',
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],    
    license="BSD2",
    entry_points={
        "console_scripts": [
            'refgenie = refgenie.refgenie:main'
        ],
    },
    keywords="bioinformatics, sequencing, ngs",
    package_data={"refgenie": [os.path.join("refgenie", "*")]},
    include_package_data=True,
    url='http://refgenie.databio.org',
    author=u'Nathan Sheffield, Vince Reuter, Michal Stolarczyk',
    **extra
)