import os
import sys

from setuptools import setup

# Ordinary dependencies
DEPENDENCIES = []
with open("requirements/requirements-all.txt", "r") as reqs_file:
    for line in reqs_file:
        if not line.strip():
            continue
        # DEPENDENCIES.append(line.split("=")[0].rstrip("<>"))
        DEPENDENCIES.append(line)

# Additional keyword arguments for setup()
extra = {"install_requires": DEPENDENCIES}

with open("refgenie/_version.py", "r") as versionfile:
    version = versionfile.readline().split()[-1].strip("\"'\n")

with open("README.md") as f:
    long_description = f.read()

setup(
    name="refgenie",
    packages=["refgenie"],
    version=version,
    description="Refgenie creates a standardized folder structure for reference genome files and indexes. "
    "You can download pre-built genomes or build your own for any fasta file",
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license="BSD2",
    entry_points={
        "console_scripts": [
            "refgenie = refgenie.__main__:main",
            "import_igenome = refgenie.add_assets_igenome:main",
        ],
    },
    keywords="bioinformatics, sequencing, ngs",
    package_data={"refgenie": [os.path.join("refgenie", "*")]},
    include_package_data=True,
    url="http://refgenie.databio.org",
    author=u"Nathan Sheffield, Vince Reuter, Michal Stolarczyk",
    **extra
)
