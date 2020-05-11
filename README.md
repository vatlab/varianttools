[![PyPI version](https://badge.fury.io/py/variant-tools.svg)](https://badge.fury.io/py/variant-tools)

# Variant Tools

A command line tool for the manipulation, annotation, and analysis of genetic variants
from next-generation sequencing studies.

# Installation

If you are using a conda environment, you can install variant tools with command

```
conda install variant_tools -c bioconda -c conda-forge
```
Option `-c conda-forge` is required to enforce the use of `conda-forge` version of dependencies (e.g. `boost-cpp`) over their counterpoarts in the base channel.

Otherwise, you can try to install it through `pip`

```
pip install variant_tools
```

You will need to install

* `libboost`
* `gsl`
* `numpy`
* `Cython`
* `hdf5`
* `blosc`

* A C++ compiler such as `gcc`

which, in a conda environment, could be installed with command

```
conda install -c conda-forge boost-cpp gsl numpy cython blosc
```

Finally, you can checkout the latest version of variant tools from github, and install it with command

```
python setup.py install
```
with the aforementioned packages installed.

# Documentation

Please refer to [Variant Tools documentation](https://vatlab.github.io/vat-docs/) for details.
