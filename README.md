# eqstats

A python package for performing statistical tests on earthquake catalogs

## Overview

This package contains code useful for performing statistical tests on earthquake catalogs. This includes
performing the statistical tests themselves, as well as tools for generating random catalogs in order
to calibrate tests. The tests and methods are outlined in this (paper)[https://doi.org/10.1002/2014JB011777].
If you use this code in any further work, please cite this paper as well as this package.

## Installation

To install the package, from the directory where you have downloaded the source code enter the following
into the shell:

    $ python setup.py install
    
This should copy the necessary files into your Python installation. The code requires `numpy`, `scipy`,
and `statsmodels`, all of which can be installed using `pip`. I have mainly used the code with
Python 3.

## Usage

The package contains two submodules: `catalogs` and `stattests`.

* `catalogs` has functions for generating random event times, event magnitudes that follow a
  Gutenberg-Richter magnitude/frequency distribution, and event times that follow an Omori-like
  aftershock time decay. These can be combined to create synthetic catalogs that are random
  and/or contain aftershocks and can approximate real seismicity for calibrating statistical
  tests.
* `stattests` contains the various functions for computing test statistics on a given catalog.
  Each test requires inputs specific to that test, so it is up to the user to compute the correct
  inputs. Most tests require either recurrence or occurrence times of events, and some require
  hyperparameters such as number of bins. One test (`bigtrig_test`) requires event magnitudes
  as well as a magnitude threshold, as it makes the assumption that large events lead to a
  detectable increase in seismicity. This submodule also has functions for computing p-values
  and detection power.
