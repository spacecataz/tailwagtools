# TailWagTools

This is a basic package for performing investigations into the possibility that
asymmetric ionospheric outflow can change the deflection of the magnetotail
(or, rather, "wag" the tail.)  It currently provides tools to examine
observed magnetic field values from Cluster against empirical fields from
the Tsyganenko model.  More functionality is being added.

This package is the result of undergraduate research at the University of Texas
at Arlington.

## Installation
Ensure that **tailwag.py** can be found via `PYTHONPATH`.

## Prerequisites
Prereqs are the Spacepy package and its prereqs:
Python, Numpy, Matplotlib, Scipy, and the NASA CDF library.

|Package | Version |
|--------|---------|
|Spacepy | 0.2.0|
|Python | 3.X  |
|Numpy  | 1.16.X |
|Matplotlib | 3.1.X |
|Scipy | 1.3.X |
|NASA CDF | >=3.6.X |
|Pandas | 1.0.1 |


Older versions may work, but I'm not about to test that today.

## Usage
The main module is **tailwag.py**.  This module provides functions and classes
to perform the requisite analysis.  Scripts to complete the analysis can be
built with this module.

If **tailwag.py** is executed as a script, it will perform a short visual
test of functionality.

## Included Data
The **sample_data** directory contains enough data to test functionality.

|File  | Contents|
|------|--------|
|example_cluster.cdf | Cluster orbit and magnetic field data |
