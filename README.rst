====
flem
====


.. image:: https://img.shields.io/pypi/v/flem.svg
        :target: https://pypi.python.org/pypi/flem

.. image:: https://img.shields.io/travis/johnjarmitage/flem.svg
        :target: https://travis-ci.org/johnjarmitage/flem

.. image:: https://readthedocs.org/projects/flem/badge/?version=latest
        :target: https://flem.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://mybinder.org/badge_logo.svg
 :target: https://mybinder.org/v2/gh/johnjarmitage/flem/master?filepath=executable_article


.. image:: https://github.com/johnjarmitage/io-page/blob/master/static/images/flem.gif


A simple diffusive landscape evolution model


* Free software: MIT license
* Documentation: https://flem.readthedocs.io.

Requirements
------------
This library requires an installation of fenics
(https://fenicsproject.org/). It also requires gdal
(https://gdal.org/download.html). These can be installed through
conda using the requirements_conda.txt. There is a further
requirements.txt for libraries that are then installed through pip.


Features
--------

This is a set of functions written in python to solve for sediment
transport. At the base it solves the concentrative-diffusive equations
described by Smith & Bretherton (1972) [1]. These are solved using a
simple finite element scheme using the fenics library. Surface run-off
is routed either from model node-to-node or cell-to-cell (see
Armitage, 2019) [2].

The model can be started either with a initial condition of a uniform
elevation with some noise added, or a SRTM 30m DEM defined by west,
south, east, north coordinates.

[1] - https://doi.org/10.1029/WR008i006p01506

[2] - https://doi.org/10.5194/esurf-7-67-2019

List of things to do:

1. Add the choice to change precipitation rates
2. Add the choice for boundary conditions
