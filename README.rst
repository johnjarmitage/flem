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


.. image:: https://github.com/johnjarmitage/io-page/blob/master/static/images/flem.gif


A simple diffusive landscape evolution model


* Free software: MIT license
* Documentation: https://flem.readthedocs.io.

Requirements
------------
This library requires an installation of fenics (https://fenicsproject.org/)

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

List of things to do:

1. Add the choice to change precipitation rates
2. Add the choice for boundary conditions
