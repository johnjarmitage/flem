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

Installation
------------

- flem requires Python 3.7
- flem requires fenics, gdal, and a bit more. Fenics is best installed using conda.
  Therefore before installing first get yourself
  `Anaconda <https://www.anaconda.com/distribution/#download-section>`_ (the 3.7 version) or
  if you prefer it light, `miniconda <https://www.anaconda.com/distribution/#download-section>`_.
- create a directory of your choice and create an ``environment.yml`` file containing the
  following:

::

  name: flem
  channels:
    - conda-forge
    - defaults
  dependencies:
    # flem requires
    # need to be specific for mshr and fenics
    - fenics=2019.1.0=py37_1
    - mshr=2019.1.0=py37h7596e34_1000
    - gdal
    - peakutils
    - matplotlib
    - scipy
    - pip
    - pip:
      # flem requires
      - flem
      - elevation
  prefix: /srv/conda

- from the terminal run: ``conda env create -f environment.yml``
- check out this `notebook <https://github.com/johnjarmitage/flem-examples>`_
  for how to run flem.
- or see `run_models.py <https://github.com/johnjarmitage/flem/blob/master/run_models.py>`_
  for a more clunky example.

What is flem?
-------------

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
