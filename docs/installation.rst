.. highlight:: shell

============
Installation
============

From sources
------------

The sources for flem can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

  $ git clone git://github.com/johnjarmitage/flem

Or download the `tarball`_:

.. code-block:: console

  $ curl  -OL https://github.com/johnjarmitage/flem/tarball/master

Dependencies
------------

Once you have a copy of the source you will then need to install the
dependencies. `Fenics <https://fenicsproject.org/>`_ is best installed
from Conda, as is `GDAL <https://gdal.org/>`_. Therefore provided are two
requirements files. First install `Anaconda <https://www.anaconda.com/>`_.
Then from within a terminal:

.. code-block:: console

  $ conda create -n flem python=3.6
  $ conda activate flem

Then install the dependencies via conda:

.. code-block:: console

  $ conda install -c conda-forge --file requirements_conda.txt

Finally install the the remaining dependencies via PyPI:

.. code-block:: console

  $ conda install pip
  $ pip install -r requirements.txt

Now you are good to go. Check out the Notebook in the tests directory for
ideas.

.. _Github repo: https://github.com/johnjarmitage/flem
.. _tarball: https://github.com/johnjarmitage/flem/tarball/master
