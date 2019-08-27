"""

flem
@author: armitage

"""

import numpy as np
import elevation as elv
import os
from osgeo import gdal
from fenics import FunctionSpace, Function, Point
from mshr import Rectangle, generate_mesh
from scipy import interpolate


def read_dem(bounds, res):
    """
    Function to read in a DEM from SRTM amd interplolate it onto a dolphyn mesh. This function uses the python package
    'elevation' (http://elevation.bopen.eu/en/stable/) and the gdal libraries.

    I will assume you want the 30m resolution SRTM model.

    :param bounds: west, south, east, north coordinates
    :return u_n, lx, ly: the elevation interpolated onto the dolphyn mesh and the lengths of the domain
    """
    west, south, east, north = bounds

    # Create a temporary file to store the DEM and go get it using elevation
    dem_path = 'tmp.tif'
    output = os.getcwd() + '/' + dem_path
    elv.clip(bounds=bounds, output=output, product='SRTM1')

    # read in the DEM into a numpy array
    gdal_data = gdal.Open(output)
    data_array = gdal_data.ReadAsArray().astype(np.float)

    # The DEM is 30m per pixel, so lets make a array for x and y at 30 m
    ny, nx = np.shape(data_array)
    lx = nx*30
    ly = ny*30

    x, y = np.meshgrid(np.linspace(0, lx/ly, nx), np.linspace(0, 1, ny))

    # Create mesh and define function space
    domain = Rectangle(Point(0, 0), Point(lx/ly, 1))
    mesh = generate_mesh(domain, res)
    V = FunctionSpace(mesh, 'P', 1)
    u_n = Function(V)

    # Get the global coordinates
    gdim = mesh.geometry().dim()
    gc = V.tabulate_dof_coordinates().reshape((-1, gdim))

    # Interpolate elevation into the initial condition
    elevation = interpolate.griddata((x.flatten(), y.flatten()), data_array.flatten(), (gc[:, 0], gc[:, 1]),
                                     method='nearest')
    u_n.vector()[:] = elevation

    # remove tmp DEM
    os.remove(output)

    return u_n, lx, ly, mesh, V
