"""
flem
@author armitage
"""

from flem.read_dem import *
from fenics import Constant, Point, FunctionSpace, interpolate, SubDomain, DirichletBC, near
from mshr import Rectangle, generate_mesh
import numpy as np


def initialise(dem, bounds, res):
    """
    Function to initialise the model space

    :param dem: 0 = no input DEM; 1 = input DEM
    :param bounds: west, south, east, north - if no DEM [0, 0, lx, ly]; if DEM [west, south, east, north]
    :param res: model resolution along the y-axis.
    :return model_space, u_n, mesh, V, bc:
    """

    if dem == 0:
        lx = bounds[1]
        ly = bounds[3]

        # Create mesh and define function space
        domain = Rectangle(Point(0, 0), Point(lx / ly, ly / ly))
        mesh = generate_mesh(domain, res)

        V = FunctionSpace(mesh, 'P', 1)

        # Define initial value
        u_D = Constant(0)
        eps = 10 / ly
        u_n = interpolate(u_D, V)
        u_n.vector().set_local(u_n.vector().get_local() + eps * np.random.random(u_n.vector().size()))

    if dem == 1:
        u_n, lx, ly, mesh, V = read_dem(bounds, res)

    # boundary conditions

    class East(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[0], lx/ly)

    class West(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[0], 0.0)

    class North(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[1], ly/ly)

    class South(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[1], 0.0)

    # Should make this into an option!

    bc = [DirichletBC(V, u_n, West()),
          DirichletBC(V, u_n, East())]

    # def boundary(x, on_boundary):
    #   return on_boundary
    # bc = DirichletBC(V, u_n, boundary)

    model_space = [lx, ly, res]

    return model_space, u_n, mesh, V, bc
