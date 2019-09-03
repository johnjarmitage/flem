"""
FEniCS program: Smith & Bretherton (1972) equations

  u'= nabla((1+De*q_w).nabla(u)) + f  in the domain
  u = 0                               on the boundaries
  u = 0                               at t = 0

  u = elevation
  f = uplift

@author: armitage
"""

from __future__ import print_function
from fenics import Constant,Point,FunctionSpace,interpolate,SubDomain,DirichletBC,near,TrialFunction,TestFunction,\
    Expression,dx,dot,grad,lhs,rhs,solve,assemble,File
from mshr import Rectangle,generate_mesh
import numpy as np
from flem.flow_func import *
import matplotlib.pyplot as plt
import peakutils # https://zenodo.org/badge/latestdoi/102883046


def solve_flem(model_space, physical_space, flow, u_n, mesh, V, bc, dt, num_steps, out_time, plot, statistics, name):
    """
    Solve for landscape evolution

    This function does hte hard work. First the model domain is created. Then we loop through time and solve the
    diffusion equation to solve for landscape evolution. Output can be saved as vtk files at every "out_time" specified.
    Plots using fenics inbuilt library can be visualised at every "plot_time"

    This function returns a 1d numpy array of time, sediment flux and if statistics is turned on a 2d numpy array of
    the final wavelength of the landscape.

    :param model_space: list of domain variables, [lx,ly,res]
    :param physical_space: list of physical parameters, [kappa, c, nexp, alpha, U]
    :param flow: 0 = MFD node-to-node; 1 = MFD cell-to-cell; 2 = SD node-to-node; 3 = SD cell-to-cell
    :param u_n: elevation function
    :param mesh: dolphyn mesh
    :param V: fenics functionspace
    :param bc: boundary conditions
    :param dt: time step size in years
    :param num_steps: number of time steps
    :param out_time: time steps to output vtk files (0=none)
    :param plot: plot sediment flux (0=off,1=on)
    :param statistics: output statistics of landscape (0=off,1=on)
    :param name: directory name for output vtk files
    :return: sed_flux, time, wavelength
    """

    # Domain dimensions
    lx = model_space[0]
    ly = model_space[1]

    # Physical parameters
    kappa = physical_space[0]           # diffusion coefficient
    c = physical_space[1]               # discharge transport coefficient
    nexp = physical_space[2]            # discharge exponent
    alpha = physical_space[3]           # precipitation rate
    De = c*pow(alpha*ly, nexp)/kappa
    uamp = physical_space[4]*ly/kappa   # uplift

    dt = dt*kappa/(ly*ly)               # time step size

    sed_flux = np.zeros(num_steps)      # array to store sediment flux
    time = np.zeros(num_steps)

    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)
    f = Constant(uamp)

    # 0 = MFD node-to-node; 1 = MFD cell-to-cell; 2 = SD node-to-node; 3 = SD cell-to-cell
    if flow == 0:
        q_n = mfd_nodenode(mesh, V, u_n, De, nexp)
    if flow == 1:
        q_n = mfd_cellcell(mesh, V, u_n, De, nexp)
    if flow == 2:
        q_n = sd_nodenode(mesh, V, u_n, De, nexp)
    if flow == 3:
        q_n = sd_cellcell(mesh, V, u_n, De, nexp)

    F = u*v*dx + dt*q_n*dot(grad(u), grad(v))*dx - (u_n + dt*f)*v*dx
    a, L = lhs(F), rhs(F)

    # Solution and sediment flux
    u = Function(V)
    q_s = Expression('u0 + displ - u1', u0=u_n, displ=Constant(uamp*dt), u1=u, degree=2)

    # Iterate
    t = 0
    i = 0
    for n in range(num_steps):

        # This needs to become an option!
        # Double rain fall
        # if n == 501:
        #   alpha = 2
        #   De    = c*pow(alpha*ly,nexp)/kappa

        # Update current time
        t += dt

        # Compute solution
        solve(a == L, u, bc)

        # Calculate sediment flux
        sed_flux[i] = assemble(q_s*dx(mesh))
        time[i] = t
        i += 1

        # Update previous solution
        u_n.assign(u)

        # Update flux
        # 0 = MFD node-to-node; 1 = MFD cell-to-cell; 2 = SD node-to-node; 3 = SD cell-to-cell
        if flow == 0:
            q = mfd_nodenode(mesh, V, u_n, De, nexp)
        if flow == 1:
            q = mfd_cellcell(mesh, V, u_n, De, nexp)
        if flow == 2:
            q = sd_nodenode(mesh, V, u_n, De, nexp)
        if flow == 3:
            q = sd_cellcell(mesh, V, u_n, De, nexp)
        q_n.assign(q)

        # Output solutions
        if out_time != 0:
            if np.mod(n, out_time) == 0:
                filename = '%s/u_solution_%d.pvd' % (name, n)
                vtkfile = File(filename)
                vtkfile << u
                filename = '%s/q_solution_%d.pvd' % (name, n)
                vtkfile = File(filename)
                vtkfile << q

    # Post processing
    if plot != 0:
        plt.plot(time*1e-6*ly*ly/kappa, sed_flux/dt*kappa, 'k', linewidth=2)
        plt.xlabel('Time (Myr)')
        plt.ylabel('Sediment Flux (m^2/yr)')
        sedname = '%s/sed_flux_%d.svg' % (name, model_space[2])
        plt.savefig(sedname, format='svg')
        plt.clf()

    if out_time != 0:
        # Output last elevation
        filename = '%s/u_solution_%d_%d.pvd' % (name, model_space[2], n)
        vtkfile = File(filename)
        u.rename("elv", "elevation")
        vtkfile << u

        # Output last water flux
        filename = '%s/q_solution_%d_%d.pvd' % (name, model_space[2], n)
        vtkfile = File(filename)
        q.rename("flx", "flux")
        vtkfile << q

    # Calculate valley spacing from peak to peak in water flux
    tol = 0.001  # avoid hitting points outside the domain
    y = np.linspace(0 + tol, 1 - tol, 100)
    x = np.linspace(0.01, lx/ly-0.01, 20)
    # x = np.linspace(0.01, 0.21, 20)
    wavelength = np.zeros(len(x))
    if statistics != 0:
        i = 0
        for ix in x:
            points = [(ix, y_) for y_ in y]  # 2D points
            q_line = np.array([q(point) for point in points])

            indexes = peakutils.indexes(q_line, thres=0.05, min_dist=5)
            if len(indexes) > 1:
                wavelength[i] = sum(np.diff(y[indexes]))/(len(indexes)-1)
            else:
                wavelength[i] = 0
            i += 1

        if plot != 0:
            plt.plot(y*1e-3*ly, q_line*kappa/ly, 'k', linewidth=2)
            plt.plot(y[indexes]*1e-3*ly, q_line[indexes]*kappa/ly, '+r')
            plt.xlabel('Distance (km)')
            plt.ylabel('Water Flux (m/yr)')
            watername = '%s/water_flux_spacing_%d.svg' % (name, model_space[2])
            plt.savefig(watername, format='svg')
            plt.clf()

    return sed_flux, time, wavelength

