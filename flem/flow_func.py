"""

flem
@author: armitage

"""

from fenics import Function, vertex_to_dof_map, vertices, Edge, dof_to_vertex_map, cells, facets
import numpy as np


def sd_cellcell(mesh, V, u_n, De, nexp):
    """
    SD cell-to-cell
    Flow routing from cell-to-cell based on the steepest route of descent

    :param mesh: mesh object generated using mshr (fenics)
    :param V: finite element function space
    :param u_n: solution (trial function) for water flux
    :param De: dimensionless diffusion coefficient
    :param nexp: water flux exponent
    :return:
    """

    # get a map of neighbours (thanks google!)
    tdim = mesh.topology().dim()
    mesh.init(tdim - 1, tdim)
    cell_neighbors = np.array([sum((list(filter(lambda ci: ci != cell.index(),
                                           facet.entities(tdim)))
                                    for facet in facets(cell)), [])
                               for cell in cells(mesh)])
    
    # first get the elevation area of each element
    dofmap = V.dofmap()
    elevation = []
    area = []
    xm = []
    ym = []
    for cell in cells(mesh):
        cellnodes = dofmap.cell_dofs(cell.index())
        elevation.append(sum(u_n.vector()[cellnodes])/3)
        area.append(cell.volume())
        p = cell.midpoint()
        xm.append(p.x())
        ym.append(p.y())
    elevation = np.array(elevation)
    area = np.array(area)
    xm = np.array(xm)
    ym = np.array(ym)

    # now sort the vector of elevations by decending topography
    ind = np.argsort(-elevation)
    sorted_neighbors = cell_neighbors[ind]

    # determine length between elements
    steep_len = []
    for cell in cells(mesh):
        xh = xm[cell.index()]
        yh = ym[cell.index()]
        
        neicells = cell_neighbors[cell.index()]
        tnei = elevation[neicells]
        imin = np.argmin(tnei)
        ncell = neicells[imin]
        xn = xm[ncell]
        yn = ym[ncell]
  
        steep_len.append(np.sqrt((xh-xn)*(xh-xn)+(yh-yn)*(yh-yn)))
    steep_len = np.array(steep_len)

    flux = area/steep_len

    # determine flux from highest to lowest cells
    for cell in cells(mesh):
        neicells = sorted_neighbors[cell.index()]
        tnei = elevation[neicells]
        imin = np.argmin(tnei)
        ncell = neicells[imin]

        if elevation[ncell] < elevation[ind[cell.index()]]:
            flux[ncell] = flux[ncell] + flux[ind[cell.index()]]

    # interpolate to the nodes
    gc = mesh.coordinates()

    flux_node = np.zeros(len(gc))
    for cell in cells(mesh):
        cellnodes = dofmap.cell_dofs(cell.index())
        for nod in cellnodes:
            flux_node[nod] = flux_node[nod] + flux[cell.index()]/3

    q = Function(V)
    q.vector()[:] = 1 + De*pow(flux_node, nexp)

    return q


def mfd_cellcell(mesh, V, u_n, De, nexp):
    """
    MFD cell-to-cell
    Flow routing distributed from cell-to-cell

    :param mesh: mesh object generated using mshr (fenics)
    :param V: finite element function space
    :param u_n: solution (trial function) for water flux
    :param De: dimensionless diffusion coefficient
    :param nexp: water flux exponent
    :return:
    """

    # get a map of neighbours (thanks google!)
    tdim = mesh.topology().dim()
    mesh.init(tdim - 1, tdim)
    cell_neighbors = np.array([sum((list(filter(lambda ci: ci != cell.index(),
                                           facet.entities(tdim)))
                                    for facet in facets(cell)), [])
                               for cell in cells(mesh)])
    
    # first get the elevation area of each element
    dofmap = V.dofmap()
    elevation = []
    area = []
    xm = []
    ym = []
    for cell in cells(mesh):
        cellnodes = dofmap.cell_dofs(cell.index())
        elevation.append(sum(u_n.vector()[cellnodes])/3)
        area.append(cell.volume())
        p = cell.midpoint()
        xm.append(p.x())
        ym.append(p.y())
    elevation = np.array(elevation)
    area = np.array(area)
    xm = np.array(xm)
    ym = np.array(ym)

    # now sort the vector of elevations by decending topography
    ind = np.argsort(-elevation)
    sorted_neighbors = cell_neighbors[ind]

    # determine length between elements
    steep_len = []
    for cell in cells(mesh):
        xh = xm[cell.index()]
        yh = ym[cell.index()]
        
        neicells = cell_neighbors[cell.index()]
        tnei = elevation[neicells]
        imin = np.argmin(tnei)
        ncell = neicells[imin]
        xn = xm[ncell]
        yn = ym[ncell]
  
        steep_len.append(np.sqrt((xh-xn)*(xh-xn)+(yh-yn)*(yh-yn)))
    steep_len = np.array(steep_len)

    flux = area/steep_len

    # determine flux from highest to lowest cells
    for cell in cells(mesh):
        neicells = sorted_neighbors[cell.index()]
        tnei = elevation[neicells]
        imin = np.argmin(tnei)
        ncell = neicells[imin]
                
        weight = np.zeros(len(neicells))
        i = 0
        for neicell in neicells:
            weight[i] = elevation[ind[cell.index()]] - elevation[neicell]
            # downhill only
            if weight[i] < 0:
              weight[i] = 0
            i += 1
            
        # weight flux by the sum of the lengths down slope
        if max(weight) > 0:
            weight = weight/sum(weight)
        else:
            weight[:] = 0
        i = 0
        for neicell in neicells:
            flux[neicell] = flux[neicell] + flux[ind[cell.index()]]*weight[i]
            i += 1

    # interpolate to the nodes
    gc = mesh.coordinates()

    flux_node = np.zeros(len(gc))
    for cell in cells(mesh):
        cellnodes = dofmap.cell_dofs(cell.index())
        for nod in cellnodes:
            flux_node[nod] = flux_node[nod] + flux[cell.index()]/3

    q = Function(V)
    q.vector()[:] = 1 + De*pow(flux_node, nexp)

    return q


def mfd_nodenode(mesh, V, u_n, De, nexp):
    """
    MFD node-to-node
    Flow routing distributed from node-to-node

    :param mesh: mesh object generated using mshr (fenics)
    :param V: finite element function space
    :param u_n: solution (trial function) for water flux
    :param De: dimensionless diffusion coefficient
    :param nexp: water flux exponent
    :return:
    """

    # get the global coordinates
    gdim = mesh.geometry().dim()
#    if dolfin.dolfin_version() == '1.6.0':
#      dofmap = V.dofmap()
#      gc = dofmap.tabulate_all_coordinates(mesh).reshape((-1,gdim))
#    else:
    gc = V.tabulate_dof_coordinates().reshape((-1,gdim))
    vtd = vertex_to_dof_map(V)
     
    # first get the elevation of each vertex
    elevation = np.zeros(len(gc))
    elevation = u_n.compute_vertex_values(mesh)
          
    # loop to get the local flux
    mesh.init(0,1)
    flux = np.zeros(len(gc))
    neighbors = []
    for v in vertices(mesh):
        idx = v.index()
        
        # get the local neighbourhood
        neighborhood = [Edge(mesh, i).entities(0) for i in v.entities(1)]
        neighborhood = np.array(neighborhood).flatten()
          
        # Remove own index from neighborhood
        neighborhood = neighborhood[np.where(neighborhood != idx)[0]]
        neighbors.append(neighborhood)
                        
        # get location
        xh = v.x(0)
        yh = v.x(1)
        
        # get distance to neighboring vertices
        length = np.zeros(len(neighborhood))
        weight = np.zeros(len(neighborhood))
        i = 0
        for vert in neighborhood:
            nidx = vtd[vert]
            xn = gc[nidx,0]
            yn = gc[nidx,1]
            length[i] = np.sqrt((xh-xn)*(xh-xn)+(yh-yn)*(yh-yn))
            flux[vert] = length[i]
#            weight[i] = elevation[idx] - elevation[vert]
#            # downhill only
#            if weight[i] < 0:
#              weight[i] = 0
#            i += 1
#
#        # weight flux by the sum of the lengths down slope
#        if max(weight) > 0:
#            weight = weight/sum(weight)
#        else:
#            weight[:] = 0
#        i = 0
#        for vert in neighborhood:
#            flux[vert] = flux[vert] + length[i]*weight[i]
#            i += 1
            
    # sort from top to botton
    sortedidx = np.argsort(-elevation)
    
    # accumulate fluxes from top to bottom
    for idx in sortedidx:
        neighborhood = neighbors[idx]
        weight = np.zeros(len(neighborhood))
        i = 0
        for vert in neighborhood:
            weight[i] = elevation[idx] - elevation[vert]
            # downhill only
            if weight[i] < 0:
              weight[i] = 0
            i += 1
            
        # weight flux by the sum of the lengths down slope
        if max(weight) > 0:
            weight = weight/sum(weight)
        else:
            weight[:] = 0
        i = 0
        for vert in neighborhood:
            flux[vert] = flux[vert] + flux[idx]*weight[i]
            i += 1

    # calculate the diffusion coefficient
    q0 = 1 + De*pow(flux,nexp)
    q = Function(V)
    q.vector()[:] = q0[dof_to_vertex_map(V)]

    return q


def sd_nodenode(mesh, V, u_n, De, nexp):
    """
    SD node-to-node
    Flow routing from node-to-node based on the steepest route of descent

    :param mesh: mesh object generated using mshr (fenics)
    :param V: finite element function space
    :param u_n: solution (trial function) for water flux
    :param De: dimensionless diffusion coefficient
    :param nexp: water flux exponent
    :return:
    """

    # get the global coordinates
    gdim = mesh.geometry().dim()
    #if dolfin.dolfin_version() == '1.6.0':
    #  dofmap = V.dofmap()
    #  gc = dofmap.tabulate_all_coordinates(mesh).reshape((-1,gdim))
    #else:
    gc = V.tabulate_dof_coordinates().reshape((-1,gdim))
    vtd = vertex_to_dof_map(V)
     
    # first get the elevation of each vertex
    elevation = np.zeros(len(gc))
    elevation = u_n.compute_vertex_values(mesh)
          
    # loop to get the local flux
    mesh.init(0,1)
    flux = np.zeros(len(gc))
    neighbors = []
    for v in vertices(mesh):
        idx = v.index()
        
        # get the local neighbourhood
        neighborhood = [Edge(mesh, i).entities(0) for i in v.entities(1)]
        neighborhood = np.array(neighborhood).flatten()
          
        # Remove own index from neighborhood
        neighborhood = neighborhood[np.where(neighborhood != idx)[0]]
        neighbors.append(neighborhood)
                        
        # get location
        xh = v.x(0)
        yh = v.x(1)
        
        # get distance to neighboring vertices
        length = np.zeros(len(neighborhood))
        weight = np.zeros(len(neighborhood))
        i = 0
        for vert in neighborhood:
            nidx = vtd[vert]
            xn = gc[nidx,0]
            yn = gc[nidx,1]
            length[i] = np.sqrt((xh-xn)*(xh-xn)+(yh-yn)*(yh-yn))
            flux[vert] = length[i]
#            weight[i] = elevation[idx] - elevation[vert]
#            # downhill only
#            if weight[i] < 0:
#              weight[i] = 0
#            i += 1
#
#        # find steepest slope
#        steepest = len(neighborhood)+2
#        if max(weight) > 0:
#            steepest = np.argmax(weight)
#        else:
#            weight[:] = 0
#        i = 0
#        for vert in neighborhood:
#            if i == steepest:
#                weight[i] = 1
#            else:
#                weight[i] = 0
#            flux[vert] = flux[vert] + length[i]*weight[i]
#            i += 1
            
    # sort from top to botton
    sortedidx = np.argsort(-elevation)
    
    # accumulate fluxes from top to bottom
    for idx in sortedidx:
        neighborhood = neighbors[idx]
        weight = np.zeros(len(neighborhood))
        i = 0
        for vert in neighborhood:
            weight[i] = elevation[idx] - elevation[vert]
            # downhill only
            if weight[i] < 0:
              weight[i] = 0
            i += 1

        # find steepest slope
        steepest = len(neighborhood)+2
        if max(weight) > 0:
            steepest = np.argmax(weight)
        else:
            weight[:] = 0
        i = 0
        for vert in neighborhood:
            if i == steepest:
                weight[i] = 1
            else:
                weight[i] = 0
            flux[vert] = flux[vert] + flux[idx]*weight[i]
            i += 1
            
    # calculate the diffusion coefficient
    q0 = 1 + De*pow(flux,nexp)
    q = Function(V)
    q.vector()[:] = q0[dof_to_vertex_map(V)]

    return q
