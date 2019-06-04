# -*- coding: utf-8 -*-

"""Top-level package for flem."""

__author__ = """John Armitage"""
__email__ = 'john.joseph.armitage@gmail.com'
__version__ = '0.1.0'


from .flow_func import sd_cellcell,mfd_cellcell,sd_nodenode,mfd_nodenode
from flem.flow_func import sd_cellcell,mfd_cellcell,sd_nodenode,mfd_nodenode
from .read_dem import read_dem
from flem.read_dem import read_dem
from .initialise import initialise
from flem.initialise import initialise
from .solve_flem import solve_flem

__all__ = ["sd_cellcell","mfd_cellcell","sd_nodenode","mfd_nodenode","solve_flem","initialise","read_dem"]
