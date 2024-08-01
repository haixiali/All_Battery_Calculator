import math
import time

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.sparse import bsr_matrix
from Gmsh_Mesh_Reader import Gmsh_Mesh_Reader

import Node
import Element_D1O2


def f(x, y):
    return y * np.sin(x)


"""
Geometry -> Mesh -> Domain

Physics
"""

if __name__ == "__main__":
    gmsh_mesh = Gmsh_Mesh_Reader("./t1.msh")
    print(gmsh_mesh)
    domain = gmsh_mesh.read_to_region_dimension1()
    print(domain)
