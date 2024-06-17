import math
import time

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.sparse import bsr_matrix
from GmshMesh import GmshMesh


def f(x, y):
    return y * np.sin(x)


if __name__ == "__main__":
    gmsh_mesh = GmshMesh("./t1.msh")
    print(gmsh_mesh)
