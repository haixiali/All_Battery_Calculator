import math
import time

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.sparse import bsr_matrix

if __name__ == "__main__":
    sps_mat = []
    num_mat = 2
    num_node = 1000_000
    for i in range(num_mat):
        sps_mat.append(sp.sparse.rand(num_node, num_node, density=0.0001, format="bsr"))
        print(f"{i / num_mat * 100:.2f}%", end=", ")

    print(f"Done!")

    sps_mat = sps_mat * 10_000

    time_start = time.time()
    i = 0
    j = 0
    s = sp.sparse.bsr_matrix((num_node, num_node))
    for mat in sps_mat:
        s = s + mat
        if j > len(sps_mat) / 20:
            j = 0
            print(f"{i / len(sps_mat) * 100:.2f}%", end=", ")

        j = j + 1
        i = i + 1

    print(f"time for summation: {time.time() - time_start}")
