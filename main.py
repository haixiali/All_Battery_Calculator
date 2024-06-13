import math
import time

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

np.set_printoptions(linewidth=800, formatter={"float": lambda x: f"{x:10.4g}"})

quadrature_weights = np.array(
    [
        0.2729250867779006,
        0.2628045445102467,
        0.2628045445102467,
        0.2331937645919905,
        0.2331937645919905,
        0.1862902109277343,
        0.1862902109277343,
        0.1255803694649046,
        0.1255803694649046,
        0.0556685671161737,
        0.0556685671161737,
    ]
)

quadrature_xi = np.array(
    [
        0.00000000000000000,
        -0.2695431559523450,
        0.26954315595234500,
        -0.5190961292068118,
        0.51909612920681180,
        -0.7301520055740494,
        0.73015200557404940,
        -0.8870625997680953,
        0.88706259976809530,
        -0.9782286581460570,
        0.97822865814605700,
    ]
)

"""
quadrature_weights = np.array([0.347854845, 0.652145155, 0.652145155, 0.347854845])

quadrature_xi = np.array([-0.861136312, -0.339981044, 0.339981044, 0.861136312])
"""

quadrature_xi_xi = quadrature_xi * quadrature_xi


# domain functions
def b_fun(x):
    return 0


def c_fun(x):
    return 1


def f_fun(x):
    return x


b_vfun = np.vectorize(b_fun)
c_vfun = np.vectorize(c_fun)
f_vfun = np.vectorize(f_fun)


def assemble_local_matrix(x_ele, robin_bc=None):
    x0 = x_ele[0]
    x2 = x_ele[1]
    x1 = 0.5 * (x0 + x2)

    a_xi = 0.5 * (x0 + x2) - x1
    b_xi = 0.5 * (x2 - x0)
    c_xi = x1

    x_quad = a_xi * quadrature_xi_xi + b_xi * quadrature_xi + c_xi
    jac = 2 * a_xi * quadrature_xi + b_xi

    aa_mat = np.diag(quadrature_weights * jac)
    bb_N = b_vfun(x_quad)
    bb_mat = np.diag(quadrature_weights * bb_N * jac)
    cc_N = c_vfun(x_quad)
    cc_mat = np.diag(quadrature_weights * cc_N * jac)
    f_vec = f_vfun(x_quad)
    f_mat = np.transpose([quadrature_weights * f_vec * jac])

    z01 = x0 - x1
    z12 = x1 - x2
    z20 = x2 - x0
    z_neg_rep = -1 / (z01 * z12 * z20)
    z_vec = z_neg_rep * np.array([z12, z20, z01])

    a_N = np.transpose([z_vec])
    b_N = -a_N * np.array([[x1 + x2, x0 + x2, x0 + x1]]).T
    c_N = a_N * np.array([[x1 * x2, x0 * x2, x0 * x1]]).T

    N = np.outer(a_N, x_quad * x_quad) + np.outer(b_N, x_quad) + np.outer(c_N, np.ones(x_quad.size))
    Np = np.outer(2 * a_N, x_quad) + np.outer(b_N, np.ones(x_quad.size))

    A = np.matmul(np.matmul(Np, aa_mat), Np.T)
    B = np.matmul(np.matmul(N, bb_mat), Np.T)
    C = np.matmul(np.matmul(N, cc_mat), N.T)
    K = A + B + C
    F = np.matmul(N, f_mat)

    """ print(f"-------- element: {i_ele}")
    print(f"mat A: \n{A}")
    print(f"mat B: \n{B}")
    print(f"mat C: \n{C}")
    print(f"mat K: \n{K}")
    print(f"mat F: \n{F}") """

    if robin_bc is not None:
        xi_bc = -1 if robin_bc[0] == 0 else 1
        x_bc = a_xi * xi_bc * xi_bc + b_xi * xi_bc + c_xi
        N_bc = np.outer(a_N, x_bc * x_bc) + np.outer(b_N, x_bc) + np.outer(c_N, np.ones(x_bc.size))

        rho = xi_bc * robin_bc[2]  # xi_bc handles the sign of R term
        R = rho * np.matmul(N_bc, N_bc.T)
        K = K + R

        S = robin_bc[3] * N_bc
        F = F + S

    return K, F


if __name__ == "__main__":

    s0 = 0
    s1 = 1
    # B.C. on left-most & right-most nodes, [loc, type, rho, f], loc=0 for left, rho=0 for Dirichlet B.C.
    boundary_conditions = [[0, 0, 0, 1],
                           [1, 2, 0, 1]]

    num_ele = 10_000

    # create elements
    num_node = 2 * num_ele + 1
    h_ele = (s1 - s0) / num_ele

    x_ele = np.array([np.arange(s0, s1, h_ele), np.arange(h_ele, s1 + h_ele, h_ele)]).T
    x_node = np.arange(s0, s1 + h_ele * 0.1, h_ele * 0.5)
    idx_node = np.array([np.arange(0, num_node - 1, 2), np.arange(1, num_node, 2), np.arange(2, num_node + 1, 2)]).T

    print(f"-------- assembling local matrix, element number: {num_ele}")
    time_start = time.time()

    K_loc = np.zeros((num_ele, 3, 3))
    F_loc = np.zeros((num_ele, 3, 1))

    for i_ele in range(num_ele):

        robin_bc = None
        if i_ele == 0 and boundary_conditions[0][1] == 2:
            robin_bc = boundary_conditions[0]
        if i_ele == num_ele - 1 and boundary_conditions[1][1] == 2:
            robin_bc = boundary_conditions[1]

        K_ele, F_ele = assemble_local_matrix(x_ele[i_ele], robin_bc)

        # print(f"---- element: {i_ele: 6d}")
        # print(f"mat K: \n{K_ele}")
        # print(f"mat F: \n{F_ele}")

        K_loc[i_ele] = K_ele
        F_loc[i_ele] = F_ele

    time_end = time.time()
    print(f"time elapsed: {time_end - time_start:.6f}")

    print(f"-------- assembling global matrix")
    time_start = time.time()

    K = sp.sparse.csr_matrix((num_node, num_node))
    F = sp.sparse.csr_matrix((num_node, 1))

    # K_glb = np.zeros((num_node, num_node, num_node))
    # F_glb = np.zeros((num_node, num_node, 1))
    time_sum1 = 0
    time_sum2 = 0
    time_sum3 = 0
    time_sorting = 0
    time_adding = 0
    for i_ele in range(num_ele):
        idx_ele = idx_node[i_ele]

        time_sorting0 = time.time()

        F_row_ptr = [0]
        K_row_ptr = [0]
        # row_ptr = np.array([0])
        for idx in idx_ele:
            gap = idx - len(F_row_ptr) + 1
            if gap > 0:
                # row_ptr = np.concatenate((row_ptr, [row_ptr[-1]] * gap))
                F_row_ptr = F_row_ptr + [F_row_ptr[-1]] * gap
                K_row_ptr = K_row_ptr + [K_row_ptr[-1]] * gap

            # row_ptr = np.concatenate((row_ptr, [row_ptr[-1] + 1]))
            F_row_ptr = F_row_ptr + [F_row_ptr[-1] + 1]
            K_row_ptr = K_row_ptr + [K_row_ptr[-1] + 3]

        # row_ptr = np.concatenate((row_ptr, [row_ptr[-1]] * (num_node - len(row_ptr) + 1)))
        F_row_ptr = F_row_ptr + [F_row_ptr[-1]] * (num_node - len(F_row_ptr) + 1)
        K_row_ptr = K_row_ptr + [K_row_ptr[-1]] * (num_node - len(K_row_ptr) + 1)

        time_sorting1 = time.time()

        F_row_ptr_arr = np.array(F_row_ptr, dtype=float, copy=False)
        F_col_idx = [0, 0, 0]
        F_ele_loc = F_loc[i_ele].flatten()

        K_row_ptr_arr = np.array(K_row_ptr, dtype=float, copy=False)
        K_col_idx = [idx_ele[0], idx_ele[1], idx_ele[2]] * 3
        K_ele_loc = K_loc[i_ele].flatten()

        time_sorting2 = time.time()

        F_ele = sp.sparse.csr_matrix((F_ele_loc, F_col_idx, F_row_ptr_arr), shape=(num_node, 1), dtype=float)
        K_ele = sp.sparse.csr_matrix((K_ele_loc, K_col_idx, K_row_ptr_arr), shape=(num_node, num_node), dtype=float)

        time_sorting3 = time.time()

        time_sum1 = time_sum1 + time_sorting1 - time_sorting0
        time_sum2 = time_sum2 + time_sorting2 - time_sorting1
        time_sum3 = time_sum3 + time_sorting3 - time_sorting2

        time_sorting = time_sorting + time_sorting3 - time_sorting0

        K = K + K_ele
        F = F + F_ele

        time_adding = time_adding + time.time() - time_sorting3

        # print(f"---- element: {i_ele: 6d}")

    time_end = time.time()
    print(
        f"time elapsed: {time_end - time_start:.6f}, total sorting: {time_sorting:.6f}, "
        f"sorting1: {time_sum1:.6f}, sorting2: {time_sum2:.6f}, sorting3: {time_sum3:.6f}, adding: {time_adding:.6f}")

    # print(f"global mat K: \n{K.toarray()}")
    # print(f"global mat F: \n{F.toarray()}")

    print(f"-------- incorporating 'essential' boundary conditions")
    for bc in boundary_conditions:
        if bc[1] == 0:
            K = K.tocsr(copy=False)
            F = F.tocsr(copy=False)
            if bc[0] == 0:
                K = K[1:, :]
                F = F[1:, :]
                F = F - bc[3] * K[:, 0]
                K = K[:, 1:]
            if bc[0] == 1:
                K = K[0: -1, :]
                F = F[0: -1, :]
                F = F - bc[3] * K[:, -1]
                K = K[:, 0: -1]

    # print(f"global mat K: \n{K}")
    # print(f"global mat F: \n{F}")

    print(f"-------- solving global matrix")
    u = sp.sparse.linalg.spsolve(K, F)
    for bc in boundary_conditions:
        if bc[1] == 0:
            if bc[0] == 0:
                u = np.insert(u, 0, bc[3], axis=0)
            if bc[0] == 1:
                u = np.insert(u, -1, bc[3], axis=0)

    print(f"u: \n{u}")

    plt.plot(x_node, u, linestyle='None', marker='o', color='blue')

    y_ext = lambda x: x + (math.exp(2 - x) + math.exp(x)) / (1 + math.exp(2))
    y_vfunc = np.vectorize(y_ext)
    x_test = np.arange(s0, s1 + 0.001, 0.01)
    y_test = y_vfunc(x_test)
    plt.plot(x_test, y_test, linestyle='-', marker='None', color='red')
    plt.show()
