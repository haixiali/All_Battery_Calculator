import numpy as np

from Node import Node


# 2nd-order, 1-dim finite element
class FiniteElementDim1Ord2(object):
    def __init__(self, nodes: list[Node], coef_functions: list[staticmethod],
                 quadrature_xi: np.array, quadrature_weights: np.array) -> None:
        self.nodes = nodes

        self.coef_vfunc_C = np.vectorize(coef_functions[0])
        self.coef_vfunc_beta = np.vectorize(coef_functions[1])
        self.coef_vfunc_alpha = np.vectorize(coef_functions[2])
        self.coef_vfunc_f = np.vectorize(coef_functions[3])

        self.quad_xi = quadrature_xi
        self.quad_xi_square = np.multiply(quadrature_xi, quadrature_xi)  # element-wise multiplication
        self.quad_weights = quadrature_weights

        # for physical domain to natural domain mapping shape function
        self.a_xi = None
        self.b_xi = None
        self.c_xi = None

        self.quad_x = None
        self.quad_x_square = None
        self.quad_jac = None

        # for coefficient functions
        self.quad_A_C_Vect = None
        self.quad_A_beta_Vect = None
        self.quad_A_alpha_Vect = None
        self.quad_f_Vect = None

        self.quad_A_C = None
        self.quad_A_beta = None
        self.quad_A_alpha = None
        self.quad_b_f = None

        # for element shape function
        self.a_N = None
        self.b_N = None
        self.c_N = None

        # for coefficient matrix
        self.N = None
        self.NP = None

        # local matrix
        self.AM_C = None
        self.AM_beta = None
        self.AM_alpha = None
        # right-hand side b vector
        self.bV = None

        self.R = None
        self.S = None

    def x0x1x2(self):
        return self.nodes[0].x, self.nodes[1].x, self.nodes[2].x

    """
    Shape function mapping physical space to quadrature space.
    """

    def update_Nxi_NPxi_coef(self):
        x0, x1, x2 = self.x0x1x2()
        self.a_xi = 0.5 * (x0 + x2) - x1
        self.b_xi = 0.5 * (x2 - x0)
        self.c_xi = x1

    """
    Calculate quadrature coordinates on physical space and Jacobian coefficients.
    """

    def update_quadrature_x_and_jacobian(self):
        self.quad_x = self.a_xi * self.quad_xi_square + self.b_xi * self.quad_xi + self.c_xi
        self.quad_x_square = np.multiply(self.quad_x, self.quad_x)
        self.quad_jac = 2 * self.a_xi * self.quad_xi + self.b_xi

    def update_quadrature_coef_function(self):
        self.quad_A_C_Vect = self.coef_vfunc_C(self.quad_x)
        self.quad_A_beta_Vect = self.coef_vfunc_beta(self.quad_x)
        self.quad_A_alpha_Vect = self.coef_vfunc_alpha(self.quad_x)
        self.quad_f_Vect = self.coef_vfunc_f(self.quad_x)

        self.quad_A_C = np.diag(-1 * self.quad_A_C_Vect * self.quad_weights * self.quad_jac)
        self.quad_A_beta = np.diag(self.quad_A_beta_Vect * self.quad_weights * self.quad_jac)
        self.quad_A_alpha = np.diag(self.quad_A_alpha_Vect * self.quad_weights * self.quad_jac)
        self.quad_b_f = np.transpose([self.quad_f_Vect * self.quad_weights * self.quad_jac])

    def update_N_and_NP_coef(self):
        x0, x1, x2 = self.x0x1x2()
        z01 = x0 - x1
        z12 = x1 - x2
        z20 = x2 - x0
        alpha = -1 / (z01 * z12 * z20)
        alpha_vec = alpha * np.array([z12, z20, z01])

        self.a_N = np.transpose([alpha_vec])
        self.b_N = -1 * self.a_N * np.array([[x1 + x2, x0 + x2, x0 + x1]]).T
        self.c_N = self.a_N * np.array([[x1 * x2, x0 * x2, x0 * x1]]).T

    def update_quadrature_N_and_NP(self):
        self.N = (np.outer(self.a_N, self.quad_x_square) +
                  np.outer(self.b_N, self.quad_x) +
                  np.outer(self.c_N, np.ones(self.quad_x.size)))
        self.NP = np.outer(2 * self.a_N, self.quad_x) + np.outer(self.b_N, np.ones(self.quad_x.size))

    def update_coefficient_matrix(self):
        self.AM_C = np.matmul(np.matmul(self.NP, self.quad_A_C), self.NP.T)
        self.AM_beta = np.matmul(np.matmul(self.N, self.quad_A_beta), self.NP.T)
        self.AM_alpha = np.matmul(np.matmul(self.N, self.quad_A_alpha), self.N.T)

        self.bV = np.matmul(self.N, self.quad_b_f)

    def update_boundary(self):
        pass
