import numpy as np

from Node import NodeDim1


# 2nd-order, 1-dim finite element
class FiniteElementOrd2Dim1(object):
    def __init__(self, nodes: list[NodeDim1], coef_functions: list[staticmethod],
                 quadrature_xi: np.array, quadrature_weights: np.array) -> None:
        self.nodes = nodes
        self.coef_vfunc_C = np.vectorize(coef_functions[0])
        self.coef_vfunc_beta = np.vectorize(coef_functions[1])
        self.coef_vfunc_alpha = np.vectorize(coef_functions[2])
        self.coef_vfunc_f = np.vectorize(coef_functions[3])
        self.quadrature_xi = quadrature_xi
        self.quadrature_xi_square = np.multiply(quadrature_xi, quadrature_xi)  # element-wise multiplication
        self.quadrature_weights = quadrature_weights

        # for physical domain to natural domain mapping shape function
        self.a_xi = None
        self.b_xi = None
        self.c_xi = None

        self.quadrature_x = None
        self.quadrature_x_square = None
        self.quadrature_jacobian = None

        # for coefficient functions
        self.quadrature_A_C = None
        self.quadrature_A_beta = None
        self.quadrature_A_alpha = None
        self.quadrature_f = None

        self.quadrature_MA_C = None
        self.quadrature_MA_beta = None
        self.quadrature_MA_alpha = None
        self.quadrature_Vb_f = None

        # for element shape function
        self.a_N = None
        self.b_N = None
        self.c_N = None

        # for coefficient matrix
        self.quadrature_N = None
        self.quadrature_NP = None

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

    def update_Nxi_NPxi_coef(self):
        x0, x1, x2 = self.x0x1x2()
        self.a_xi = 0.5 * (x0 + x2) - x1
        self.b_xi = 0.5 * (x2 - x0)
        self.c_xi = x1

    def update_quadrature_x_and_jac(self):
        self.quadrature_x = self.a_xi * self.quadrature_xi_square + self.b_xi * self.quadrature_xi + self.c_xi
        self.quadrature_x_square = np.multiply(self.quadrature_x, self.quadrature_x)
        self.quadrature_jacobian = 2 * self.a_xi * self.quadrature_xi + self.b_xi

    def update_quadrature_coef_function(self):
        self.quadrature_A_C = self.coef_vfunc_C(self.quadrature_x)
        self.quadrature_A_beta = self.coef_vfunc_beta(self.quadrature_x)
        self.quadrature_A_alpha = self.coef_vfunc_alpha(self.quadrature_x)
        self.quadrature_f = self.coef_vfunc_f(self.quadrature_x)

        self.quadrature_MA_C = np.diag(-1 * self.quadrature_A_C * self.quadrature_weights * self.quadrature_jacobian)
        self.quadrature_MA_beta = np.diag(self.quadrature_A_beta * self.quadrature_weights * self.quadrature_jacobian)
        self.quadrature_MA_alpha = np.diag(self.quadrature_A_alpha * self.quadrature_weights * self.quadrature_jacobian)
        self.quadrature_Vb_f = np.transpose([self.quadrature_f * self.quadrature_weights * self.quadrature_jacobian])

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
        self.quadrature_N = (np.outer(self.a_N, self.quadrature_x_square) +
                             np.outer(self.b_N, self.quadrature_x) +
                             np.outer(self.c_N, np.ones(self.quadrature_x.size)))
        self.quadrature_NP = (np.outer(2 * self.a_N, self.quadrature_x) +
                              np.outer(self.b_N, np.ones(self.quadrature_x.size)))

    def update_coefficient_matrix(self):
        self.AM_C = np.matmul(np.matmul(self.quadrature_NP, self.quadrature_MA_C), self.quadrature_NP.T)
        self.AM_beta = np.matmul(np.matmul(self.quadrature_N, self.quadrature_MA_beta), self.quadrature_NP.T)
        self.AM_alpha = np.matmul(np.matmul(self.quadrature_N, self.quadrature_MA_alpha), self.quadrature_N.T)

        self.bV = np.matmul(self.quadrature_N, self.quadrature_Vb_f)
