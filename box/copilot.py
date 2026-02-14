import numpy as np
from scipy.linalg import solveh_banded
import matplotlib.pyplot as plt


class Beam1D:
    """
    Euler-Bernoulli beam FEM with:
    - 2 DOF per node (w, theta)
    - full stiffness matrix for reactions
    - symmetric band matrix for solving
    - supports: fixed, simple
    Does not work, AI is not yet good enough for cases like this
    """

    def __init__(self, E, I, L, n_elem):
        self.E = E
        self.I = I
        self.L = L
        self.n_elem = n_elem

        self.x = np.linspace(0.0, L, n_elem + 1)
        self.Le = self.x[1] - self.x[0]

        self.ndof = 2 * (n_elem + 1)
        self.bw = 4  # half-bandwidth (max DOF offset = 3)

        self.K_full = np.zeros((self.ndof, self.ndof))
        self.F_full = np.zeros(self.ndof)

        self.fixed_dofs = set()

        self._build_element_matrices()

    def _build_element_matrices(self):
        L = self.Le
        E = self.E
        I = self.I

        self.ke = (E * I / L**3) * np.array([
            [12,     6*L,   -12,     6*L],
            [6*L,  4*L**2, -6*L,   2*L**2],
            [-12,   -6*L,    12,    -6*L],
            [6*L,  2*L**2, -6*L,   4*L**2]
        ])

        self.fe_uniform = (L / 12.0) * np.array([6, 3*L, 6, -3*L])

    # ------------------------------------------------------------
    # Loads
    # ------------------------------------------------------------

    def add_uniform_load(self, q, e):
        idx = np.array([2*e, 2*e+1, 2*e+2, 2*e+3])
        load = q * self.fe_uniform
        self.F_full[idx] += load

    def add_nodal_force(self, node, w=0.0, m=0.0):
        self.F_full[2*node] += w
        self.F_full[2*node+1] += m

    # ------------------------------------------------------------
    # Supports
    # ------------------------------------------------------------

    def fix_support(self, node):
        """Fixed support: w = 0, theta = 0"""
        self.fixed_dofs.add(2*node)
        self.fixed_dofs.add(2*node + 1)

    def simple_support(self, node):
        """Simple support: w = 0, theta free"""
        self.fixed_dofs.add(2*node)

    # ------------------------------------------------------------
    # Assembly
    # ------------------------------------------------------------

    def assemble_stiffness(self):
        self.K_full[:, :] = 0.0
        for e in range(self.n_elem):
            idx = np.array([2*e, 2*e+1, 2*e+2, 2*e+3])
            for a in range(4):
                for b in range(4):
                    self.K_full[idx[a], idx[b]] += self.ke[a, b]

    def apply_boundary_conditions(self):
        for dof in self.fixed_dofs:
            self.K_full[dof, :] = 0.0
            self.K_full[:, dof] = 0.0
            self.K_full[dof, dof] = 1.0
            self.F_full[dof] = 0.0

    # ------------------------------------------------------------
    # Band matrix builder
    # ------------------------------------------------------------

    def build_band_matrix(self):
        Kb = np.zeros((self.bw, self.ndof))
        for i in range(self.ndof):
            for j in range(i, min(i + self.bw, self.ndof)):
                Kb[j - i, i] = self.K_full[i, j]
        return Kb

    # ------------------------------------------------------------
    # Solve
    # ------------------------------------------------------------

    def solve(self):
        self.assemble_stiffness()
        self.apply_boundary_conditions()
        Kb = self.build_band_matrix()
        self.U = solveh_banded(Kb, self.F_full, lower=False)
        return self.U

    # ------------------------------------------------------------
    # Post-processing
    # ------------------------------------------------------------

    def get_displacement(self):
        return self.x, self.U[0::2]

    def get_rotation(self):
        return self.x, self.U[1::2]

    def get_reactions(self):
        R = self.K_full @ self.U - self.F_full
        return {dof: R[dof] for dof in sorted(self.fixed_dofs)}

if __name__ == "__main__":
    E = 210e9
    I = 5e-6
    L = 3.0
    a = 0.5
    q = 1000.0
    n_elem = 80

    beam = Beam1D(E, I, L, n_elem)

    # uniform load
    for e in range(n_elem):
        beam.add_uniform_load(q, e)

    # fixed support at x=0
    beam.fix_support(0)

    # simple support at x=a
    i = np.argmin(np.abs(beam.x - a))
    beam.simple_support(i)

    U = beam.solve()
    x, w = beam.get_displacement()
    reactions = beam.get_reactions()

    print("Support reactions:")
    for dof, R in reactions.items():
        print(f"  DOF {dof:3d}: {R: .3e}")

    plt.plot(x, w)
    plt.axvline(a, color='k', ls='--')
    plt.xlabel("x [m]")
    plt.ylabel("w [m]")
    plt.title("Beam deflection with fixed + simple support")
    plt.show()