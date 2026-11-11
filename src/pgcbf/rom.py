"""Structure-preserving Petrov-Galerkin ROM (Section sec:methodology).

Implements eq:pg_coords, eq:reduced_dynamics, eq:quad_split, Lemma lem:h2a_char and
Proposition prop:dahlquist quantities for the class-(H) eigenmode basis.
"""
from __future__ import annotations

import numpy as np

from .lattice import Lattice


class ROM:
    def __init__(self, lat: Lattice, K_m: int):
        self.lat = lat
        self.K_m = int(K_m)
        self.K = 2 * self.K_m
        Nm, Km = lat.Nm, self.K_m
        self.Phi_q = lat.modes[:, :Km]
        self.Phi_q_perp = lat.modes[:, Km:]
        Zq = np.zeros((Nm, Km))
        self.Phi = np.block([[self.Phi_q, Zq], [Zq, self.Phi_q]])
        Zt = np.zeros((Nm, Nm - Km))
        self.Phi_perp = np.block([[self.Phi_q_perp, Zt], [Zt, self.Phi_q_perp]])  # basis of V_perp
        # eq:pg_coords (general Petrov-Galerkin form, evaluated, not assumed)
        Pf = lat.Pfull
        self.P0 = self.Phi.T @ Pf @ self.Phi
        self.Zmap = np.linalg.solve(self.P0, self.Phi.T @ Pf)  # z(x) = Zmap x
        self.Pperp = np.eye(lat.N) - self.Phi @ self.Zmap
        # reduced operators
        self.Pe = self.Phi.T @ lat.Pfe @ self.Phi  # P := Phi^T P_eps Phi  (eq:quad_split)
        self.P = self.Pe
        P0inv = np.linalg.inv(self.P0)
        self.J0 = P0inv @ self.Phi.T @ Pf @ lat.J @ Pf @ self.Phi @ P0inv
        self.R0 = P0inv @ self.Phi.T @ Pf @ lat.R @ Pf @ self.Phi @ P0inv
        self.Ar0 = (self.J0 - self.R0) @ self.P0
        self.Arperp = self.Zmap @ lat.Af @ self.Pperp
        self.g_r = self.Zmap @ lat.G
        self.d_r = self.Zmap @ lat.D
        self.A_adv_r = self.Zmap @ lat.A_adv @ self.Phi  # linear part of F_nl,r
        # advective cross-coupling block A_{adv,r,perp} := P0^{-1} Phi^T P_full A_adv P_perp = Zmap A_adv P_perp
        # (Lemma lem:delta_cert_classH; linear part of the F_nl contribution to Delta_perp in eq:reduced_dynamics);
        # nonzero iff the advection band K_f exceeds K_m.  Review round: P_full (manuscript), not P_full,eps.
        self.A_adv_rperp = self.Zmap @ lat.A_adv @ self.Pperp
        self.Omega_r2 = self.Phi_q.T @ lat.Omega2 @ self.Phi_q
        self.omega_r = np.sqrt(np.diag(self.Omega_r2))

    # ------------------------------------------------------------ coordinates
    def z_of_x(self, X):
        return self.Zmap @ X

    def x_perp(self, X):
        return self.Pperp @ X

    def F_nl_r(self, Z):
        """Implements eq:reduced_dynamics: F_nl,r(z) = P0^{-1} Phi^T P_f F_nl(Phi z) (eq:reduced_dynamics)."""
        return self.Zmap @ self.lat.F_nl(self.Phi @ Z)

    def f_r(self, Z):
        return self.Ar0 @ Z + self.F_nl_r(Z)

    def hPhi(self, Z, E_safe):
        """Implements eq:quad_split (Lemma lem:quad_split): h(Phi z) = E_safe - z^T P z (eq:quad_split)."""
        if Z.ndim == 1:
            return E_safe - Z @ self.P @ Z
        return E_safe - np.einsum("ip,ij,jp->p", Z, self.P, Z)

    def Pnorm(self, Z):
        if Z.ndim == 1:
            return np.sqrt(Z @ self.P @ Z)
        return np.sqrt(np.einsum("ip,ij,jp->p", Z, self.P, Z))

    def tail_Pnorm(self, X):
        Xp = self.Pperp @ X
        if Xp.ndim == 1:
            return np.sqrt(Xp @ self.lat.Pfe @ Xp)
        return np.sqrt(np.einsum("ip,ij,jp->p", Xp, self.lat.Pfe, Xp))

    # ------------------------------------------------------------ hypotheses
    def h2a_residual(self) -> float:
        """Implements (H2a) of thm:main_cbf:  ||Phi^T P_eps P_perp||_2 (must vanish)."""
        return np.linalg.norm(self.Phi.T @ self.lat.Pfe @ self.Pperp, 2)

    def h2a_invariance_residual(self) -> float:
        """Implements Lemma lem:h2a_char: ||Omega^2 Phi_q - Phi_q S||_2."""
        S = self.Omega_r2
        return np.linalg.norm(self.lat.Omega2 @ self.Phi_q - self.Phi_q @ S, 2)

    def h2b_margin(self) -> float:
        """Implements (H2b) of thm:main_cbf:  min_{v in Im(I-Phi_q Phi_q^T)} v^T R_p v / ||v||^2 - eps (must be > 0)."""
        Vt = self.Phi_q_perp
        return np.linalg.eigvalsh(Vt.T @ self.lat.Rp @ Vt).min() - self.lat.eps
