"""Full-order class-(H) telegraph reaction-diffusion-advection (RDA) lattice.

Plant (eq:full_system with Sigma = 0 unless stochastic):
    q' = p
    p' = -(nu L_h + alpha I) q - gamma p - c Pi_f D_h q - kappa (q - sin q) + B u + B_w w

* Omega_full^2 = nu L_h + alpha I  (L_h: Dirichlet FD Laplacian, mass-normalised)
* J_full canonical, R_full = blkdiag(0, gamma I)  (class (H), R_p = gamma I)
* F_nl(x) = A_adv x + F_sin(x),  A_adv = [[0, 0], [-c W_f D_h W_f, 0]],  F_sin(x) = [0 ; -kappa (q - sin q)]
  (Assumption assump:generic).  A_full = (J_full - R_full) P_full stays the pure damped-wave operator.
  D_h is the skew central-difference advection operator and W_f = Phi_f diag(w) Phi_f^T the graded
  band-limiting filter on the first K_f = n_adv_modes lattice sine modes (w = adv_mode_weights; w = 1 gives
  the orthogonal projector Pi_f).  With K_f > K_m the advection couples the ROM span to the tail
  (A_{adv,r,perp} != 0); see design_note.md for why the unfiltered lattice advection is uncertifiable under an
  L2-ball ODD (Remark rem:rda_cutoff) and why the cutoff modes carry a reduced weight.

All arrays are float64.  No numerical constant other than the configuration
inputs appears here.
"""
from __future__ import annotations

import numpy as np


def dirichlet_laplacian(n: int, h: float) -> np.ndarray:
    """Second-order central FD Laplacian (positive definite), Remark rem:tail_rates(ii)."""
    L = np.zeros((n, n))
    idx = np.arange(n)
    L[idx, idx] = 2 / h**2
    L[idx[:-1], idx[:-1] + 1] = -1 / h**2
    L[idx[:-1] + 1, idx[:-1]] = -1 / h**2
    return L


def central_difference(n: int, h: float) -> np.ndarray:
    """Skew-symmetric central difference d/dxi with homogeneous Dirichlet data."""
    Dm = np.zeros((n, n))
    idx = np.arange(n - 1)
    Dm[idx, idx + 1] = 1 / (2 * h)
    Dm[idx + 1, idx] = -1 / (2 * h)
    return Dm


def gaussian_profiles(xi: np.ndarray, centers, width: float, h: float) -> np.ndarray:
    """Spatially smooth distributed input profiles, columns scaled by sqrt(h) so that
    the Euclidean lattice norm approximates the L2(0,ell) norm (resolution independent)."""
    cols = []
    for c0 in centers:
        b = np.exp(-((xi - c0) ** 2) / (2 * width**2))
        b = b / np.sqrt(h * np.sum(b**2))  # unit L2 norm of the profile
        cols.append(np.sqrt(h) * b)
    return np.stack(cols, axis=1)


class Lattice:
    """Class-(H) momentum-damped RDA lattice with all full-order operators."""

    def __init__(self, cfg: dict):
        p = cfg["plant"]
        self.cfg = cfg
        self.Nm = int(p["N_m"])
        self.N = 2 * self.Nm
        self.ell = float(p["ell"])
        self.nu = float(p["nu"])
        self.alpha = float(p["alpha"])
        self.gamma = float(p["gamma"])
        self.c = float(p["c"])
        self.kappa = float(p["kappa"])
        self.n_f = int(p["n_adv_modes"])
        self.h = self.ell / (self.Nm + 1)
        self.xi = self.h * np.arange(1, self.Nm + 1)
        self.eps = self.gamma / 2  # eq:eps_choice
        Nm, h = self.Nm, self.h
        self.Lh = dirichlet_laplacian(Nm, h)
        self.Dh = central_difference(Nm, h)
        self.Omega2 = self.nu * self.Lh + self.alpha * np.eye(Nm)
        # exact eigen-decomposition of the mass-normalised stiffness (sorted ascending)
        lam, V = np.linalg.eigh(self.Omega2)
        # fix the sign convention of each mode deterministically (first nonzero entry > 0)
        sgn = np.sign(V[0])  # sine modes never vanish at the first node
        self.modes = V * sgn
        self.omega2 = lam
        self.omega = np.sqrt(lam)
        self.omega_min = self.omega[0]  # smallest full-order natural frequency
        I, Z = np.eye(Nm), np.zeros((Nm, Nm))
        self.J = np.block([[Z, I], [-I, Z]])
        self.Rp = self.gamma * I
        self.R = np.block([[Z, Z], [Z, self.Rp]])
        self.Pfull = np.block([[self.Omega2, Z], [Z, I]])
        self.Pfe = np.block([[self.Omega2, self.eps * I], [self.eps * I, I]])  # eq:eps_choice
        self.Af = (self.J - self.R) @ self.Pfull
        # band-limited advection (linear part of F_nl)
        Phi_f = self.modes[:, : self.n_f]
        self.Pi_f = Phi_f @ Phi_f.T
        # graded band-limiting weights W_f = Phi_f diag(w) Phi_f^T: weight 1 on the first n_adv_full_weight modes,
        # the disclosed cutoff-mode weight w = adv_cutoff_weight on modes n_adv_full_weight+1 .. K_f
        # (w = 1: plain orthogonal projector Pi_f)
        nfull = int(p.get("n_adv_full_weight", self.n_f))
        wc = float(p.get("adv_cutoff_weight", 1))
        w = np.where(np.arange(self.n_f) < nfull, 1, wc).astype(float)
        self.adv_w = w
        self.adv_cutoff_weight = wc
        self.W_f = (Phi_f * w) @ Phi_f.T
        # band-limited large-scale advection of the large-scale field: Pi_f D_h Pi_f
        self.Cadv = -self.c * self.W_f @ self.Dh @ self.W_f  # acts on q, enters the p-equation
        self.A_adv = np.block([[Z, Z], [self.Cadv, Z]])
        self.A_lin = self.Af + self.A_adv
        # inputs
        # smooth, spectrally band-limited actuator / disturbance profiles (Gaussian patches
        # passed through a Butterworth roll-off in the sine-mode index), unit L2 norm
        self.B = self._smooth(gaussian_profiles(self.xi, p["actuator_centers"], float(p["actuator_width"]), h),
                              p["actuator_band"])
        self.Bw = self._smooth(gaussian_profiles(self.xi, p["disturbance_centers"], float(p["disturbance_width"]), h),
                               p["disturbance_band"])
        self.m = self.B.shape[1]
        self.nw = self.Bw.shape[1]
        self.G = np.vstack([np.zeros((Nm, self.m)), self.B])
        self.D = np.vstack([np.zeros((Nm, self.nw)), self.Bw])

    def _smooth(self, Bc, band):
        """Butterworth spectral roll-off f_k = (1 + (k/k_c)^{2p})^{-1/2} in the lattice sine-mode
        index k (band = [k_c, p]); profiles renormalised to unit L2 norm."""
        kc, order = float(band[0]), int(band[1])
        k = np.arange(1, self.Nm + 1)
        f = 1 / np.sqrt(1 + (k / kc) ** (2 * order))
        out = (self.modes * f) @ (self.modes.T @ Bc)
        return out / np.linalg.norm(out, axis=0) * np.linalg.norm(Bc, axis=0)

    # ------------------------------------------------------------------ physics
    @property
    def peclet(self) -> float:
        """Implements Remark rem:rda_cutoff: Pe := c ell / nu (Remark rem:rda_cutoff)."""
        return self.c * self.ell / self.nu

    def nonlinear(self, X: np.ndarray) -> np.ndarray:
        """Genuinely nonlinear part of F_nl: [0; -kappa (q - sin q)] (X is N or N x P)."""
        q = X[: self.Nm]
        out = np.zeros_like(X)
        out[self.Nm :] = -self.kappa * (q - np.sin(q))
        return out

    def F_nl(self, X: np.ndarray) -> np.ndarray:
        """Full residual F_nl(x) of Assumption assump:generic (F_nl(0) = 0)."""
        return self.A_adv @ X + self.nonlinear(X)

    def F(self, X: np.ndarray) -> np.ndarray:
        """Drift F(x) = A_full x + F_nl(x)."""
        return self.Af @ X + self.F_nl(X)

    def h_barrier(self, X: np.ndarray, E_safe: float) -> np.ndarray:
        """Implements eq:h_def  h(x) = E_safe - x^T P_full,eps x (vectorised over columns)."""
        if X.ndim == 1:
            return E_safe - X @ self.Pfe @ X
        return E_safe - np.einsum("ip,ij,jp->p", X, self.Pfe, X)
