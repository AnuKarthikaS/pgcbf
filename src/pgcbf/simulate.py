"""Sampled-data closed loop with exact multi-rate latch timing (eq:multirate_hold).

Full-order plant x' = A x + n(x) + G u + D w, A = A_full + A_adv (linear part of F_nl),
n(x) the genuinely nonlinear residual.  Integrated with the second-order exponential
time-differencing scheme ETD2RK (Cox & Matthews); the stiff linear part is propagated
exactly by matrix exponentials, so the step is not limited by lambda_N.

Each sampling period [t_k, t_k + T_s) is split *exactly* at tau_A and tau_B into three
segments; the integrator never steps across a latch instant, so the control switch is
exact (independent of the integrator resolution):
    u = u_{B,k-1} on [t_k, t_k+tau_A),  u_{A,k} on [t_k+tau_A, t_k+tau_B),  u_{B,k} on [t_k+tau_B, t_{k+1}),
with u_{B,-1} = Psi_B(z_0).
"""
from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.linalg import expm


def phi_matrices(A: np.ndarray, h: float):
    """Exponential (ETD) integrator kernels, numerical method of the validation suite: exp(Ah), phi_1(Ah), phi_2(Ah) from one augmented exponential (exact, no cancellation):
    expm([[Ah, I, 0], [0, 0, I], [0, 0, 0]]) = [[e^{Ah}, phi_1(Ah), phi_2(Ah)], ...]."""
    n = A.shape[0]
    M = np.zeros((3 * n, 3 * n))
    M[:n, :n] = A * h
    M[:n, n:2 * n] = np.eye(n)
    M[n:2 * n, 2 * n:] = np.eye(n)
    Ex = expm(M)
    return Ex[:n, :n], Ex[:n, n:2 * n], Ex[:n, 2 * n:]


def van_loan_cov(A: np.ndarray, S: np.ndarray, h: float, rank_rel: float):
    """Exact covariance C = int_0^h e^{As} S S^T e^{A^T s} ds (Van Loan); returns L with L L^T = C."""
    n = A.shape[0]
    M = np.zeros((2 * n, 2 * n))
    M[:n, :n] = -A * h
    M[:n, n:] = S @ S.T * h
    M[n:, n:] = A.T * h
    Ex = expm(M)
    C = Ex[n:, n:].T @ Ex[:n, n:]
    C = (C + C.T) / 2
    w, V = np.linalg.eigh(C)
    keep = w > rank_rel * w.max()  # drop the numerically-null (roundoff) directions
    return V[:, keep] * np.sqrt(w[keep])


def _clean(M, tol):
    """Drop entries below tol * max|M| (roundoff of the orthogonal change of coordinates)."""
    M = np.array(M, copy=True)
    M[np.abs(M) <= tol * np.abs(M).max()] = 0
    return M


def _sparse(M):
    """Exact sparse copy (only exact zeros are dropped)."""
    return sparse.csr_matrix(M)


class Simulator:
    """Full-order sampled-data simulator.  Internally uses the orthonormal lattice eigen-
    coordinates y = T^T x, T = blkdiag(V, V) (an exact change of coordinates of the full
    N-dimensional state), in which exp(Ah), phi_1, phi_2 are exactly sparse."""

    def __init__(self, lat, rom, policy, cfg, tau_A=None, tau_B=None, T_s=None, substeps=None,
                 adversary="worst", eps_adv=None, controller=True, noise=None, E_safe=None, scheme="etd2"):
        self.lat, self.rom, self.pol = lat, rom, policy
        tm = cfg["timing"]
        self.Ts = float(tm["T_s"] if T_s is None else T_s)
        self.tA = float(tm["tau_A"] if tau_A is None else tau_A)
        self.tB = float(tm["tau_B"] if tau_B is None else tau_B)
        self.sub = [int(s) for s in (cfg["sim"]["substeps"] if substeps is None else substeps)]
        self.adv = adversary
        self.eps_adv = float(cfg["safety"]["eps_adv"] if eps_adv is None else eps_adv)
        self.controller = controller
        self.E = float(cfg["safety"]["E_safe"] if E_safe is None else E_safe)
        self.scheme = scheme
        Nm = lat.Nm
        V = lat.modes
        T = np.zeros((lat.N, lat.N))
        T[:Nm, :Nm] = V
        T[Nm:, Nm:] = V
        self.T = T
        tol = float(cfg["tolerances"]["modal_clean"])
        self.A = _clean(T.T @ lat.A_lin @ T, tol)  # block pattern: roundoff-level entries dropped
        self.G = T.T @ lat.G
        self.Dm = T.T @ lat.D
        self.Pfe = _sparse(_clean(T.T @ lat.Pfe @ T, tol))
        self.Zmap = _clean(rom.Zmap @ T, tol)
        self.Pperp = _sparse(_clean(T.T @ rom.Pperp @ T, tol))
        self.V = V
        Sig = None if noise is None else T.T @ noise
        seg = [self.tA, self.tB - self.tA, self.Ts - self.tB]
        self.segments = []
        cache = {}
        for i, (L, ns) in enumerate(zip(seg, self.sub)):
            if L <= 0:
                continue
            h = L / ns
            key = (i, ns)
            if h not in [c[0] for c in cache.values()]:
                E, p1, p2 = phi_matrices(self.A, h)
                Lc = van_loan_cov(self.A, Sig, h, float(cfg["tolerances"]["cov_rank_rel"])) if Sig is not None else None
                cache[key] = (h, _sparse(E), _sparse(p1), _sparse(p2), Lc)
            else:
                cache[key] = next(c for c in cache.values() if c[0] == h)
            self.segments.append((i, ns) + cache[key])
        self.Pd = rom.P @ rom.d_r

    # ------------------------------------------------------------ forcing
    def w_star(self, Y):
        """Implements Remark rem:hji_game: w* = -eps_adv L_dr h~^T/||.|| = eps_adv d_r^T P z/||d_r^T P z||."""
        if self.adv is None or self.eps_adv == 0:
            return np.zeros((self.lat.nw, Y.shape[1]))
        v = self.Pd.T @ (self.Zmap @ Y)
        nrm = np.linalg.norm(v, axis=0)
        out = np.zeros_like(v)
        nz = nrm > 0
        out[:, nz] = self.eps_adv * v[:, nz] / nrm[nz]
        return out

    def nonlin(self, Y):
        """Modal coordinates of n(x) = [0; -kappa (q - sin q)] (nodal evaluation, exact transforms)."""
        out = np.zeros_like(Y)
        if self.lat.kappa != 0:
            q = self.V @ Y[: self.lat.Nm]
            out[self.lat.Nm:] = self.V.T @ (-self.lat.kappa * (q - np.sin(q)))
        return out

    def N(self, Y, U):
        return self.nonlin(Y) + self.G @ U + self.Dm @ self.w_star(Y)

    def h_full(self, Y):
        return self.E - np.sum(Y * (self.Pfe @ Y), axis=0)

    def tail(self, Y):
        Yp = self.Pperp @ Y
        return np.sqrt(np.sum(Yp * (self.Pfe @ Yp), axis=0))

    # ------------------------------------------------------------ main loop
    def run(self, X0, n_periods, rng=None, record="substep", stop_on_exit=False, node_hook=None):
        """record: 'substep' (every integrator node incl. latch instants) or 'period'.

        node_hook (Experiment 4b, read-only): if given, called before every integrator step as
        node_hook(k, i, s, Y, U, Z_k, W) with k the period, i the latch segment (0: [0, tau_A), 1: [tau_A, tau_B),
        2: [tau_B, T_s)), s = t - t_k the offset of the step's left end point, Y the modal state, U the applied
        (right-continuous) input, Z_k the latched sample z(t_k) and W = w*(Y).  It does not alter the integration."""
        lat, rom, pol = self.lat, self.rom, self.pol
        Y = self.T.T @ np.array(X0, dtype=float)
        P = Y.shape[1]
        zeroU = np.zeros((lat.m, P))
        uBprev = pol.psi_B(self.Zmap @ Y) if self.controller else zeroU
        ts, hs, hts, tails, us = [0 * self.Ts], [self.h_full(Y)], [rom.hPhi(self.Zmap @ Y, self.E)], [self.tail(Y)], []
        z_s, uA_h, uB_h = [], [], []
        out_of_odd = np.zeros(P, bool)
        for k in range(n_periods):
            Z = self.Zmap @ Y
            out_of_odd |= np.any(np.abs(Z) > pol.D_ODD / 2, axis=0)
            if self.controller:
                uA, uB = pol.psi_A(Z), pol.psi_B(Z)
            else:
                uA = uB = zeroU
            z_s.append(Z); uA_h.append(uA); uB_h.append(uB)
            inputs = (uBprev, uA, uB)
            tseg = k * self.Ts
            soff = 0 * self.Ts
            for (i, ns, h, E, p1, p2, Lc) in self.segments:
                U = inputs[i]
                for s in range(ns):
                    if node_hook is not None:
                        node_hook(k, i, soff + s * h, Y, U, Z, self.w_star(Y))
                    Nx = self.N(Y, U)
                    Aa = E @ Y + p1 @ (h * Nx)
                    if Lc is not None:  # exponential Euler-Maruyama with exact OU covariance
                        Y = Aa + Lc @ rng.standard_normal((Lc.shape[1], P))
                    elif self.scheme == "etd1":
                        Y = Aa
                    else:  # ETD2RK
                        Y = Aa + p2 @ (h * (self.N(Aa, U) - Nx))
                    if record == "substep":
                        ts.append(tseg + (s + 1) * h)
                        hs.append(self.h_full(Y)); hts.append(rom.hPhi(self.Zmap @ Y, self.E))
                        tails.append(self.tail(Y)); us.append(U)
                tseg += ns * h
                soff = (self.tA, self.tB, self.Ts)[i]  # end of segment i = offset of the next one
            uBprev = uB
            if record == "period":
                ts.append((k + 1) * self.Ts)
                hs.append(self.h_full(Y)); hts.append(rom.hPhi(self.Zmap @ Y, self.E))
                tails.append(self.tail(Y)); us.append(uB)
            if stop_on_exit and np.all(np.min(np.array(hs), axis=0) < 0):
                break
        return {"t": np.array(ts), "h": np.array(hs), "hPhi": np.array(hts), "tail": np.array(tails),
                "u": np.array(us), "z_samples": np.array(z_s), "uA": np.array(uA_h), "uB": np.array(uB_h),
                "X_final": self.T @ Y, "out_of_odd": out_of_odd}


def initial_conditions(lat, rom, cert, n, rng, level_fraction, tail_fraction, n_tail_modes, z_norm2=None):
    """Implements (H7) / eq:eps_basis_def: initial states: z(0) with ||z(0)||_P^2 = level_fraction*(E - eps_b^2 - dh), i.e. on
    the level h~_r = (1 - level_fraction)*(E - eps_b^2 - dh) of C~_r (level_fraction = 1: boundary of C~_r), or with
    ||z(0)||_P^2 = z_norm2 if given (S6: z_norm2 = E - h_act/2 puts z(0) on h_r = h_act/2); random P-direction; and a
    random tail x_perp(0) in the first n_tail_modes tail modes with ||x_perp(0)||_Pfe = tail_fraction * M_perp/lambda_{K+1}."""
    E, eb, dh = cert["E_safe"], cert["eps_basis_P"], cert["Delta_h_dip"]
    lev = (E - eb**2 - dh) * level_fraction if z_norm2 is None else z_norm2
    from .certificates import psd_isqrt
    G = rng.standard_normal((rom.K, n))
    Z = psd_isqrt(rom.P) @ G
    Z = Z / rom.Pnorm(Z) * np.sqrt(lev)
    Bp = rom.Phi_perp
    nt = rom.lat.Nm - rom.K_m
    sel = np.r_[np.arange(n_tail_modes), nt + np.arange(n_tail_modes)]
    Y = np.zeros((Bp.shape[1], n))
    Y[sel] = rng.standard_normal((sel.size, n))
    Xp = Bp @ Y
    nrm = np.sqrt(np.einsum("ip,ij,jp->p", Xp, lat.Pfe, Xp))
    Xp = Xp / nrm * (tail_fraction * cert["tail_init_bound"])
    return rom.Phi @ Z + Xp
