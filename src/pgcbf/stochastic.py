"""Corollary cor:stochastic_w2 (Ito sampled-data safety, Wasserstein-2 horizon) for
additive, state-independent diffusion Sigma (so Sigma_r, Sigma_perp are constant).

Implements eq:delta_stoch, eq:T_TTL and the OU tail envelope eps_basis,sigma.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.stats import beta as beta_dist

from .certificates import psd_sqrt, opnorm, exp_weights, phi_hold, two_branch


def noise_matrix(lat, rom, ncfg: dict) -> np.ndarray:
    """Diffusion Sigma (N x n_noise) acting on the momentum equation.

    type 'profile': sigma * D (smooth, band-limited profiles of the disturbance channel)
    type 'tail_white': sigma * [0; I - Phi_q Phi_q^T] (spatially white small-scale forcing,
                       orthogonal to the ROM span; Sigma_r = 0 exactly)
    type 'white': sigma * [0; I]."""
    s = float(ncfg["sigma"])
    Nm = lat.Nm
    Z = np.zeros((Nm, Nm))
    if ncfg["type"] == "profile":
        return s * lat.D
    if ncfg["type"] == "tail_white":
        return s * np.vstack([Z, np.eye(Nm) - rom.Phi_q @ rom.Phi_q.T])
    if ncfg["type"] == "white":
        return s * np.vstack([Z, np.eye(Nm)])
    raise ValueError(ncfg["type"])


def corollary(cfg: dict, cert: dict, sysd: dict, Sigma: np.ndarray, T: float, alpha: float,
              ignore_tail_diffusion: bool = False) -> dict:
    """Implements eq:delta_stoch and eq:T_TTL (Corollary cor:stochastic_w2). Evaluate every quantity of Corollary cor:stochastic_w2 for horizon T and level alpha.

    ignore_tail_diffusion=True gives the 'crude' bound (sigma_perp := 0) used as a foil."""
    lat, rom = sysd["lat"], sysd["rom"]
    sc = cfg["stochastic"]
    E, Ts, g0 = cert["E_safe"], cert["T_s"], cert["gamma0"]
    tA, tB = cert["tau_A"], cert["tau_B"]
    lam, L_h, g_h, d_h = cert["lambda_K1"], cert["L_h"], cert["g_h"], cert["d_h"]
    L_perp, eC, eA, L_Tz = cert["L_perp"], cert["eps_CROWN"], cert["eps_A"], cert["L_Tz"]
    u_max, L_Lfh, L_Lgh = cert["u_max"], cert["L_Lfh"], cert["L_Lgh"]
    K = rom.K
    P = rom.P
    Ph = psd_sqrt(P)
    r_W2 = float(sc["r_W2"])
    Sr = rom.Zmap @ Sigma
    sig_r = opnorm(Sr)
    Pfeh = psd_sqrt(lat.Pfe)
    sig_perp2 = 0 if ignore_tail_diffusion else float(np.linalg.norm(Pfeh @ rom.Pperp @ Sigma, "fro") ** 2)
    out = {"sigma_r": sig_r, "sigma_perp2": sig_perp2, "T": T, "alpha": alpha}

    def consts(Tp):
        """Drift/diffusion-enlarged constants at horizon Tp (all replacements of the Corollary)."""
        eps_adv_T = cert_eps_adv + r_W2 * Tp
        M_perp_T = cert["M_perp"] + cert["PperpD_norm"] * r_W2 * Tp
        eb_T = M_perp_T / lam
        dsig = (sig_perp2 ** (1 / 2)) * math.sqrt(2 * math.log(4 * (1 + lam * Tp / 2) / alpha) / lam) if sig_perp2 > 0 else 0
        eb_s = eb_T + dsig
        Mz_T = cert["M_z"] + opnorm(rom.d_r) * r_W2 * Tp + L_perp * (eb_s - cert["eps_basis_P"])
        du_lag = L_Tz * Mz_T * Ts
        L_hold = (L_Lfh + L_Lgh * u_max + g0 * L_h) * Mz_T
        dh = g_h * (eA - eC) * (tB - tA) + g_h * du_lag * tA + L_hold * Ts**2 / 2
        return dict(eps_adv_T=eps_adv_T, M_perp_T=M_perp_T, eb_T=eb_T, Delta_sigma=dsig, eb_sigma=eb_s,
                    M_z_T=Mz_T, du_lag_T=du_lag, L_hold_T=L_hold, dh_T=dh)

    cert_eps_adv = float(cfg["safety"]["eps_adv"])
    c = consts(T)
    out.update(c)
    # T_perp: sup{T' : eb_sigma(T')^2 + dh(T') + L_h M_z(T') T_s < E}
    def ok(Tp):
        cc = consts(Tp)
        return cc["eb_sigma"] ** 2 + cc["dh_T"] + L_h * cc["M_z_T"] * Ts < E
    if not ok(Ts):
        T_perp = 0
    else:
        lo, hi = Ts, Ts
        while ok(hi) and hi < float(sc["T_perp_search_max"]):
            lo, hi = hi, 2 * hi
        if ok(hi):
            T_perp = math.inf
        else:
            for _ in range(int(sc["bisection_iters"])):
                mid = (lo + hi) / 2
                lo, hi = (mid, hi) if ok(mid) else (lo, mid)
            T_perp = lo
    out["T_perp"] = T_perp
    # diffusion constants for constant Sigma_r
    d_ito = float(np.linalg.norm(Ph @ Sr, "fro") ** 2)
    sigma_h = 2 * math.sqrt(E) * opnorm(Ph @ Sr)
    N_T = max(1, int(round(T / Ts)))
    r_BM = sig_r * math.sqrt(2 * K * Ts * math.log(2**3 * K * N_T / alpha))
    out.update(delta_Ito=d_ito, sigma_h=sigma_h, r_BM=r_BM, N_T=N_T)
    # replaced shell / certificate (eq:shell, eq:teacher_cbf with eps_basis,sigma)
    w1, w2, w3 = exp_weights(g0, Ts, tA, tB)
    phi = phi_hold(g0 * Ts, float(cfg["tolerances"]["phi_series_below"]))
    room = E - c["eb_sigma"] ** 2 - c["dh_T"]
    LMT = L_h * c["M_z_T"] * Ts
    h_act = LMT + float(cfg["safety"]["h_act_fraction"]) * (room - LMT)
    r2 = room - h_act
    out.update(h_act_sigma=h_act, r_min2_sigma=r2)
    h_act_minus = h_act - L_h * (c["M_z_T"] * Ts + r_BM)
    out["h_act_minus"] = h_act_minus
    dc = two_branch(cert["c_M"], r2, h_act, c["dh_T"], g0)
    out["delta_cert_sigma"] = dc["delta_cert"]
    terms = {"tail": L_h * L_perp * c["eb_sigma"], "quant": g_h * (w2 * eA + (1 - w2) * eC),
             "lag": w1 * g_h * c["du_lag_T"], "adv": d_h * cert_eps_adv, "hold": phi * c["L_hold_T"] * Ts}
    dts = sum(terms.values())
    out["delta_total_sigma"] = dts
    # eq:delta_stoch
    dpw = (g_h * max(eA, eC + c["du_lag_T"]) - g_h * (w2 * eA + (1 - w2) * eC + w1 * c["du_lag_T"])
           + (1 - phi) * c["L_hold_T"] * Ts + (L_Lfh + L_Lgh * u_max + g0 * L_h + g_h * L_Tz) * r_BM)
    h0 = min(float(sc["h0_fraction"]) * h_act_minus, h_act) if h_act_minus > 0 else float(sc["h0_fraction"]) * h_act
    d_stoch = dpw + g0 * h_act + d_ito + math.log(2 / alpha) * sigma_h**2 / (2 * h0)
    d_drift = float(sc["delta_drift_fraction"]) * dc["delta_cert"] if dc["delta_cert"] > 0 else float(sc["delta_drift_fraction"])
    out.update(Delta_delta_pw=dpw, h0=h0, delta_stoch=d_stoch, delta_drift=d_drift)
    # Kushner exponent, interior drift bound (sound: ||grad h~|| <= L_h, ||f+gu+dw+Delta_perp|| <= M_z(T))
    theta = math.log(2 / alpha) / h0
    mu_int = L_h * c["M_z_T"]
    c_int = theta * (mu_int + d_ito + theta * sigma_h**2 / 2) * math.exp(-theta * (h_act_minus + c["dh_T"]))
    beta0 = (E - float(sc["z0_level"]) * (E - c["eb_sigma"] ** 2 - c["dh_T"])) - c["eb_sigma"] ** 2 - c["dh_T"] + c["dh_T"]
    out.update(theta_K=theta, mu_int=mu_int, c_int=c_int, beta0=beta0)
    T1 = d_drift / (d_h * r_W2) if r_W2 > 0 else math.inf
    T2 = alpha * (1 - math.exp(-theta * (beta0 - h0))) / (2 * c_int) if c_int > 0 else math.inf
    T_TTL = min(T1, T2, T_perp)
    out.update(T_TTL=T_TTL, T_TTL_drift=T1, T_TTL_kushner=T2)
    conds = {"margin": dts + d_stoch + d_drift <= dc["delta_cert"],
             "beta0_gt_h0": beta0 > h0,
             "h_act_minus_pos": h_act_minus > 0,
             "r_min2_pos": r2 > 0,
             "T_le_TTL": T <= T_TTL,
             "theorem": bool(cert["theorem_valid"])}
    out["conditions"] = conds
    out["valid"] = bool(all(conds.values()))
    out["bound"] = alpha
    return out


def clopper_pearson(k: int, n: int, conf: float):
    """Exact two-sided Clopper-Pearson interval for a binomial proportion."""
    a = (1 - conf) / 2
    lo = 0 if k == 0 else float(beta_dist.ppf(a, k, n - k + 1))
    hi = 1 if k == n else float(beta_dist.ppf(1 - a, k + 1, n - k))
    return lo, hi
