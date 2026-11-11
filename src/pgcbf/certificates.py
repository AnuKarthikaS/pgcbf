"""Every certificate quantity of Theorem thm:main_cbf, computed at run time in float64.

Each function names the equation label it implements.  No numerical constant other
than configuration inputs and formula integers appears in this module.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.linalg import eigh, sqrtm

from .lattice import Lattice
from .rom import ROM
from .controller import Teacher, QuantizedNet, Policy, learn_rounding
from .crown import TwinCrown, E_q


class CertificateRefused(AssertionError):
    """Raised when a hypothesis of Theorem thm:main_cbf (or a corollary) fails."""


# =============================================================== linear algebra
def psd_sqrt(P):
    w, V = np.linalg.eigh(P)
    return (V * np.sqrt(w)) @ V.T


def psd_isqrt(P):
    w, V = np.linalg.eigh(P)
    return (V / np.sqrt(w)) @ V.T


def opnorm(A):
    return float(np.linalg.norm(A, 2)) if A.size else 0


def mu_2P(A, P):
    """Implements the logarithmic norm of eq:mu_p0/eq:mu_pe: mu_{2,P}(A) = lambda_max(P^{1/2} A P^{-1/2} sym) (Prop. prop:dahlquist)."""
    S, Si = psd_sqrt(P), psd_isqrt(P)
    B = S @ A @ Si
    return float(np.linalg.eigvalsh((B + B.T) / 2).max())


def gen_eig_min(M, P):
    """lambda_min of the pencil (M, P): min z^T M z / z^T P z."""
    return float(eigh((M + M.T) / 2, P, eigvals_only=True).min())


# =============================================================== class-(H) checks
def eps_choice_check(lat: Lattice):
    """Implements eq:eps_choice: 0 < gamma < 2 omega_min, eps = gamma/2, P_full,eps > 0."""
    return {"gamma_lt_2omega_min": bool(0 < lat.gamma < 2 * lat.omega_min),
            "lambda_min_Pfe": float(np.linalg.eigvalsh(lat.Pfe).min()),
            "omega_min": float(lat.omega_min)}


def dahlquist(rom: ROM):
    """Implements eq:mu_pe (Proposition prop:dahlquist): Q_{gamma/2} = gamma P_{gamma/2}, mu_{2,P}(A_r0) = -gamma/2."""
    P, A = rom.P, rom.Ar0
    Q = -(P @ A + A.T @ P)
    g = rom.lat.gamma
    return {"mu_2P": mu_2P(A, P), "Q_minus_gammaP": opnorm(Q - g * P) / opnorm(g * P), "mu_R": -mu_2P(A, P)}


def tail_rate(rom: ROM):
    """Implements eq:tail_dissipation: largest lambda_{K+1} with -x^T(P A + A^T P)x >= 2 lambda ||x||_P^2 on V_perp."""
    lat = rom.lat
    B = rom.Phi_perp
    Pt = B.T @ lat.Pfe @ B
    Qt = -(B.T @ (lat.Pfe @ lat.Af + lat.Af.T @ lat.Pfe) @ B)
    return gen_eig_min(Qt, Pt) / 2


def pencil_rate(eps, omega, gamma):
    """Implements Lemma lem:pencil_rate / Remark rem:tail_rates: modal pencil rate r(eps, omega) = gamma/2 - omega|gamma-2eps|/(2 sqrt(omega^2-eps^2))."""
    return gamma / 2 - omega * abs(gamma - 2 * eps) / (2 * np.sqrt(omega**2 - eps**2))


# =============================================================== main certificate
def build(cfg: dict, train: bool = True):
    """Construct plant, ROM, Teacher, Student and policy from the configuration."""
    lat = Lattice(cfg)
    rom = ROM(lat, cfg["rom"]["K_m"])
    c = cfg["controller"]
    P = rom.P
    portmap = rom.g_r.T @ P
    lamminP = float(np.linalg.eigvalsh(P).min())
    E = float(cfg["safety"]["E_safe"])
    D_ODD = (1 + float(cfg["safety"]["D_ODD_margin"])) * 2 * math.sqrt(E / lamminP)
    teacher = Teacher(lat.m, int(c["H"]), int(c["L"]), float(c["theta"]), float(c["L_T"]), int(cfg["seeds"]["teacher"]))
    K = rom.K
    Y = opnorm(portmap) * math.sqrt(K) * D_ODD / 2
    offsets, hist = None, None
    if c["rounding"] == "learned" and train:
        rng = np.random.default_rng(int(cfg["seeds"]["train_data"]))
        Ztr = rng.uniform(-D_ODD / 2, D_ODD / 2, (K, int(c["train"]["samples"])))
        offsets, hist = learn_rounding(teacher, int(c["b_B"]), Y, portmap @ Ztr, c["train"], int(cfg["seeds"]["train"]))
    qnet = QuantizedNet(teacher, int(c["b_B"]), Y, offsets)
    D_v = port_box(rom, E)
    pol = Policy(teacher, qnet, portmap, float(c["u_max"]), D_ODD, int(c["b_A"]), deadzone=bool(c["deadzone"]),
                 lut_space=c["lut_space"], D_v=D_v)
    return {"lat": lat, "rom": rom, "teacher": teacher, "qnet": qnet, "policy": pol, "D_ODD": D_ODD, "Y": Y, "D_v": D_v,
            "train_hist": hist, "offsets": offsets}


def port_box(rom: ROM, E: float) -> float:
    """Implements the S3 port-space Path A (Section sec:problem, LUT Path; Remark rem:lut_variants):
    D_v := 2 sup_{z in C_r0} ||g_r^T P z||_inf, closed form on C_r0 = {||z||_P^2 <= E_safe}:
    D_v = 2 sqrt(E_safe) max_j sqrt((g_r^T P g_r)_jj)  (Cauchy-Schwarz in the P-inner product, attained)."""
    return 2 * math.sqrt(E) * math.sqrt(float(np.max(np.diag(rom.g_r.T @ rom.P @ rom.g_r))))


def exp_weights(g0, Ts, tA, tB):
    """Implements eq:exp_weights."""
    den = math.expm1(g0 * Ts)
    w1 = math.expm1(g0 * tA) / den
    w2 = (math.exp(g0 * tB) - math.exp(g0 * tA)) / den
    w3 = (math.exp(g0 * Ts) - math.exp(g0 * tB)) / den
    return w1, w2, w3


def phi_hold(xi, series_below):
    """Implements the hold factor phi of eq:delta_total_exact (integration-by-parts identity after eq:exact_exp_integral): phi(xi) = ((xi-1)e^xi + 1)/(xi(e^xi - 1)) (footnote series near 0)."""
    if abs(xi) < series_below:
        return 1 / 2 + xi / (2 * math.factorial(3)) - xi**3 / math.factorial(2 * 3)  # exact Taylor coefficients
    return ((xi - 1) * math.exp(xi) + 1) / (xi * math.expm1(xi))


# =============================================================== per-quantity functions
def odd_containment(lat: Lattice, rom: ROM, E, odd_margin, D_ODD):
    """Implements eq:H0 (and the second half of (H0)): rho_ODD and E_safe < lambda_min(Pfe) rho^2, E < D^2 lambda_min(P)/4."""
    lamminPfe = float(np.linalg.eigvalsh(lat.Pfe).min())
    lamminP = float(np.linalg.eigvalsh(rom.P).min())
    rho = (1 + odd_margin) * math.sqrt(E / lamminPfe)
    return {"rho_ODD": rho, "D_ODD": D_ODD, "H0": bool(E < lamminPfe * rho**2 and E < D_ODD**2 * lamminP / 4),
            "lambda_min_P": lamminP}


def sin_lipschitz(radius):
    """Lipschitz modulus of q - sin q on |q| <= radius: 1 - cos(radius) (2 beyond pi); used in assump:generic."""
    return (1 - math.cos(radius)) if radius < math.pi else 2


def lipschitz_Fnl(lat: Lattice, rho_odd):
    """Implements Assumption assump:generic: L_Fnl = ||A_adv|| + kappa (1 - cos rho_ODD) on X_ODD."""
    return opnorm(lat.A_adv) + lat.kappa * sin_lipschitz(rho_odd)


def momentum_injector(lat: Lattice):
    return np.vstack([np.zeros((lat.Nm, lat.Nm)), np.eye(lat.Nm)])


def M_perp_bound(lat: Lattice, rom: ROM, rho_odd, u_max, eps_adv):
    """Implements eq:M_perp: sound triangle bound of sup ||Pfe^{1/2} P_perp (F + G u + D w)|| on X_ODD."""
    PP = psd_sqrt(lat.Pfe) @ rom.Pperp
    Ep = momentum_injector(lat)
    lin = opnorm(PP @ (lat.Af @ rom.Phi @ rom.Zmap + lat.A_adv)) * rho_odd
    nl = opnorm(PP @ Ep) * lat.kappa * rho_odd**3 / math.factorial(3)  # ||q - sin q||_2 <= ||q||_2^3 / 3!
    u = opnorm(PP @ lat.G) * u_max
    w = opnorm(PP @ lat.D) * eps_adv
    return {"M_perp": lin + nl + u + w, "M_perp_lin": lin, "M_perp_nl": nl, "M_perp_u": u, "M_perp_w": w,
            "PperpD_norm": opnorm(PP @ lat.D)}


def eps_basis_P(M_perp, lam):
    """Implements eq:eps_basis_P: eps_basis_P = M_perp / lambda_{K+1}."""
    return M_perp / lam


def L_perp_P(lat: Lattice, rom: ROM, rho_odd):
    """Implements eq:Lperp_def: sharp linear part + sin-nonlinearity part of L_perp_P (Pfe-norm on V_perp).

    Delta_perp(z, x_perp) = A_{r,perp} x_perp + P0^{-1} Phi^T P_full [F_nl(Phi z + x_perp) - F_nl(Phi z)] with
    A_{r,perp} = 0 in class (H) (eigenmode invariance) and F_nl = A_adv x + F_sin(x); the linear part is therefore
    the advective block A_{adv,r,perp} = Zmap A_adv P_perp, the remainder is the Lipschitz bound of F_sin."""
    Bp = rom.Phi_perp
    Wt = psd_isqrt(Bp.T @ lat.Pfe @ Bp)  # x_perp = Bp y, ||x_perp||_Pfe = ||Wt^{-1} y||
    lin = opnorm(rom.Zmap @ (lat.Af + lat.A_adv) @ Bp @ Wt)
    qsel = np.hstack([np.eye(lat.Nm), np.zeros((lat.Nm, lat.Nm))])
    nl = opnorm(rom.Zmap @ momentum_injector(lat)) * lat.kappa * sin_lipschitz(rho_odd) * opnorm(qsel @ Bp @ Wt)
    return {"L_perp": lin + nl, "L_perp_lin": lin, "L_perp_nl": nl,
            "A_adv_rperp_norm": opnorm(rom.A_adv_rperp), "A_rperp_full_norm": opnorm(rom.Arperp)}


def reduced_constants(lat: Lattice, rom: ROM, E, u_max, eps_adv, L_perp, eps_b):
    """Implements eq:reduced_constants: L_h, g_h, d_h, M_z, L_{L_f h}, L_{L_g h} over C_r0 = {||z||_P^2 <= E}."""
    P = rom.P
    Ph, Pih = psd_sqrt(P), psd_isqrt(P)
    lammaxP = float(np.linalg.eigvalsh(P).max())
    L_h = 2 * math.sqrt(lammaxP * E)
    g_h = 2 * math.sqrt(E) * opnorm(Ph @ rom.g_r)
    d_h = 2 * math.sqrt(E) * opnorm(Ph @ rom.d_r)
    Ar = rom.Ar0 + rom.A_adv_r
    qr = np.hstack([np.eye(rom.K_m), np.zeros((rom.K_m, rom.K_m))])
    s_q = opnorm(qr @ Pih) * math.sqrt(E)  # sup ||z_q||_2 on C_r0
    Zp = rom.Zmap @ momentum_injector(lat)
    n_r_sup = opnorm(Zp) * lat.kappa * s_q**3 / math.factorial(3)
    M_z = (opnorm(Ar @ Pih) * math.sqrt(E) + n_r_sup + opnorm(rom.g_r) * u_max
           + opnorm(rom.d_r) * eps_adv + L_perp * eps_b)
    sym = P @ Ar + Ar.T @ P
    # grad_z(-z^T sym z) = -2 sym z, so its sup over C_r0 is 2 ||sym P^{-1/2}|| sqrt(E) (factor 2 restored in the
    # review round; the previous value under-estimated Lip(L_f h) by about 2x: soundness defect R-L1)
    L_Lfh = (2 * opnorm(sym @ Pih) * math.sqrt(E) + 2 * opnorm(P) * n_r_sup
             + 2 * opnorm(Zp) * lat.kappa * sin_lipschitz(s_q) * math.sqrt(lammaxP * E))
    L_Lgh = 2 * opnorm(rom.g_r.T @ P)
    return {"L_h": L_h, "g_h": g_h, "d_h": d_h, "M_z": M_z, "L_Lfh": L_Lfh, "L_Lgh": L_Lgh, "s_q": s_q}


def teacher_constants(rom: ROM, teacher: Teacher, u_max, D_ODD):
    """Implements eq:LTz (L_Tz = L_T ||g_r^T P||), Lemma lem:sector (eta_T) and eq:sat_cond ((H-sat))."""
    K, m = rom.K, rom.lat.m
    L_Tz = teacher.L_T * opnorm(rom.g_r.T @ rom.P)
    Gport = rom.g_r.T @ rom.P @ rom.P @ rom.g_r
    lhs = L_Tz * math.sqrt(K) * D_ODD / 2
    gmin = float(np.linalg.eigvalsh(Gport).min())
    return {"L_T": teacher.L_T, "L_Tz": L_Tz, "eta_T": teacher.eta, "u_max": u_max, "H_sat_lhs": lhs,
            "G_port_min_eig": gmin, "H_sat": bool(lhs <= u_max), "H3": bool(lhs <= u_max and gmin > 0 and teacher.H // 2 >= m)}


def eps_A_triangle(eC, eCv, L_lut, K_lut, D_lut, b_A):
    """Implements eq:eps_A_triangle: eps_A = eps_CROWN + L sqrt(K) D / 2^{b_A+1} (LUT index rounding)."""
    return max(eC, eCv) + L_lut * math.sqrt(K_lut) * D_lut / 2 ** (b_A + 1)


def codesign_window(K_lut, D_lut, M_lut, Ts, b_A, B_max, S_log2, b_out, m):
    """Implements eq:codesign_window and eq:lut_capacity in the port-space form (S3: m (b_A+1) <= B_max,
    m (b_A+1) + ceil(log2(m b_out)) <= log2 S_bits, lower bound floor(log2(sqrt(m) D_v/(2 M_v T_s))) + 1,
    upper bound floor(B_max/m) - 1); the state-space variant passes K, D_ODD, M_z instead."""
    lo = math.floor(math.log2(math.sqrt(K_lut) * D_lut / (2 * M_lut * Ts))) + 1
    hi = B_max // K_lut - 1
    addr = K_lut * (b_A + 1)
    store = addr + math.ceil(math.log2(m * b_out))
    H5 = bool(b_A >= 1 and addr <= B_max and store <= S_log2)  # (H5): b_A >= 1 satisfies eq:lut_capacity
    return {"codesign_lo": lo, "codesign_hi": hi, "codesign_window": bool(lo <= b_A <= hi), "lut_address_bits": addr,
            "lut_storage_bits": store, "H5": H5}


def delta_dip(g_h, eA, eC, du_lag, tA, tB, L_hold, Ts):
    """Implements eq:delta_dip: Delta h_dip = g_h (eps_A - eps_C)(tau_B - tau_A) + g_h du_lag tau_A + L_hold T_s^2 / 2."""
    return g_h * (eA - eC) * (tB - tA) + g_h * du_lag * tA + L_hold * Ts**2 / 2


def delta_total_exact(L_h, L_perp, eps_b, g_h, w1, w2, eA, eC, du_lag, d_h, eps_adv, phi, L_hold, Ts, drop_quant=False):
    """Implements eq:delta_total_exact (term by term)."""
    terms = {"tail": L_h * L_perp * eps_b,
             "quant": 0 * g_h if drop_quant else g_h * (w2 * eA + (1 - w2) * eC),
             "lag": w1 * g_h * du_lag,
             "adv": d_h * eps_adv,
             "hold": phi * L_hold * Ts}
    return sum(terms.values()), terms


def shell(E, eps_b, dh, L_h, M_z, Ts, h_act_fraction):
    """Implements eq:eps_basis_def / eq:shell / (H4): h_act in (L_h M_z T_s, E - eps_b^2 - dh), r_min^2 = room - h_act."""
    LMT = L_h * M_z * Ts
    room = E - eps_b**2 - dh
    h_act = LMT + h_act_fraction * (room - LMT)
    r2 = room - h_act
    return {"H4_slack": room - LMT, "h_act": h_act, "r_min2": r2, "r_min": math.sqrt(r2) if r2 > 0 else math.nan,
            "H4": bool(room > LMT and r2 > 0), "C_tilde_room": room}


def sector_modulus(lat: Lattice, rom: ROM, s_q):
    """L_sec := sup ||P^{1/2} grad F_nl,r P^{-1/2}|| on C_r0 (Prop. prop:deadzone_stability)."""
    Ph, Pih = psd_sqrt(rom.P), psd_isqrt(rom.P)
    qr = np.hstack([np.eye(rom.K_m), np.zeros((rom.K_m, rom.K_m))])
    Zp = rom.Zmap @ momentum_injector(lat)
    return opnorm(Ph @ rom.A_adv_r @ Pih) + opnorm(Ph @ Zp) * lat.kappa * sin_lipschitz(s_q) * opnorm(qr @ Pih)


def rss_ball(rom: ROM, eA, eC, du_lag, mu_R, L_sec, r_min):
    """Implements eq:rss_ball and eq:deadzone_cond ((H-dz) r_ss < r_min)."""
    eq_bar = max(eA, eC + du_lag) + du_lag
    lammaxP = float(np.linalg.eigvalsh(rom.P).max())
    r_ss = math.sqrt(lammaxP) * opnorm(rom.g_r) * eq_bar / (mu_R - L_sec) if mu_R > L_sec else math.inf
    return {"eps_q_bar": eq_bar, "L_sec": L_sec, "r_ss": r_ss, "Hdz": bool(mu_R > L_sec and r_ss < r_min)}


def certificate(cfg: dict, sysd: dict | None = None, train: bool = True, ablate: tuple = ()) -> dict:
    """Evaluate every constant of Section sec:methodology / Theorem thm:main_cbf (Implements eq:delta_total_exact,
    eq:control_authority and the hypothesis ledger (H0)-(H7) of Theorem thm:main_cbf, plus (H-dz) of
    Proposition prop:deadzone_stability under the key Hdz).  theorem_valid = all of (H0)-(H7) and the timing order
    and the codesign window (the latter is optional in the paper; this implementation requires it).

    ablate: names of terms deliberately removed to build the *naive* (unsound) certificates of
    the ablation study: "tail" (eps_basis_P := 0, L_perp := 0) and "quant" (w_2 eps_A + (1-w_2) eps_CROWN := 0)."""
    if sysd is None:
        sysd = build(cfg, train=train)
    lat, rom, teacher, qnet, pol = sysd["lat"], sysd["rom"], sysd["teacher"], sysd["qnet"], sysd["policy"]
    sf, ct, tm, hw = cfg["safety"], cfg["controller"], cfg["timing"], cfg["hardware"]
    E, u_max, eps_adv, g0 = float(sf["E_safe"]), float(ct["u_max"]), float(sf["eps_adv"]), float(sf["gamma0"])
    Ts, tA, tB = float(tm["T_s"]), float(tm["tau_A"]), float(tm["tau_B"])
    K, m = rom.K, lat.m
    out = {"N": lat.N, "K": K, "m": m, "Pe": lat.peclet, "E_safe": E, "eps_adv": eps_adv}
    out.update(eps_choice_check(lat))
    out.update(odd_containment(lat, rom, E, float(sf["odd_margin"]), sysd["D_ODD"]))
    rho = out["rho_ODD"]
    D_ODD = sysd["D_ODD"]
    out["L_Fnl"] = lipschitz_Fnl(lat, rho)
    out["K_f"] = lat.n_f
    out["adv_cutoff_weight"] = float(lat.adv_w[-1])
    # (H2)
    dq = dahlquist(rom)
    out.update(mu_2P_Ar0=dq["mu_2P"], Q_identity_relerr=dq["Q_minus_gammaP"], mu_R=dq["mu_R"])
    out["H2a_residual"] = rom.h2a_residual()
    out["H2a_invariance_residual"] = rom.h2a_invariance_residual()
    out["H2b_margin"] = rom.h2b_margin()
    lam = tail_rate(rom)
    out["lambda_K1"] = lam
    out["H2"] = bool(out["H2a_residual"] <= float(cfg["tolerances"]["h2a"]) * opnorm(lat.Pfe)
                     and out["H2b_margin"] > 0 and lam > 0 and out["gamma_lt_2omega_min"])
    # tail
    out.update(M_perp_bound(lat, rom, rho, u_max, eps_adv))
    eps_b = 0 * lam if "tail" in ablate else eps_basis_P(out["M_perp"], lam)
    out["eps_basis_P"] = eps_b
    out.update(L_perp_P(lat, rom, rho))
    if "tail" in ablate:
        out["L_perp"] = 0 * out["L_perp"]
    L_perp = out["L_perp"]
    # reduced constants, teacher
    out.update(reduced_constants(lat, rom, E, u_max, eps_adv, L_perp, eps_b))
    out.update(teacher_constants(rom, teacher, u_max, D_ODD))
    # Prop. prop:twin_crown
    tc = TwinCrown(qnet, pol.portmap)
    eC_crown, _ = tc.eps_crown(D_ODD, K, int(cfg["crown"]["splits_per_dim"]))
    Eq = E_q(qnet, sysd["Y"])
    eC = min(eC_crown, Eq)
    if pol.lut_space == "port":  # Path-A entries are Psi_B(v_q), v_q in the port box: certify that box too
        eCv, _ = TwinCrown(qnet, np.eye(m)).eps_crown(pol.D_lut, m, int(cfg["crown"]["splits_per_dim"]))
        out["eps_CROWN_portbox"] = eCv
    else:
        eCv = eC
    out.update(eps_CROWN_backward=eC_crown, E_q=Eq, eps_CROWN=eC, b_B=qnet.b_B)
    # LUT
    b_A = pol.b_A
    M_z, L_Tz, g_h = out["M_z"], out["L_Tz"], out["g_h"]
    du_lag = L_Tz * M_z * Ts
    if pol.lut_space == "state":
        K_lut, D_lut, L_lut, M_lut = K, D_ODD, L_Tz, M_z
    else:  # Remark rem:lut_variants(i): K -> m, D_ODD -> D_v, M_z -> M_v, L_Tz -> L_T
        K_lut, D_lut, L_lut, M_lut = m, pol.D_lut, teacher.L_T, opnorm(rom.g_r.T @ rom.P) * M_z
    eA = eps_A_triangle(eC, eCv, L_lut, K_lut, D_lut, b_A)
    out.update(eps_A=eA, du_lag=du_lag, b_A=b_A, D_lut=D_lut, K_lut=K_lut)
    out.update(codesign_window(K_lut, D_lut, M_lut, Ts, b_A, int(hw["B_max"]), int(hw["S_bits_log2"]), int(hw["b_out"]), m))
    out["M_v"] = M_lut
    # margins
    w1, w2, w3 = exp_weights(g0, Ts, tA, tB)
    phi = phi_hold(g0 * Ts, float(cfg["tolerances"]["phi_series_below"]))
    L_hold = (out["L_Lfh"] + out["L_Lgh"] * u_max + g0 * out["L_h"]) * M_z
    dh = delta_dip(g_h, eA, eC, du_lag, tA, tB, L_hold, Ts)
    dtot, terms = delta_total_exact(out["L_h"], L_perp, eps_b, g_h, w1, w2, eA, eC, du_lag, out["d_h"], eps_adv, phi, L_hold, Ts,
                                    drop_quant="quant" in ablate)
    out.update(w1=w1, w2=w2, w3=w3, phi=phi, L_hold=L_hold, Delta_h_dip=dh, delta_total=dtot,
               **{"dt_" + k: v for k, v in terms.items()}, gamma0=g0, T_s=Ts, tau_A=tA, tau_B=tB,
               timing_ok=bool(0 < tA < tB < Ts) or bool(tA == tB == Ts))
    out.update(shell(E, eps_b, dh, out["L_h"], M_z, Ts, float(sf["h_act_fraction"])))
    cm = delta_cert_bound_classH(rom, teacher, E, lat.kappa, out["s_q"])
    out.update(lambda_min_MP=cm["lambda_min_MP"], s_nl=cm["s_nl"], c_M=cm["c_M"])
    dc = two_branch(cm["c_M"], out["r_min2"], out["h_act"], dh, g0)
    out.update(delta_cert=dc["delta_cert"], delta_cert_branch=dc["branch"])
    out["H6"] = bool(dtot <= out["delta_cert"])  # (H6) eq:control_authority
    out["H6_margin"] = out["delta_cert"] - dtot
    out["tail_share"] = terms["tail"] / dtot  # tail-coupling share L_h L_perp eps_basis_P / delta_total
    out["null_space_rate"] = null_space_rate(rom)
    out.update(rss_ball(rom, eA, eC, du_lag, dq["mu_R"], sector_modulus(lat, rom, out["s_q"]), out["r_min"]))
    out["tail_init_bound"] = out["M_perp"] / lam  # (H7): ||x_perp(0)||_Pfe <= M_perp / lambda_{K+1}
    out["H7"] = bool(out["C_tilde_room"] > 0 and math.isfinite(out["tail_init_bound"]))  # (H7) initial set non-empty
    out["H1"] = bool(all(math.isfinite(out[k]) for k in ("L_Fnl", "M_perp", "L_perp")))  # (H1) Assumption 1 over X_ODD
    out["theorem_valid"] = bool(all(out[h] for h in THEOREM_HYPS) and out["timing_ok"] and out["codesign_window"])
    out["prop3_valid"] = bool(out["Hdz"] and out["H3"])
    return out


# Theorem thm:main_cbf hypotheses (H0)-(H7); (H-dz) belongs to Proposition prop:deadzone_stability only.
THEOREM_HYPS = ("H0", "H1", "H2", "H3", "H4", "H5", "H6", "H7")


def two_branch(c, r2, h_act, dh, g0):
    """Implements eq:delta_cert_closed (two-branch form, applied to the sound class-(H) slope c_M). Two-branch infimum of an affine-in-rho lower bound (eq:delta_cert_closed form)."""
    b1 = c * r2 + g0 * h_act
    b2 = c * (r2 + h_act + dh) - g0 * dh
    return {"delta_cert": min(b1, b2), "branch": 1 if b1 <= b2 else 2, "b1": b1, "b2": b2}


def delta_cert_bound_classH(rom: ROM, teacher: Teacher, E, kappa, s_q):
    """Implements Lemma lem:delta_cert_classH (class-(H) certified margin, m <= K), the sound lower bound of
    eq:teacher_cbf, with Lemma lem:sector and eq:mu_pe:

    L_f h + L_g h T >= z^T M z - 2|z^T P n_r(z)| >= (lambda_min(M,P) - 2 s_nl) ||z||_P^2 = c_M ||z||_P^2,
    M := gamma P - (P A_adv,r + A_adv,r^T P) + 2 eta_T P g_r g_r^T P,  ||P^{1/2} n_r(z)|| <= s_nl ||z||_P
    (Q = gamma P by eq:mu_pe; sector Lemma lem:sector under (H-sat)); delta_cert by the two-branch infimum."""
    P = rom.P
    g = rom.lat.gamma
    Q = -(P @ rom.Ar0 + rom.Ar0.T @ P)
    S = P @ rom.A_adv_r + rom.A_adv_r.T @ P
    M = Q - S + 2 * teacher.eta * P @ rom.g_r @ rom.g_r.T @ P
    lam_MP = gen_eig_min(M, P)
    Ep = np.vstack([np.zeros((rom.lat.Nm, rom.lat.Nm)), np.eye(rom.lat.Nm)])
    Pih = psd_isqrt(P)
    qr = np.hstack([np.eye(rom.K_m), np.zeros((rom.K_m, rom.K_m))])
    # ||P^{1/2} n_r(z)|| <= ||P^{1/2} Zmap Ep|| kappa ||z_q||^3/6 <= ... * s_q^2 ||qr P^{-1/2}|| ||z||_P / 6
    s_nl = opnorm(psd_sqrt(P) @ rom.Zmap @ Ep) * kappa * s_q**2 * opnorm(qr @ Pih) / math.factorial(3)
    return {"lambda_min_MP": lam_MP, "s_nl": s_nl, "c_M": lam_MP - 2 * s_nl, "gamma": g, "M": M}


def null_space_rate(rom: ROM):
    """Implements eq:null_space_cond: min of z^T Q z / z^T P z on {g_r^T P z = 0} (= gamma in class (H))."""
    P = rom.P
    Q = -(P @ rom.Ar0 + rom.Ar0.T @ P)
    Nmat = null_basis(rom.g_r.T @ P)
    if Nmat.shape[1] == 0:
        return math.inf
    return gen_eig_min(Nmat.T @ Q @ Nmat, Nmat.T @ P @ Nmat)


def null_basis(A):
    u, s, vt = np.linalg.svd(A)
    rank = int(np.sum(s > s.max() * A.shape[1] * np.finfo(float).eps))
    return vt[rank:].T


def bracket(rom: ROM, pol: Policy, Z, E, eps_b, dh, g0):
    """Exact bracket L_f h~ + L_g h~ T(z) + gamma0 h~(z) of eq:teacher_cbf at columns of Z."""
    P = rom.P
    PZ = P @ Z
    fr = rom.f_r(Z)
    T = pol.teacher(Z)
    htil = rom.hPhi(Z, E) - eps_b**2 - dh
    return -2 * np.sum(PZ * fr, axis=0) - 2 * np.sum(PZ * (rom.g_r @ T), axis=0) + g0 * htil


def sample_shell(rom: ROM, r2, h_act, dh, n, rng):
    """Uniform-in-level samples of C_act = {r_min^2 <= ||z||_P^2 <= r_min^2 + h_act + dh}."""
    K = rom.K
    G = rng.standard_normal((K, n))
    Pih = psd_isqrt(rom.P)
    U = Pih @ G
    U = U / rom.Pnorm(U)
    lev = r2 + (h_act + dh) * rng.random(n)
    # include both endpoints explicitly
    lev[: n // 4] = r2
    lev[n // 4: n // 2] = r2 + h_act + dh
    return U * np.sqrt(lev)


def refuse_if_invalid(cert: dict, required=None, what="Theorem thm:main_cbf"):
    """The certificate refuses (raises) if any hypothesis (H0)-(H7), the timing order or the codesign window fails."""
    required = THEOREM_HYPS if required is None else required
    bad = [h for h in required if not cert.get(h, False)]
    if not cert.get("timing_ok", True):
        bad.append("timing")
    if not cert.get("codesign_window", True):
        bad.append("codesign window")
    if bad:
        raise CertificateRefused(f"{what}: hypotheses {bad} fail")
    return True


# =============================================================== class (D) closed-form lemma
def classD_closed_form(c, r2, h_act, dh, g0):
    """Implements eq:delta_cert_closed (Lemma lem:delta_cert_closed): min of the two branches, branch 1 iff c >= gamma0."""
    return two_branch(c, r2, h_act, dh, g0)
