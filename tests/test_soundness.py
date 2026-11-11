"""Soundness tests: Proposition prop:twin_crown (>= 1e5 samples, monotone in b_B), the two-branch
delta_cert vs brute force (class H and class D), and Theorem thm:main_cbf (a)-(c) on a short run."""
import copy
import itertools

import numpy as np

from pgcbf.certificates import build, bracket, certificate, classD_closed_form, sample_shell
from pgcbf.crown import TwinCrown
from pgcbf.simulate import Simulator, initial_conditions


def _gap(s, Z):
    pol = s["policy"]
    return np.linalg.norm(pol.psi_B(Z) - pol.teacher(Z), axis=0)


def test_twin_crown_sound_1e5(system):
    cfg, s, c = system
    n = int(cfg["experiments"]["crown_samples"])
    assert n >= 100000
    rng = np.random.default_rng(int(cfg["seeds"]["mc"]))
    D, K = s["D_ODD"], s["rom"].K
    Z = np.hstack([rng.uniform(-D / 2, D / 2, (K, n)), np.array(list(itertools.product([-D / 2, D / 2], repeat=K))).T])
    assert c["eps_CROWN"] >= _gap(s, Z).max()


def test_twin_crown_monotone_in_bB(cfg):
    eps = []
    for b in cfg["experiments"]["crown_b_B_sweep"]:
        cf = copy.deepcopy(cfg)
        cf["controller"]["b_B"] = int(b)
        cf["controller"]["rounding"] = "nearest"
        s = build(cf, train=False)
        eps.append(TwinCrown(s["qnet"], s["policy"].portmap).eps_crown(s["D_ODD"], s["rom"].K, int(cf["crown"]["splits_per_dim"]))[0])
    assert all(a >= b for a, b in zip(eps, eps[1:]))


def test_delta_cert_below_brute_force(system):
    cfg, s, c = system
    rng = np.random.default_rng(int(cfg["seeds"]["mc"]))
    n = int(cfg["experiments"]["crown_samples"])
    Z = sample_shell(s["rom"], c["r_min2"], c["h_act"], c["Delta_h_dip"], n, rng)
    br = bracket(s["rom"], s["policy"], Z, c["E_safe"], c["eps_basis_P"], c["Delta_h_dip"], c["gamma0"])
    assert c["delta_total"] <= c["delta_cert"] <= br.min()


def test_classD_two_branch_matches_affine_brute_force(cfg):
    d = cfg["experiments"]["classD"]
    r2, ha, dh = float(d["r_min2"]), float(d["h_act"]), float(d["dh"])
    rho = np.linspace(r2, r2 + ha + dh, int(cfg["experiments"]["shell_samples"]))
    for g0 in d["gamma0_values"]:
        for cM in (float(g0) / 2, float(g0) * 2):  # both branches
            cf = classD_closed_form(cM, r2, ha, dh, float(g0))
            affine = cM * rho + float(g0) * (r2 + ha - rho)  # c ||z||^2 + gamma0 h~, h~ = r2 + h_act - rho on C_act
            assert cf["branch"] == (1 if cM >= float(g0) else 2)  # Lemma lem:delta_cert_closed
            assert abs(cf["delta_cert"] - affine.min()) <= float(cfg["tolerances"]["roundoff_rel"]) * abs(affine).max()


def test_theorem1_short_run(system):
    cfg, s, c = system
    sc = cfg["sim"]
    X0 = initial_conditions(s["lat"], s["rom"], c, int(sc["n_seeds"]), np.random.default_rng(int(cfg["seeds"]["ic"])),
                            float(sc["ic_level_fraction"]), float(sc["ic_tail_fraction"]), int(sc["ic_n_tail_modes"]))
    sim = Simulator(s["lat"], s["rom"], s["policy"], cfg, adversary="worst")
    n = int(cfg["experiments"]["latch_periods"]) * int(cfg["experiments"]["latch_grid_steps"])
    r = sim.run(X0, n, record="substep")
    assert r["h"].min() >= 0
    assert r["tail"].max() <= c["eps_basis_P"] * (1 + float(cfg["tolerances"]["roundoff_rel"]))


def _Lfh(rom, P, Z):
    """L_{f_r} h~_r(z) = grad h~_r(z)^T f_r(z) = -2 z^T P f_r(z), columnwise (h~_r = E_safe - z^T P z)."""
    return -2 * np.einsum("ip,ij,jp->p", Z, P, rom.f_r(Z))


def test_Lfh_lipschitz_dominates_finite_differences(system):
    """Regression test for soundness defect R-L1 (missing factor 2 in grad of the quadratic form in L_{L_f h}).

    The certified constant Lip(L_f h) of eq:reduced_constants must dominate the empirical difference quotient
    |L_f h(z1) - L_f h(z2)| / ||z1 - z2||_2 over pairs in C_r0 = {z^T P z <= E_safe}: random pairs plus
    finite-difference pairs at the boundary maximiser of ||(P A_r + A_r^T P) z||_2 along the gradient direction.
    The quotient must also exceed the pre-fix (factor-2-less) value, so the test is not vacuous: it fails on the
    old code."""
    from pgcbf.certificates import psd_isqrt, opnorm
    cfg, s, c = system
    rom, E = s["rom"], c["E_safe"]
    P = rom.P
    rng = np.random.default_rng(int(cfg["seeds"]["mc"]))
    n = int(cfg["experiments"]["crown_samples"]) // 10
    G1, G2 = rng.standard_normal((rom.K, n)), rng.standard_normal((rom.K, n))
    Z1 = G1 / rom.Pnorm(G1) * np.sqrt(E * rng.random(n))
    Z2 = G2 / rom.Pnorm(G2) * np.sqrt(E * rng.random(n))
    q_rand = np.abs(_Lfh(rom, P, Z1) - _Lfh(rom, P, Z2)) / np.linalg.norm(Z1 - Z2, axis=0)
    # targeted pairs: z* = sqrt(E) P^{-1/2} v (v top right singular vector of sym P^{-1/2}), step along grad
    Ar = rom.Ar0 + rom.A_adv_r
    sym = P @ Ar + Ar.T @ P
    Pih = psd_isqrt(P)
    _, _, Vt = np.linalg.svd(sym @ Pih)
    q_fd = []
    for sgn in (1, -1):
        zs = sgn * np.sqrt(E) * (Pih @ Vt[0])
        d = sym @ zs
        d /= np.linalg.norm(d)
        for hstep in np.sqrt(E) * np.logspace(-6, -3, 4):
            z1 = zs * (1 - hstep / np.sqrt(E))  # strictly inside C_r0
            z2 = z1 + hstep * d
            if float(z2 @ P @ z2) > E:
                z2 = z1 - hstep * d
            assert float(z1 @ P @ z1) <= E and float(z2 @ P @ z2) <= E
            q_fd.append(abs(_Lfh(rom, P, z1[:, None])[0] - _Lfh(rom, P, z2[:, None])[0]) / hstep)
    q_max = max(q_rand.max(), max(q_fd))
    tol = float(cfg["tolerances"]["roundoff_rel"])
    assert q_max <= c["L_Lfh"] * (1 + tol), (q_max, c["L_Lfh"])
    L_Lfh_prefix = c["L_Lfh"] - opnorm(sym @ Pih) * np.sqrt(E)  # the value without the factor 2 (defect R-L1)
    assert q_max > L_Lfh_prefix, (q_max, L_Lfh_prefix)
