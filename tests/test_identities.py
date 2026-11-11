"""Identity tests: Proposition prop:dahlquist (eq:mu_pe), Lemma lem:pencil_rate, Lemma lem:quad_split,
(H2a), exp_weights partition of unity, phi hold factor (closed form vs quadrature)."""
import numpy as np
from scipy.integrate import quad

from pgcbf.certificates import exp_weights, mu_2P, pencil_rate, phi_hold


def test_prop2_mu_equals_minus_gamma_half(system):
    cfg, s, c = system
    lat = s["lat"]
    tol = float(cfg["tolerances"]["roundoff_rel"])
    assert abs(c["mu_2P_Ar0"] + lat.gamma / 2) <= tol * lat.gamma * lat.N
    assert c["Q_identity_relerr"] <= tol * lat.N
    assert abs(mu_2P(lat.Af, lat.Pfe) + lat.gamma / 2) <= tol * lat.gamma * lat.N * lat.omega2.max()


def test_tail_rate_equals_pencil_rate(system):
    cfg, s, c = system
    lat, rom = s["lat"], s["rom"]
    tol = float(cfg["tolerances"]["roundoff_rel"])
    pr = pencil_rate(lat.eps, lat.omega[rom.K_m:], lat.gamma).min()
    assert abs(c["lambda_K1"] - pr) <= tol * lat.N * lat.omega2.max()
    assert c["lambda_K1"] > 0


def test_quad_split_1e12(system):
    cfg, s, c = system
    lat, rom = s["lat"], s["rom"]
    tol = float(cfg["tolerances"]["roundoff_rel"])
    rng = np.random.default_rng(int(cfg["seeds"]["mc"]))
    n = int(cfg["stochastic"]["n_paths"])
    E = c["E_safe"]
    Z = rng.standard_normal((rom.K, n))
    Z = Z / rom.Pnorm(Z) * np.sqrt(E * rng.random(n))
    Xp = rom.Pperp @ rng.standard_normal((lat.N, n))
    Xp = Xp / rom.tail_Pnorm(Xp) * c["eps_basis_P"] * rng.random(n)
    lhs = lat.h_barrier(rom.Phi @ Z + Xp, E)
    rhs = rom.hPhi(Z, E) - np.einsum("ip,ij,jp->p", Xp, lat.Pfe, Xp)
    assert np.max(np.abs(lhs - rhs)) <= tol * np.abs(lhs).max()
    assert c["H2a_residual"] <= float(cfg["tolerances"]["h2a"]) * np.linalg.norm(lat.Pfe, 2)


def test_exp_weights_and_phi(system):
    cfg, s, c = system
    tol = float(cfg["tolerances"]["roundoff_rel"])
    g0, Ts, tA, tB = c["gamma0"], c["T_s"], c["tau_A"], c["tau_B"]
    w = exp_weights(g0, Ts, tA, tB)
    assert abs(sum(w) - 1) <= tol and min(w) > 0
    Z = (np.exp(g0 * Ts) - 1) / g0
    ref = [quad(lambda t: np.exp(g0 * t), a, b, epsabs=0, epsrel=tol)[0] / Z for a, b in ((0, tA), (tA, tB), (tB, Ts))]
    assert np.allclose(w, ref, rtol=np.sqrt(tol), atol=0)
    phi_ref = quad(lambda t: t * np.exp(g0 * t), 0, Ts, epsabs=0, epsrel=tol)[0] / Z / Ts
    assert abs(c["phi"] - phi_ref) <= np.sqrt(tol) * phi_ref
    assert abs(phi_hold(g0 * Ts, float(cfg["tolerances"]["phi_series_below"])) - c["phi"]) == 0
