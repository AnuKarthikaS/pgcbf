"""Experiment 8 - numerical integrity.

(i)  Stiffness: spectral radius of the lattice operator vs. the explicit-Euler/RK4 stability
     limits; the ETD2RK exponential integrator propagates the stiff linear part exactly.
(ii) Convergence: step halving on the nonlinear open-loop lattice; error ratio -> 4 (order 2).
(iii) Exact latch timing: with a linear plant (kappa = 0) and piecewise-constant inputs that switch
     at tau_A, tau_B inside every period, the simulator matches the closed-form augmented-matrix
     exponential solution to round-off, while switching at the nearest integrator grid point
     (step-resolution latching) produces an O(step) error.
(iv) Modal sparse propagation == dense nodal ETD2RK (validates the change of coordinates).
(v)  Bit-for-bit determinism: two identical closed-loop runs produce identical SHA-256 hashes."""
import math
import matplotlib.ticker
import numpy as np
from scipy.linalg import expm
from _common import Ledger, load_config, override, write_table, FIG, array_hash
from pgcbf.certificates import build, certificate
from pgcbf.simulate import Simulator, initial_conditions, phi_matrices
from pgcbf import plotting as pl

L = Ledger("exp8_numerics")
cfg = load_config("default")
ex, sc = cfg["experiments"], cfg["sim"]
sysd = build(cfg, train=True)
c = certificate(cfg, sysd)
lat, rom, pol = sysd["lat"], sysd["rom"], sysd["policy"]
X0 = initial_conditions(lat, rom, c, int(sc["ic_n_tail_modes"]), np.random.default_rng(int(cfg["seeds"]["ic"])),
                        float(sc["ic_level_fraction"]), float(sc["ic_tail_fraction"]), int(sc["ic_n_tail_modes"]))
# ---------------- (i) stiffness
ev = np.linalg.eigvals(lat.A_lin)
rad = np.abs(ev).max()
Ts = float(cfg["timing"]["T_s"])
Tc0 = float(ex["conv_T_s"])
amp_euler = np.abs(1 + Tc0 * ev).max()
L.check("stiffness: at the convergence-test step explicit Euler is unstable (max|1 + h lambda| > 1) "
        "while the exponential integrator is used", amp_euler > 1, spectral_radius=rad, stiffness_ratio=rad / np.abs(ev).min(),
        step=Tc0, explicit_euler_amplification=amp_euler, rad_times_Ts=rad * Ts)
# ---------------- (ii) convergence (open loop, nonlinear, single-rate segments)
Tc, Hc = float(ex["conv_T_s"]), float(ex["conv_horizon"])
cfc = override(cfg, timing={"T_s": Tc, "tau_A": Tc, "tau_B": Tc})
nper = int(round(Hc / Tc))
amp = np.sqrt(c["E_safe"] / c["lambda_min_Pfe"])  # large-amplitude state (nonlinearity active), ||x||_2 = sqrt(E/lambda_min)
Xc = X0 / np.linalg.norm(X0, axis=0) * amp


def run_sub(nsub, scheme="etd2"):
    sim = Simulator(lat, rom, pol, cfc, substeps=[nsub, nsub, nsub], adversary=None, controller=False, scheme=scheme)
    return sim.run(Xc, nper, record="period")["X_final"]


ref = run_sub(int(ex["conv_ref_substeps"]))
subs = [int(s) for s in ex["conv_substeps"]]
errs = [np.linalg.norm(run_sub(s) - ref) / np.linalg.norm(ref) for s in subs]
errs1 = [np.linalg.norm(run_sub(s, "etd1") - ref) / np.linalg.norm(ref) for s in subs]
ratios = [a / b for a, b in zip(errs, errs[1:])]
order = np.log2(ratios)
L.check("ETD2RK convergence: observed order in [1.8, 2.2] on step halving (error ratio ~ 4)",
        np.all(np.abs(order - 2) <= 1 / 5), ratios=ratios, orders=order)
L.check("ETD1 (exponential Euler) for reference: observed order ~ 1", np.all(np.abs(np.log2(np.array(errs1[:-1]) / np.array(errs1[1:])) - 1) <= 1 / 5),
        orders1=np.log2(np.array(errs1[:-1]) / np.array(errs1[1:])))
write_table("exp8_convergence", ["substeps per period", "step", "rel. error ETD2RK", "ratio", "rel. error ETD1"],
            [[s, Tc / s, e, (ratios[i - 1] if i > 0 else None), e1] for i, (s, e, e1) in enumerate(zip(subs, errs, errs1))],
            caption="Step-halving convergence of the exponential integrator (open-loop nonlinear lattice).")
# ---------------- (iii) exact latch timing
cfl = override(cfg, plant={"kappa": 0})
s_l = build(cfl, train=True)
lat_l = s_l["lat"]
rngu = np.random.default_rng(int(cfg["seeds"]["mc"]))
uA_c, uB_c = rngu.uniform(-1, 1, (lat_l.m, 1)), rngu.uniform(-1, 1, (lat_l.m, 1))
uA_c, uB_c = uA_c * pol.u_max / np.linalg.norm(uA_c), uB_c * pol.u_max / np.linalg.norm(uB_c) / 2


class ConstPolicy:
    """u_A, u_B constant but different: every period switches u_B -> u_A at tau_A and u_A -> u_B at tau_B."""
    D_ODD = np.inf

    def psi_A(self, Z):
        return np.repeat(uA_c, Z.shape[1], axis=1)

    def psi_B(self, Z):
        return np.repeat(uB_c, Z.shape[1], axis=1)


npl = int(ex["latch_periods"])
tA, tB = float(cfg["timing"]["tau_A"]), float(cfg["timing"]["tau_B"])
sim = Simulator(lat_l, s_l["rom"], ConstPolicy(), cfl, adversary=None)
Xs = sim.run(X0[:, :1], npl, record="period")["X_final"][:, 0]


def aug_step(x, u, dt):
    n = lat_l.N
    M = np.zeros((n + 1, n + 1))
    M[:n, :n] = lat_l.A_lin * dt
    M[:n, n] = (lat_l.G @ u)[:, 0] * dt
    E = expm(M)
    return E[:n, :n] @ x + E[:n, n]


x = X0[:, 0].copy()
for k in range(npl):
    x = aug_step(x, uB_c, tA); x = aug_step(x, uA_c, tB - tA); x = aug_step(x, uB_c, Ts - tB)
err_exact = np.linalg.norm(Xs - x) / np.linalg.norm(x)
# step-resolution latching: uniform integrator grid of n_grid steps per period, switches rounded to the grid
n_grid = int(ex["latch_grid_steps"])  # a uniform grid not commensurate with tau_A, tau_B
hg = Ts / n_grid
iA, iB = int(round(tA / hg)), int(round(tB / hg))
xg = X0[:, 0].copy()
for k in range(npl):
    xg = aug_step(xg, uB_c, iA * hg); xg = aug_step(xg, uA_c, (iB - iA) * hg); xg = aug_step(xg, uB_c, Ts - iB * hg)
err_grid = np.linalg.norm(xg - x) / np.linalg.norm(x)
tol = float(cfg["tolerances"]["roundoff_rel"])
L.check("exact latch timing: simulator == closed-form augmented exponential (kappa=0) to round-off", err_exact <= tol * npl * lat.N,
        rel_err=err_exact, periods=npl)
L.check("step-resolution latching would be O(step) wrong (discriminating test)", err_grid > err_exact * npl * lat.N,
        rel_err_grid=err_grid, grid_step=hg, tau_A=tA, tau_B=tB)
# ---------------- (iv) modal sparse vs dense nodal ETD2RK (closed loop, worst-case adversary)
nmod = int(ex["latch_periods"])
simm = Simulator(lat, rom, pol, cfg, adversary="worst")
Xm = simm.run(X0, nmod, record="period")["X_final"]
A = lat.A_lin
segs = [tA, tB - tA, Ts - tB]
mats = [phi_matrices(A, h) for h in segs]
Pd = rom.P @ rom.d_r


def Nn(X, U):
    v = Pd.T @ (rom.Zmap @ X)
    w = simm.eps_adv * v / np.linalg.norm(v, axis=0)
    return lat.nonlinear(X) + lat.G @ U + lat.D @ w


Xn = X0.copy()
uprev = pol.psi_B(rom.Zmap @ Xn)
for k in range(nmod):
    Z = rom.Zmap @ Xn
    ins = (uprev, pol.psi_A(Z), pol.psi_B(Z))
    for (E, p1, p2), h, U in zip(mats, segs, ins):
        N0 = Nn(Xn, U)
        Aa = E @ Xn + p1 @ (h * N0)
        Xn = Aa + p2 @ (h * (Nn(Aa, U) - N0))
    uprev = ins[2]
err_modal = np.linalg.norm(Xm - Xn) / np.linalg.norm(Xn)
L.check("modal sparse propagation == dense nodal ETD2RK (closed loop)", err_modal <= np.sqrt(tol), rel_err=err_modal)
# ---------------- (v) determinism
def det_run():
    s = build(cfg, train=True)
    sm = Simulator(s["lat"], s["rom"], s["policy"], cfg, adversary="worst")
    r = sm.run(X0, nmod, record="substep")
    return array_hash(r["h"], r["X_final"], *s["qnet"].Wq)


h1, h2 = det_run(), det_run()
L.check("bit-for-bit determinism: two runs give identical SHA-256", h1 == h2, sha256=h1)
write_table("exp8_numerics", ["check", "value"], [[r"spectral radius $\max|\lambda|$", rad], ["explicit Euler amplification at conv. step", amp_euler], ["stiffness ratio", rad / np.abs(ev).min()],
            ["latch: rel. error exact", err_exact], ["latch: rel. error grid-rounded", err_grid], ["modal vs nodal rel. error", err_modal],
            ["ETD2RK observed orders", "$" + ",\\ ".join(f"{o:.3f}" for o in order) + "$"], ["determinism SHA-256 (prefix)", r"\texttt{" + h1[:16] + "}"]],
            caption="Numerical-integrity checks.")
fig, ax = pl.figure(1, 1, 3.5, 2.4)
st = [Tc / s for s in subs]
ax[0, 0].loglog(st, errs, "o-", label="ETD2RK")
ax[0, 0].loglog(st, errs1, "s-", label="ETD1")
ax[0, 0].loglog(st, errs[0] * (np.array(st) / st[0]) ** 2, "k:", label="slope 2")
ax[0, 0].set_xlabel(r"step $T_s/n_{\mathrm{sub}}$"); ax[0, 0].set_ylabel("relative error"); ax[0, 0].legend(loc="lower right")
ax[0, 0].set_title("Step-halving convergence")
ax[0, 0].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax[0, 0].set_xticks(st, [f"{x:.3g}" for x in st])
pl.save(fig, FIG, "exp8_convergence")
L.data.update(errs=errs, errs1=errs1, err_exact=err_exact, err_grid=err_grid, err_modal=err_modal, sha=h1)
L.finish()
