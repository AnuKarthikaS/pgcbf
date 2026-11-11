"""Experiment 2 - Proposition prop:twin_crown: soundness of the twin-CROWN bound
(bound >= max |Psi_B - T| on >= 1e5 random samples of Z_ODD, plus the box vertices),
monotonicity in b_B, scaling ~ 2^{-b_B}, and learned (Adam) vs nearest rounding."""
import itertools
import numpy as np
from _common import Ledger, load_config, override, write_table, FIG
from pgcbf.certificates import build, certificate
from pgcbf.crown import TwinCrown, E_q
from pgcbf import plotting as pl

L = Ledger("exp2_crown")
cfg = load_config("default")
ex = cfg["experiments"]
n = int(ex["crown_samples"])
rng = np.random.default_rng(int(cfg["seeds"]["mc"]))


def true_gap(sysd, Z):
    pol = sysd["policy"]
    return np.linalg.norm(pol.psi_B(Z) - pol.teacher(Z), axis=0)


def samples(sysd, K):
    D = sysd["D_ODD"]
    Z = rng.uniform(-D / 2, D / 2, (K, n))
    V = np.array(list(itertools.product([-D / 2, D / 2], repeat=K))).T  # vertices
    return np.hstack([Z, V])


rows, sweep = [], []
for b in ex["crown_b_B_sweep"]:
    for rounding in ("nearest", "learned"):
        cf = override(cfg, controller={"b_B": int(b), "rounding": rounding})
        sysd = build(cf, train=True)
        rom = sysd["rom"]
        tc = TwinCrown(sysd["qnet"], sysd["policy"].portmap)
        eC, _ = tc.eps_crown(sysd["D_ODD"], rom.K, int(cf["crown"]["splits_per_dim"]))
        Eq = E_q(sysd["qnet"], sysd["Y"])
        Z = samples(sysd, rom.K)
        g = true_gap(sysd, Z)
        sweep.append({"b_B": int(b), "rounding": rounding, "eps_crown": eC, "E_q": Eq, "max_gap": g.max(), "mean_gap": g.mean()})
        L.check(f"b_B={b} {rounding}: twin-CROWN bound >= max sampled |Psi_B - T| ({Z.shape[1]} samples)", eC >= g.max(),
                eps_crown=eC, max_gap=g.max(), ratio=eC / g.max())
        rows.append([b, rounding, eC, Eq, g.max(), g.mean(), eC / g.max()])
for rounding in ("nearest", "learned"):
    e = [s["eps_crown"] for s in sweep if s["rounding"] == rounding]
    L.check(f"{rounding}: eps_CROWN monotone non-increasing in b_B", all(a >= b for a, b in zip(e, e[1:])), seq=e)
    bb = np.array([s["b_B"] for s in sweep if s["rounding"] == rounding], float)
    slope = np.polyfit(bb, np.log2(e), 1)[0]
    L.check(f"{rounding}: eps_CROWN scales as 2^(-b_B) (log2 slope within [-1.25,-0.75])",
            abs(slope + 1) <= 1 / 4, slope=slope)
# port-box (Path-A table entries) at the deployed b_B
sysd = build(cfg, train=True)
cert = certificate(cfg, sysd)
tcv = TwinCrown(sysd["qnet"], np.eye(sysd["lat"].m))
Dv = sysd["policy"].D_lut
V = rng.uniform(-Dv / 2, Dv / 2, (sysd["lat"].m, n))
from pgcbf.controller import project_U
gv = np.linalg.norm(project_U(sysd["qnet"].forward(V), sysd["policy"].u_max) - project_U(sysd["teacher"].forward(V), sysd["policy"].u_max), axis=0)
L.check("port box: twin-CROWN bound >= sampled gap on the Path-A table domain", cert["eps_CROWN_portbox"] >= gv.max(),
        eps=cert["eps_CROWN_portbox"], max_gap=gv.max())
nearest = [s for s in sweep if s["b_B"] == cfg["controller"]["b_B"] and s["rounding"] == "nearest"][0]
learned = [s for s in sweep if s["b_B"] == cfg["controller"]["b_B"] and s["rounding"] == "learned"][0]
L.check("learned rounding does not increase the certified eps_CROWN at the deployed b_B",
        learned["eps_crown"] <= nearest["eps_crown"], learned=learned["eps_crown"], nearest=nearest["eps_crown"])
L.data.update(sweep=sweep, deployed_eps_crown=cert["eps_CROWN"], deployed_eps_A=cert["eps_A"])
write_table("exp2_crown", ["b_B", "rounding", "eps_CROWN", "E_q", "max sampled gap", "mean gap", "bound / max gap"], rows,
            caption=f"Proposition~\\ref{{prop:twin_crown}}: twin-CROWN bound vs.\\ the sampled Student--Teacher gaps on "
                    f"${Z.shape[1]}$ points of $\\mathcal{{Z}}_{{\\mathrm{{ODD}}}}$ (${n}$ uniform samples plus all "
                    f"${Z.shape[1] - n}$ box vertices).")
fig, ax = pl.figure(1, 1, 3.5, 3.1)
for rounding, mk in (("nearest", "o"), ("learned", "s")):
    s = [x for x in sweep if x["rounding"] == rounding]
    ax[0, 0].semilogy([x["b_B"] for x in s], [x["eps_crown"] for x in s], mk + "-", label=f"twin-CROWN, {rounding}")
    ax[0, 0].semilogy([x["b_B"] for x in s], [x["max_gap"] for x in s], mk + ":", label=f"max gap, {rounding}")
    ax[0, 0].semilogy([x["b_B"] for x in s], [x["E_q"] for x in s], mk + "--", alpha=0.6, label=f"$E_q$, {rounding}")
ax[0, 0].set_xlabel(r"$b_B$"); ax[0, 0].set_ylabel(r"$\|\Psi_B(z)-\hat{\mathcal{T}}(z)\|_2$"); ax[0, 0].legend(loc="upper center", bbox_to_anchor=(1 / 2, -0.22), ncol=2)
ax[0, 0].set_title(r"Twin-CROWN soundness vs. bit-width $b_B$")
pl.save(fig, FIG, "exp2_crown")
L.finish()
