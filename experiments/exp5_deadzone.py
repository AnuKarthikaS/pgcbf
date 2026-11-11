"""Experiment 5 - Proposition prop:deadzone_stability (practical stability to B_ss) and
ablation (d): no deadzone / (H-dz) violated.  Config R1 (Pe = 1, mu_R > L_sec), w = 0.

(i)  Full method (mid-tread LUT with deadzone, deployed b_A): ||z(t)||_P <= max(||z(t0)||_P, r_ss)
     for all t, trajectories enter B_ss(r_ss) subset B_int (r_ss < r_min) and the control
     becomes exactly zero (no limit cycle).
(ii) Ablation (d1), S7: at the deployed b_A and across the sweep experiments.deadzone_b_A_sweep, mid-tread with
     deadzone vs. mid-rise (no zero level).  Mid-rise never reaches u_A = 0 and keeps switching sign (chattering);
     the deadzone variant reaches u_A = 0 exactly.
(iii) Ablation (d2): b_A = 1 -> eps_A large -> r_ss >= r_min, (H-dz) fails and the certificate
     refuses Prop. 3 (prop3_valid false)."""
import numpy as np
from _common import Ledger, load_config, override, write_table, FIG
from pgcbf.certificates import build, certificate
from pgcbf.simulate import Simulator, initial_conditions
from pgcbf import plotting as pl

L = Ledger("exp5_deadzone")
base = load_config("r1_passive")
ex = base["experiments"]
sc = base["sim"]
nseed = int(sc["n_seeds"]) // int(sc["ic_n_tail_modes"])
b_dep = int(base["controller"]["b_A"])
sweep = [int(b) for b in ex["deadzone_b_A_sweep"]]
variants = [(f"deadzone, $b_A={b}$", {"b_A": b}, b, True) for b in sweep] + \
           [(f"mid-rise, $b_A={b}$", {"b_A": b, "deadzone": False}, b, False) for b in sweep] + [("deadzone, $b_A=1$", {"b_A": 1}, 1, True)]
rows, res = [], {}
X0 = None
Ts = float(base["timing"]["T_s"])
for label, ov, b, dzf in variants:
    cfg = override(base, controller=ov)
    sysd = build(cfg, train=True)
    c = certificate(cfg, sysd)
    lat, rom, pol = sysd["lat"], sysd["rom"], sysd["policy"]
    if X0 is None:
        X0 = initial_conditions(lat, rom, c, nseed, np.random.default_rng(int(cfg["seeds"]["ic"])),
                                float(sc["ic_level_fraction"]), float(sc["ic_tail_fraction"]), int(sc["ic_n_tail_modes"]))
    sim = Simulator(lat, rom, pol, cfg, adversary=None)
    n = int(round(float(ex["deadzone_horizon"]) / sim.Ts))
    r = sim.run(X0, n, record="period")
    Zs = r["z_samples"]  # (n, K, P)
    zn = np.sqrt(np.einsum("nip,ij,njp->np", Zs, rom.P, Zs))
    T = np.stack([pol.teacher(Zs[k]) for k in range(Zs.shape[0])])
    lagA = np.linalg.norm(r["uA"] - T, axis=1)
    win = int(round(float(ex["deadzone_tail_window"]) / sim.Ts))
    uA_late = r["uA"][-win:]
    nonzero_frac = float(np.mean(np.any(uA_late != 0, axis=1)))
    sw = np.sum(np.diff(np.sign(uA_late), axis=0) != 0, axis=(0, 1)) / (win * sim.Ts)  # sign switches per unit time
    eff = np.sum(uA_late ** 2, axis=(0, 1)) * sim.Ts
    res[label] = dict(cert=c, zn=zn, t=np.arange(zn.shape[0]) * sim.Ts, uA=r["uA"], lag=lagA, nonzero=nonzero_frac, sw=sw,
                      eff=eff, b=b, dz=dzf, win=win)
    rows.append([("deadzone" if dzf else "mid-rise"), b, c["eps_A"], c["r_ss"], c["r_min"], c["Hdz"], zn[-win:].max(),
                 lagA[-1].min(), sw.mean(), r["h"].min()])

full = res[f"deadzone, $b_A={b_dep}$"]
c = full["cert"]
L.check(f"full method (deadzone, b_A={b_dep}): prop3_valid ((H-dz) r_ss < r_min and mu_R > L_sec)", c["prop3_valid"],
        r_ss=c["r_ss"], r_min=c["r_min"])
bound = np.maximum(full["zn"][0], c["r_ss"])
L.check("Prop3(b): ||z(t_k)||_P <= max(||z(0)||_P, r_ss) for all k, all seeds", np.all(full["zn"] <= bound * (1 + float(base["tolerances"]["roundoff_rel"]))),
        max_ratio=(full["zn"] / bound).max())
L.check("Prop3: trajectories enter B_ss(r_ss) subset B_int", np.all(full["zn"][-1] <= c["r_ss"]), final_max=full["zn"][-1].max(), r_ss=c["r_ss"])
L.check("deadzone: Path-A control exactly zero in the final window (no limit cycle)", full["nonzero"] == 0, nonzero_frac=full["nonzero"])
rate = c["mu_R"] - c["L_sec"]
for b in sweep:
    dz, mr = res[f"deadzone, $b_A={b}$"], res[f"mid-rise, $b_A={b}$"]
    tagb = " (deployed b_A)" if b == b_dep else ""
    L.check(f"(d) b_A={b}{tagb}: deadzone -> u_A exactly 0 late, no switching, zero late effort",
            dz["nonzero"] == 0 and np.all(dz["sw"] == 0) and np.all(dz["eff"] == 0), nonzero=dz["nonzero"], switches=dz["sw"].mean())
    L.check(f"(d) b_A={b}{tagb}: mid-rise (no deadzone) -> u_A never zero, persistent sign switching (chattering), effort > 0",
            mr["nonzero"] == 1 and np.all(mr["sw"] > 0) and np.all(mr["eff"] > 0), nonzero=mr["nonzero"],
            switches_per_time=mr["sw"].mean(), effort=mr["eff"].min())
    dec = dz["zn"][-1] / dz["zn"][-dz["win"]]
    L.check(f"(d) b_A={b}: deadzone -> once u = 0, ||z||_P contracts at least at rate mu_R - L_sec (Prop3(a))",
            np.all(dec <= np.exp(-rate * dz["win"] * Ts)), factor=dec.max(), predicted=np.exp(-rate * dz["win"] * Ts))
    L.check(f"(d) b_A={b}: mid-rise lag error |u_A - T(z)| persists at the LSB scale", mr["lag"][-1].min() > 0, lag=mr["lag"][-1].min())
swr = np.array([res[f"mid-rise, $b_A={b}$"]["sw"].mean() for b in sweep])
L.data.update(midrise_switches_per_time=dict(zip(map(str, sweep), swr)), deployed_b_A=b_dep,
              midrise_late_max_z={str(b): float(res[f"mid-rise, $b_A={b}$"]["zn"][-full["win"]:].max()) for b in sweep},
              deadzone_late_max_z={str(b): float(res[f"deadzone, $b_A={b}$"]["zn"][-full["win"]:].max()) for b in sweep})
b1 = res["deadzone, $b_A=1$"]["cert"]
L.check("d2: b_A = 1 violates (H-dz): r_ss >= r_min, certificate refuses Prop. 3 (prop3_valid false)",
        (not b1["Hdz"]) and (not b1["prop3_valid"]), r_ss=b1["r_ss"], r_min=b1["r_min"])
# ------------------------------------------------------------------ figures: b_A sweep (S7)
fig, ax = pl.figure(2, 1, 3.58, 2.9)
a = ax[0, 0]
a.semilogy(sweep, swr, "o-", ms=3, color="C1", label="mid-rise (no deadzone)")
a.set_xlabel(r"$b_A$"); a.set_ylabel("late sign switches per unit time"); a.set_xticks(sweep)
a.set_title(r"(d) chattering vs. $b_A$ (deadzone: 0)"); a.legend(loc="upper center", bbox_to_anchor=(1 / 2, -0.25))
a = ax[0, 1]
for dzf, col, lab in ((False, "C1", "mid-rise"), (True, "C0", "deadzone")):
    pre = "deadzone" if dzf else "mid-rise"
    a.semilogy(sweep, [res[f"{pre}, $b_A={b}$"]["zn"][-full["win"]:].max() for b in sweep], "o-" if dzf else "s--", ms=3, color=col,
               label=lab + ": late max")
a.semilogy(sweep, [res[f"deadzone, $b_A={b}$"]["cert"]["r_ss"] for b in sweep], "k:", lw=0.8, label=r"$r_{\mathrm{ss}}(b_A)$")
a.set_xlabel(r"$b_A$"); a.set_xticks(sweep); a.set_ylabel(r"$\|z(t_k)\|_P$ (raw)")
a.set_title(r"late state size vs. $b_A$"); a.legend(loc="upper center", bbox_to_anchor=(1 / 2, -0.25))
pl.save(fig, FIG, "exp5_deadzone_sweep")
fig, ax = pl.figure(2, 1, 3.58, 2.2)
for i, (pre, col) in enumerate((("deadzone", "C0"), ("mid-rise", "C1"))):
    a = ax[0, i]
    for b, ls in zip(sweep, ("-", "--", "-.", ":")):
        R = res[f"{pre}, $b_A={b}$"]
        a.semilogy(R["t"], np.maximum(R["zn"][:, 0], np.finfo(float).tiny), ls=ls, lw=0.7, label=f"$b_A={b}$")
    a.set_xlabel(r"$t$"); a.set_ylabel(r"$\|z(t_k)\|_P$ (seed 0, raw)"); a.set_title(f"{pre}: state norm, $b_A$ sweep")
    a.legend(loc="best")
pl.save(fig, FIG, "exp5_deadzone_z")
write_table("exp5_deadzone", ["LUT", "$b_A$", "eps_A", "r_ss", "r_min", "Hdz", "late $\\max\\|z\\|_P$", "final lag $\\|u_A-\\hat{\\mathcal{T}}\\|_2$",
                              "switches/time", "$\\min h$"], rows,
            caption="Proposition~\\ref{prop:deadzone_stability} and ablation (d) on R1 (disturbance $w\\equiv0$): mid-tread LUT with deadzone vs.\\ "
                    "mid-rise LUT (no zero level) at the deployed $b_A=%d$ and over the sweep $b_A\\in\\{%s\\}$; the last row is "
                    "$b_A=1$, which violates (H-dz). Final lag: $\\min$ over seeds of $\\|u_A(t_k)-\\hat{\\mathcal{T}}(z(t_k))\\|_2$ "
                    "at the last sample. For every swept $b_A$, $u_A\\equiv0$ in the late window with the deadzone and "
                    "$u_A\\neq0$ at every late sample without it (asserted)." % (b_dep, ",".join(map(str, sweep))))
L.finish()
