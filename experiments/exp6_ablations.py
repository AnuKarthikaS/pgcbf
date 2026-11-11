"""Experiment 6 - ablations (a), (b), (c), (e) on the R2 benchmark (worst-case adversary).

Each ablation removes one component; each must produce a safety violation of the full-order
plant or a certificate failure (both are asserted where they occur):
 (a) no ROM tail bookkeeping: K_m too small vs. the tail criterion (L_perp eps_basis_P dropped);
 (b) no CROWN quantization margin: coarse b_B with the g_h(w2 eps_A + (1-w2) eps_CROWN) term dropped;
 (c) no multi-rate hold: single-rate tau_A = tau_B = T_s (weights w1 = 1, w2 = w3 = 0);
 (e) adversary above the certified eps_adv."""
import numpy as np
from _common import Ledger, load_config, override, write_table, FIG
from pgcbf.certificates import build, certificate, refuse_if_invalid, CertificateRefused
from pgcbf.simulate import Simulator, initial_conditions
from pgcbf import plotting as pl
from pgcbf.texnames import hyp_list, num

L = Ledger("exp6_ablations")
base = load_config("default")
ex, sc = base["experiments"], base["sim"]
s0 = build(base, train=True)
c0 = certificate(base, s0)
X0 = initial_conditions(s0["lat"], s0["rom"], c0, int(sc["n_seeds"]), np.random.default_rng(int(base["seeds"]["ic"])),
                        float(sc["ic_level_fraction"]), float(sc["ic_tail_fraction"]), int(sc["ic_n_tail_modes"]))
H = float(ex["ablation_horizon"])


def refuses(c):
    try:
        refuse_if_invalid(c)
        return False, ""
    except CertificateRefused as e:
        msg = str(e)
        return True, [h for h in ("H0", "H1", "H2", "H3", "H4", "H5", "H6", "H7", "timing", "codesign window") if "'%s'" % h in msg]


def simulate(cfg, substeps=None, record="substep"):
    s = build(cfg, train=True)
    sim = Simulator(s["lat"], s["rom"], s["policy"], cfg, adversary="worst", substeps=substeps)
    r = sim.run(X0, int(round(H / sim.Ts)), record=record)
    return r, r["h"].min(axis=0)


rows, curves = [], {}
# ------------------------------------------------------------------ (a)
cfa = override(base, rom={"K_m": int(ex["ablation_K_m"])})
ca = certificate(cfa, train=True)
cna = certificate(cfa, train=True, ablate=("tail",))
ref, why = refuses(ca)
L.check("(a) full certificate refuses K_m too small (tail criterion)", ref, reason=why, eps_basis_P=ca["eps_basis_P"],
        dt_tail=ca["dt_tail"], delta_total=ca["delta_total"], delta_cert=ca["delta_cert"])
L.check("(a) naive certificate without tail terms would pass (H4) and (H6)", cna["H4"] and cna["H6"],
        naive_delta_total=cna["delta_total"], naive_delta_cert=cna["delta_cert"], naive_H3=cna["H3"])
r, mh = simulate(cfa)
L.check("(a) K_m too small: full-order safety violated (some seed exits C)", np.any(mh < 0), n_exit=int(np.sum(mh < 0)), min_h=mh.min())
rows.append(["(a) $K_m=%d$, no tail bookkeeping" % ex["ablation_K_m"], ref, hyp_list(why), ca["delta_total"], ca["delta_cert"], cna["delta_total"],
             cna["delta_cert"], int(np.sum(mh < 0)), mh.min()])
curves["(a)"] = r
# ------------------------------------------------------------------ (b)
cfb = override(base, controller={"b_B": int(ex["ablation_b_B"])})
cb = certificate(cfb, train=True)
cnb = certificate(cfb, train=True, ablate=("quant",))
ref, why = refuses(cb)
L.check("(b) full certificate refuses coarse b_B", ref, reason=why, eps_CROWN=cb["eps_CROWN"], dt_quant=cb["dt_quant"])
L.check("(b) naive certificate without the quantization term would accept", cnb["theorem_valid"],
        naive_delta_total=cnb["delta_total"], delta_cert=cnb["delta_cert"])
r, mh = simulate(cfb)
L.check("(b) coarse b_B: full-order safety violated", np.any(mh < 0), n_exit=int(np.sum(mh < 0)), min_h=mh.min())
rows.append(["(b) $b_B=%d$, no CROWN margin" % ex["ablation_b_B"], ref, hyp_list(why), cb["delta_total"], cb["delta_cert"], cnb["delta_total"],
             cnb["delta_cert"], int(np.sum(mh < 0)), mh.min()])
curves["(b)"] = r
# ------------------------------------------------------------------ (c)
Ts = float(base["timing"]["T_s"])
cfc = override(base, timing={"tau_A": Ts, "tau_B": Ts})
cc = certificate(cfc, train=True)
ref, why = refuses(cc)
L.check("(c) single-rate at the deployed T_s: certificate refuses (lag term w1 = 1)", ref, reason=why,
        delta_total_single=cc["delta_total"], delta_total_multi=c0["delta_total"], delta_cert=cc["delta_cert"])
L.check("(c) multi-rate hold strictly reduces delta_total and keeps (H6)", c0["H6"] and c0["theorem_valid"] and c0["delta_total"] < cc["delta_total"])
wrows = [["multi rate", c0["tau_A"], c0["tau_B"], c0["w1"], c0["w2"], c0["Delta_h_dip"], c0["dt_lag"], c0["delta_total"], c0["H6"]],
         ["single rate", cc["tau_A"], cc["tau_B"], cc["w1"], cc["w2"], cc["Delta_h_dip"], cc["dt_lag"], cc["delta_total"], cc["H6"]]]
TsV = float(ex["single_rate_Ts_violation"])
fa, fb = (float(x) for x in ex["multi_rate_fractions"])
sub = [int(x) for x in ex["violation_substeps"]]
r1, mh1 = simulate(override(base, timing={"T_s": TsV, "tau_A": TsV, "tau_B": TsV}), substeps=sub)
r2, mh2 = simulate(override(base, timing={"T_s": TsV, "tau_A": fa * TsV, "tau_B": fb * TsV}), substeps=sub)
L.check(f"(c) single-rate at T_s={TsV}: full-order safety violated", np.any(mh1 < 0), n_exit=int(np.sum(mh1 < 0)), min_h=mh1.min())
L.check(f"(c) multi-rate at the same T_s={TsV}: no exit (outside certificate scope, observed)", np.all(mh2 >= 0),
        n_exit=int(np.sum(mh2 < 0)), min_h=mh2.min())
for lab, cfx in (("single", override(base, timing={"T_s": TsV, "tau_A": TsV, "tau_B": TsV})),
                 ("multi", override(base, timing={"T_s": TsV, "tau_A": fa * TsV, "tau_B": fb * TsV}))):
    cx = certificate(cfx, train=True)
    wrows.append([f"{lab} rate, $T_s={TsV:g}$", cx["tau_A"], cx["tau_B"], cx["w1"], cx["w2"], cx["Delta_h_dip"], cx["dt_lag"], cx["delta_total"], cx["H6"]])
rows.append(["(c) single rate $\\tau_A=\\tau_B=T_s$", ref, hyp_list(why), cc["delta_total"], cc["delta_cert"], None, None, int(np.sum(mh1 < 0)), mh1.min()])
curves["(c) single"] = r1
curves["(c) multi"] = r2
write_table("exp6_multirate", ["schedule", "tau_A", "tau_B", "w1", "w2", "Delta_h_dip", "dt_lag", "delta_total", "H6"], wrows,
            caption="Ablation (c): exponential weights $w_1,w_2$, $\\Delta h_{\\rm dip}$ and $\\delta_{\\rm total}$ for multi-rate vs.\\ single-rate hold.")
# ------------------------------------------------------------------ (e)
cfe = override(base, safety={"eps_adv": float(ex["ablation_eps_adv"])})
ce = certificate(cfe, train=True)
ref, why = refuses(ce)
L.check("(e) certificate refuses eps_adv above the certified level", ref, reason=why, eps_adv=ex["ablation_eps_adv"])
lo, hi = float(base["safety"]["eps_adv"]), float(ex["ablation_eps_adv"])
for _ in range(int(base["stochastic"]["bisection_iters"]) // 3):
    mid = (lo + hi) / 2
    lo, hi = (mid, hi) if certificate(override(base, safety={"eps_adv": mid}), s0)["theorem_valid"] else (lo, mid)
eps_star = lo
r, mh = simulate(cfe)
L.check(f"(e) adversary eps_adv={ex['ablation_eps_adv']}: full-order safety violated", np.any(mh < 0), n_exit=int(np.sum(mh < 0)), min_h=mh.min())
rs, mhs = simulate(override(base, safety={"eps_adv": eps_star}))
L.check("(e) at the largest certifiable eps_adv* the worst-case adversary does not break safety", np.all(mhs >= 0),
        eps_adv_star=eps_star, min_h=mhs.min())
rows.append([f"(e) $\\epsilon_{{\\rm adv}}={ex['ablation_eps_adv']:g}$", ref, hyp_list(why), ce["delta_total"], ce["delta_cert"], None, None, int(np.sum(mh < 0)), mh.min()])
curves["(e)"] = r
L.data.update(eps_adv_star=eps_star)
write_table("exp6_ablations", ["ablation", "refused", "failing", "delta_total", "delta_cert", "naive $\\delta_{\\rm total}$", "naive $\\delta_{\\rm cert}$",
                               "exits (of $%d$)" % X0.shape[1], "$\\min_t h$"], rows,
            caption="Ablations at the deployed $w^\\star$: every removed component yields a certificate refusal and a full-order "
                    "safety violation under the worst-case adversary. Row (e): the largest certifiable adversary level is "
                    "$\\epsilon_{\\rm adv}^\\star=%s$." % num(eps_star).strip("$"))
fig, ax = pl.figure(2, 3, 3.58, 1.75)
TsC = num(TsV).strip("$")
TT = {"(a)": "(a) $K_m$ too small", "(b)": "(b) coarse $b_B$", "(c) single": "(c1) single rate ($T_s=%s$)" % TsC,
      "(c) multi": "(c2) multi rate ($T_s=%s$)" % TsC, "(e)": "(e) $\\epsilon_{\\rm adv}$ too large"}
for i, (k, r) in enumerate(curves.items()):
    a = ax[i // 2, i % 2]
    for j in range(r["h"].shape[1]):
        a.plot(r["t"], r["h"][:, j], lw=0.4)
    pl.boundary(a)
    a.set_title(TT[k] + ": raw $h$")
    a.set_xlabel("$t$"); a.set_ylabel("$h(x(t))$")
ax[2, 1].axis("off")
pl.save(fig, FIG, "exp6_ablations")
L.finish()
