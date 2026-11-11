"""Experiment 9 - disclosed sweep of the cutoff-mode advection weight w (S1; W_f = Phi_f diag(1,..,1,w,..,w) Phi_f^T).

For every w on the geometric grid experiments.w_sweep (w = 1 is the plain projector Pi_f) at the R2 Peclet number:
certificate accept/refuse with the failing hypotheses, eps_basis_P, the tail-coupling share of delta_total,
||A_{adv,r,perp}||_2, the (H6) margin, the certificate-level refusal of the ablations (a), (b), (c), (e), and the
open-loop exit of the full-order plant from the deployed R2 initial set.  The deployed w is the largest grid value
whose full certificate accepts with a positive (H6) margin while every ablation is still refused; the configuration
literal must equal it.  A second scan of the Peclet number (advection speed c) at w = 1 asserts that the
certifiable-c set and the open-loop-unstable-c set are disjoint (the Euclidean-ODD supremum defining M_perp makes
full-strength cutoff-mode advection uncertifiable wherever the open loop is unstable, Remark rem:rda_cutoff)."""
import numpy as np
from _common import write_csv, Ledger, load_config, override, write_table, FIG
from pgcbf.certificates import build, certificate, THEOREM_HYPS
from pgcbf.lattice import Lattice
from pgcbf.rom import ROM
from pgcbf.simulate import Simulator, initial_conditions
from pgcbf import plotting as pl
from pgcbf.texnames import hyp_list, num

L = Ledger("exp9_wsweep")
base = load_config("default")
ex, sc = base["experiments"], base["sim"]
s0 = build(base, train=True)  # Teacher / Student / LUT do not depend on the advection (P, g_r are advection-free)
c0 = certificate(base, s0)
X0 = initial_conditions(s0["lat"], s0["rom"], c0, int(sc["n_seeds"]), np.random.default_rng(int(base["seeds"]["ic"])),
                        float(sc["ic_level_fraction"]), float(sc["ic_tail_fraction"]), int(sc["ic_n_tail_modes"]))
Ts = float(base["timing"]["T_s"])
TsO = float(ex["wsweep_openloop_T_s"])
fa, fb = (float(x) for x in ex["multi_rate_fractions"])


def cert(cfg):
    lat = Lattice(cfg)
    return certificate(cfg, dict(s0, lat=lat, rom=ROM(lat, int(cfg["rom"]["K_m"])))), lat


def failing(c):
    bad = [h for h in THEOREM_HYPS if not c[h]]
    return bad + ([] if c["codesign_window"] else ["codesign_window"]) + ([] if c["timing_ok"] else ["timing_ok"])


def ablations_refused(cfg):
    """Certificate-level refusal of ablations (a) K_m too small, (b) coarse b_B, (c) single-rate, (e) large eps_adv."""
    out = {"a": not certificate(override(cfg, rom={"K_m": int(ex["ablation_K_m"])}), train=True)["theorem_valid"],
           "b": not certificate(override(cfg, controller={"b_B": int(ex["ablation_b_B"])}), train=True)["theorem_valid"],
           "c": not cert(override(cfg, timing={"tau_A": Ts, "tau_B": Ts}))[0]["theorem_valid"],
           "e": not cert(override(cfg, safety={"eps_adv": float(ex["ablation_eps_adv"])}))[0]["theorem_valid"]}
    return out


def open_loop_exit(cfg):
    lat = Lattice(cfg)
    rom = ROM(lat, int(cfg["rom"]["K_m"]))
    sim = Simulator(lat, rom, s0["policy"], cfg, T_s=TsO, tau_A=fa * TsO, tau_B=fb * TsO, adversary=None, controller=False)
    r = sim.run(X0, int(round(float(sc["horizon"]) / TsO)), record="period")
    return int(np.sum(r["h"].min(axis=0) < 0)), float(np.max(np.linalg.eigvals(lat.A_lin).real))


rows, recs = [], []
for w in (float(x) for x in ex["w_sweep"]):
    cfg = override(base, plant={"adv_cutoff_weight": w})
    c, lat = cert(cfg)
    ok = c["theorem_valid"]
    ab = ablations_refused(cfg) if ok else {k: True for k in "abce"}
    n_exit, growth = open_loop_exit(cfg)
    rec = dict(w=w, valid=ok, failing=failing(c), eps_b=c["eps_basis_P"], share=c["tail_share"], A=c["A_adv_rperp_norm"],
               Lp=c["L_perp"], margin=c["H6_margin"], abl=ab, n_exit=n_exit, growth=growth,
               delta_total=c["delta_total"], delta_cert=c["delta_cert"])
    recs.append(rec)
    rows.append([w, ok, hyp_list(rec["failing"]),
                 c["eps_basis_P"], c["tail_share"], c["A_adv_rperp_norm"], c["H6_margin"],
                 all(ab.values()) if ok else None, "%d/%d" % (n_exit, X0.shape[1])])
    print(f"  w={w:.4g} valid={ok} failing={rec['failing']} share={c['tail_share']:.3f} margin={c['H6_margin']:.3g} "
          f"ablations refused={ab} exits={n_exit}", flush=True)

ws = np.array([r["w"] for r in recs])
L.check("w grid is geometric, contains w = 1 and reaches down to w <= 0.005 with >= 8 values",
        ws.size >= 8 and ws.max() == 1 and ws.min() <= ex["w_sweep"][-1] and np.allclose(ws[1:] / ws[:-1], ws[1] / ws[0]),
        n=ws.size, w_min=ws.min())
L.check("plain projector w = 1 is refused at the R2 Peclet number", not recs[0]["valid"], failing=recs[0]["failing"])
L.check("open loop exits C for every seed at every w (the plant is unsafe for every swept w)",
        all(r["n_exit"] == X0.shape[1] for r in recs), exits=[r["n_exit"] for r in recs])
L.check("||A_adv,r,perp||_2 > 0 at every w (K_f > K_m) and proportional to w",
        all(r["A"] > 0 for r in recs) and np.allclose([r["A"] / r["w"] for r in recs], recs[0]["A"] / recs[0]["w"]))
cand = [r for r in recs if r["valid"] and r["margin"] > 0 and all(r["abl"].values())]
L.check("a certifiable frontier exists: some grid w accepts with (H6) margin > 0 and all ablations refused", len(cand) > 0)
wstar = max(r["w"] for r in cand)
acc = [r["w"] for r in recs if r["valid"]]
L.check("certifiable set is a down-set of the grid (every w below the frontier accepts)", all(w in acc for w in ws[ws <= min(acc)]) and
        all(r["valid"] == (r["w"] <= max(acc)) for r in recs), w_accept_max=max(acc))
L.check("deployed config weight equals the sweep-selected w* (largest admissible grid value)",
        base["plant"]["adv_cutoff_weight"] == wstar, w_star=wstar, w_config=base["plant"]["adv_cutoff_weight"])
rstar = next(r for r in recs if r["w"] == wstar)
L.check("at w*, the tail-coupling share of delta_total is material (>= wstar_min_tail_share)", rstar["share"] >= float(ex["wstar_min_tail_share"]),
        share=rstar["share"])
nxt = [r for r in recs if r["w"] > wstar]
L.data.update(w_star=wstar, w_star_share=rstar["share"], w_star_margin=rstar["margin"], w_accept_max=max(acc),
              w_first_refused_above=min(r["w"] for r in nxt) if nxt else None,
              records=[{k: v for k, v in r.items()} for r in recs])
write_csv("exp9_wsweep_deltas", ["w", "delta_total", "delta_cert", "theorem_valid", "failing", "w_star"],
          [[r["w"], r["delta_total"], r["delta_cert"], r["valid"], " ".join(r["failing"]) or "none", r["w"] == wstar] for r in recs])
write_table("exp9_wsweep", ["$w$", "theorem_valid", "failing", "eps_basis_P", "tail share", "A_adv_rperp_norm", "margin",
                            "abl.\\ refused", "exits"], rows,
            caption="Disclosed sweep of the cutoff-mode advection weight $w$ at the R2 P\\'eclet number ($w=1$: plain "
                    "projector $\\Pi_f$). The deployed $w^\\star$ is the largest grid value accepted with a positive "
                    "authority margin in \\eqref{eq:control_authority} and all ablations refused. Tail share: "
                    "$L_hL_{\\perp,P}\\varepsilon_{\\mathrm{basis},P}/\\delta_{\\mathrm{total}}$; margin: "
                    "$\\delta_{\\mathrm{cert}}-\\delta_{\\mathrm{total}}$; exits: open-loop exits of $%d$ seeds at $c=%s$; --: ablations not evaluated (the full certificate already refuses)." % (X0.shape[1], num(base["plant"]["c"]).strip("$")))
# ------------------------------------------------------------------ Peclet scan at w = 1 (plain projector)
crow, cert_c, unst_c = [], [], []
for cc in (float(x) for x in ex["c_sweep_plain"]):
    cfg = override(base, plant={"adv_cutoff_weight": 1 + 0 * cc, "c": cc})
    c, lat = cert(cfg)
    growth = float(np.max(np.linalg.eigvals(lat.A_lin).real))
    (cert_c if c["theorem_valid"] else []).append(cc)
    (unst_c if growth > 0 else []).append(cc)
    crow.append([cc, lat.peclet, c["theorem_valid"], hyp_list(failing(c)),
                 c["eps_basis_P"], c["tail_share"], growth])
L.check("plain projector (w = 1): certifiable-c set and open-loop-unstable-c set are disjoint and both non-empty",
        len(cert_c) > 0 and len(unst_c) > 0 and max(cert_c) < min(unst_c) and not set(cert_c) & set(unst_c),
        max_certifiable_c=max(cert_c) if cert_c else None, min_unstable_c=min(unst_c) if unst_c else None)
L.data.update(plain_max_certifiable_c=max(cert_c), plain_min_unstable_c=min(unst_c))
write_table("exp9_plain_peclet", ["$c$", "Pe", "theorem_valid", "failing", "eps_basis_P", "tail_share", "max Re eig$(A_{\\mathrm{full}}+A_{\\mathrm{adv}})$"], crow,
            caption="Plain projector $\\Pi_f$ ($w=1$): certificate acceptance versus open-loop instability over the "
                    "advection speed $c$. The two sets are disjoint.")
# ------------------------------------------------------------------ figure
fig, ax = pl.figure(2, 1, 3.58, 3.0)
a = ax[0, 0]
a.loglog(ws, [r["eps_b"] for r in recs], "o-", ms=3, label=r"$\varepsilon_{\mathrm{basis},P}$")
a.loglog(ws, [r["A"] for r in recs], "s-", ms=3, label=r"$\|A_{\mathrm{adv},r,\perp}\|_2$")
a.loglog(ws, [r["share"] for r in recs], "^-", ms=3, label="tail share")
a.axvline(wstar, color="k", ls=":", lw=0.8, label=r"$w^\star$")
a.set_xlabel(r"cutoff-mode weight $w$"); a.set_title("Tail constants vs. $w$"); a.legend(loc="upper center", bbox_to_anchor=(1 / 2, -0.25), ncol=2, columnspacing=1)
a = ax[0, 1]
mg = np.array([r["margin"] for r in recs])
a.semilogx(ws, mg, "o-", ms=3, label=r"$\delta_{\mathrm{cert}}-\delta_{\mathrm{total}}$")
pl.boundary(a, label="zero margin")
a.axvline(wstar, color="k", ls=":", lw=0.8, label=r"$w^\star$")
a.set_xlabel(r"cutoff-mode weight $w$"); a.set_ylabel("margin"); a.set_yscale("symlog", linthresh=1e-3)
a.set_title("Control-authority margin vs. $w$"); a.legend(loc="upper center", bbox_to_anchor=(1 / 2, -0.25), ncol=2, columnspacing=1)
pl.save(fig, FIG, "exp9_wsweep")
L.finish()
