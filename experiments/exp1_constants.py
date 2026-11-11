"""Experiment 1 - every certificate constant computed at run time, hypothesis ledger, Prop. 2
identities, Lemma quad_split, and the certificate's refusal on non-admissible configurations.

Configs: R2 (default.json, Pe = 20, open-loop unstable) and R1 (r1_passive.json, Pe = 1), both with the advection band
K_f = n_adv_modes > K_m (cutoff-mode weight w from the w-sweep of exp9).  Gate G-D (spec correction): A_full stays the
pure damped-wave operator (Q = gamma P, (H2a) by eigenmode invariance, A_{r,perp} = 0); the linear cross-coupling is
A_{adv,r,perp} = P0^{-1} Phi^T P_full A_adv P_perp (Lemma lem:delta_cert_classH), which must be nonzero in R2, and L_perp_P(R2) must be
materially larger than L_perp_P(R1)."""
import numpy as np
from scipy.linalg import expm
from _common import write_csv, Ledger, load_config, override, write_table
from pgcbf.texnames import tex, hyp_list, HYP_DESC, num
from pgcbf.certificates import (build, certificate, refuse_if_invalid, CertificateRefused, pencil_rate, mu_2P, port_box,
                                THEOREM_HYPS)

L = Ledger("exp1_constants")


def _fmtpe(x):
    return ("%.3g" % x)

cfgs = {"R2": load_config("default"), "R1": load_config("r1_passive")}
certs, systems = {}, {}
for name, cfg in cfgs.items():
    sysd = build(cfg, train=True)
    c = certificate(cfg, sysd)
    certs[name], systems[name] = c, sysd
    tol = float(cfg["tolerances"]["roundoff_rel"])
    lat, rom = sysd["lat"], sysd["rom"]
    # ---- Proposition prop:dahlquist (eq:mu_pe)
    L.check(f"{name}: Prop2 mu_2P(A_r0) = -gamma/2", abs(c["mu_2P_Ar0"] + lat.gamma / 2) <= tol * lat.gamma * lat.N,
            mu_2P=c["mu_2P_Ar0"], minus_gamma_half=-lat.gamma / 2)
    L.check(f"{name}: Prop2 Q = gamma P identity", c["Q_identity_relerr"] <= tol * lat.N, rel_err=c["Q_identity_relerr"])
    # full-order mu at eps = gamma/2 equals -gamma/2 as well (all modes underdamped)
    muF = mu_2P(lat.Af, lat.Pfe)
    L.check(f"{name}: Prop2 full-order mu_2,Pfe(A_full) = -gamma/2", abs(muF + lat.gamma / 2) <= tol * lat.gamma * lat.N * lat.omega2.max(),
            mu_full=muF)
    # tail rate lambda_{K+1} vs modal pencil formula (Remark rem:tail_rates) at eps = gamma/2
    pr = pencil_rate(lat.eps, lat.omega[rom.K_m:], lat.gamma).min()
    L.check(f"{name}: tail rate lambda_K+1 equals modal pencil rate", abs(c["lambda_K1"] - pr) <= tol * lat.N * lat.omega2.max(),
            lambda_K1=c["lambda_K1"], pencil=pr)
    # tail decay: unforced, advection-free lattice from a pure tail state decays exactly at gamma/2
    rng = np.random.default_rng(int(cfg["seeds"]["ic"]))
    xp = rom.Pperp @ rng.standard_normal(lat.N)
    ts = np.linspace(0, float(cfg["sim"]["horizon"]), int(cfg["sim"]["n_seeds"]))
    ratio = [np.sqrt((expm(lat.Af * t) @ xp) @ lat.Pfe @ (expm(lat.Af * t) @ xp) / (xp @ lat.Pfe @ xp)) / np.exp(-lat.gamma * t / 2) for t in ts]
    L.check(f"{name}: tail energy decays exactly as exp(-gamma t/2) (Q = gamma P)", np.max(np.abs(np.array(ratio) - 1)) <= np.sqrt(tol),
            max_dev=np.max(np.abs(np.array(ratio) - 1)))
    # ---- Lemma lem:quad_split under (H2a)
    E = c["E_safe"]
    Z = rng.standard_normal((rom.K, int(cfg["stochastic"]["n_paths"])))
    Z = Z / rom.Pnorm(Z) * np.sqrt(E * rng.random(Z.shape[1]))  # inside C_r0
    Xp = rom.Pperp @ rng.standard_normal((lat.N, Z.shape[1]))
    Xp = Xp / rom.tail_Pnorm(Xp) * c["eps_basis_P"] * rng.random(Z.shape[1])  # ||x_perp|| <= eps_basis_P
    lhs = lat.h_barrier(rom.Phi @ Z + Xp, E)
    rhs = rom.hPhi(Z, E) - np.einsum("ip,ij,jp->p", Xp, lat.Pfe, Xp)
    scale = np.abs(lhs).max()
    L.check(f"{name}: Lemma quad_split h(Phi z + x_perp) = h(Phi z) - ||x_perp||^2 (rel 1e-12)",
            np.max(np.abs(lhs - rhs)) <= tol * scale, max_abs_err=np.max(np.abs(lhs - rhs)), scale=scale)
    L.check(f"{name}: (H2a) Phi^T P_eps P_perp = 0", c["H2a_residual"] <= float(cfg["tolerances"]["h2a"]) * np.linalg.norm(lat.Pfe, 2),
            residual=c["H2a_residual"])
    L.check(f"{name}: Q = gamma P relative error < 1e-12 with K_f = {lat.n_f} > K_m = {rom.K_m} (A_full advection-free)",
            c["Q_identity_relerr"] < float(cfg["tolerances"]["h2a"]) and lat.n_f > rom.K_m, rel_err=c["Q_identity_relerr"])
    L.check(f"{name}: A_r,perp (A_full-based) vanishes by eigenmode invariance", c["A_rperp_full_norm"] <= tol * lat.N * lat.omega2.max(),
            norm=c["A_rperp_full_norm"])
    # S3: port box closed form D_v = 2 sqrt(E) max_j sqrt((g_r^T P g_r)_jj) vs sampled sup over C_r0
    G = rng.standard_normal((rom.K, int(cfg["experiments"]["crown_samples"])))
    Zb = G / rom.Pnorm(G) * np.sqrt(E)
    Dv_s = 2 * np.abs(rom.g_r.T @ rom.P @ Zb).max()
    gPg = rom.g_r.T @ rom.P @ rom.g_r
    j = int(np.argmax(np.diag(gPg)))
    zstar = rom.g_r[:, j] * np.sqrt(E / gPg[j, j])  # Cauchy-Schwarz maximiser on the boundary of C_r0
    Dv_att = 2 * np.abs(rom.g_r.T @ rom.P @ zstar).max()
    L.check(f"{name}: S3 D_v closed form = attained sup (round-off) >= sampled sup (within 1%)",
            abs(Dv_att - c["D_lut"]) <= tol * lat.N * c["D_lut"] and Dv_s <= c["D_lut"] * (1 + tol) and Dv_s >= c["D_lut"] * (1 - 1 / 100),
            D_v=c["D_lut"], D_v_attained=Dv_att, D_v_sampled=Dv_s)
    c["D_v_sampled"] = Dv_s
    for hyp in THEOREM_HYPS + ("H_sat", "timing_ok", "codesign_window"):
        L.check(f"{name}: hypothesis {hyp}", c[hyp])
    L.check(f"{name}: theorem_valid (all of (H0)-(H7)) and the certificate does not refuse", c["theorem_valid"] and refuse_if_invalid(c))
L.check("R1: prop3_valid ((H-dz): r_ss < r_min and mu_R > L_sec)", certs["R1"]["prop3_valid"] and certs["R1"]["Hdz"],
        r_ss=certs["R1"]["r_ss"], r_min=certs["R1"]["r_min"])
L.check("R2: (H-dz) honestly fails (L_sec > mu_R: open-loop unstable advection); Prop. 3 not claimed", not certs["R2"]["prop3_valid"],
        L_sec=certs["R2"]["L_sec"], mu_R=certs["R2"]["mu_R"])
# ---- gate G-D (S1): advective cross-coupling active in the tail
r2, r1 = certs["R2"], certs["R1"]
L.check("G-D: ||A_adv,r,perp||_2 > 0 in R2 (K_f > K_m)", r2["A_adv_rperp_norm"] > 0 and r2["K_f"] > systems["R2"]["rom"].K_m,
        A_adv_rperp=r2["A_adv_rperp_norm"], K_f=r2["K_f"])
ratio = r2["L_perp"] / r1["L_perp"]
L.check("G-D: L_perp_P(R2) materially larger than L_perp_P(R1) (ratio >= Pe ratio / 2)", ratio >= (r2["Pe"] / r1["Pe"]) / 2,
        L_perp_R2=r2["L_perp"], L_perp_R1=r1["L_perp"], ratio=ratio)
L.check("G-D: L_perp_P = linear advective part + sin-nonlinearity part (both reported)",
        all(abs(cc["L_perp"] - cc["L_perp_lin"] - cc["L_perp_nl"]) <= float(cfgs["R2"]["tolerances"]["roundoff_rel"]) * cc["L_perp"]
            for cc in (r1, r2)), lin_R2=r2["L_perp_lin"], nl_R2=r2["L_perp_nl"])
L.check("G-D: tail-coupling share of delta_total reported and material in R2 (>= 10%)", r2["tail_share"] >= 1 / 10,
        share_R2=r2["tail_share"], share_R1=r1["tail_share"])
L.data.update(L_perp_ratio=ratio)

# ---- refusals (config adjustments recorded in README)
base = cfgs["R2"]
refusals = [
    ("state-space LUT ($m = K$ variant)", override(base, controller={"lut_space": "state"})),
    ("$u_{\\max}$ below the (H-sat) left side", override(base, controller={"u_max": certs["R2"]["H_sat_lhs"] * base["experiments"]["u_max_refusal_fraction"]})),
    ("pure lattice advection (no spectral cutoff)", override(base, plant={"n_adv_modes": base["plant"]["N_m"], "adv_cutoff_weight": 1})),
    ("plain projector $\\Pi_f$ ($w = 1$)", override(base, plant={"adv_cutoff_weight": 1})),
    ("previous sampling period $T_s=%s$ at $w^\\star$" % num(base["experiments"]["wave2_T_s"]).strip("$"),
     override(base, timing={k: v * base["experiments"]["wave2_T_s"] / base["timing"]["T_s"] for k, v in base["timing"].items()})),
    ("$T_s$, $\\tau_A$, $\\tau_B$ $\\times %d$ (hold term)" % base["experiments"]["Ts_refusal_factor"], override(base, timing={k: v * base["experiments"]["Ts_refusal_factor"] for k, v in base["timing"].items()})),
]
ref_rows = []
for label, cf in refusals:
    cc = certificate(cf, train=False)
    try:
        refuse_if_invalid(cc)
        refused, msg = False, ""
    except CertificateRefused as e:
        refused, msg = True, str(e)
    bad = [h for h in THEOREM_HYPS if not cc[h]] + (["codesign_window"] if not cc["codesign_window"] else []) \
        + (["timing_ok"] if not cc["timing_ok"] else [])
    L.check(f"certificate refuses: {label}", refused, reason=msg)
    ref_rows.append([label, refused, hyp_list(bad), cc["delta_total"], cc["delta_cert"], "[%d, %d]" % (cc["codesign_lo"], cc["codesign_hi"])])
write_table("exp1_refusals", ["configuration", "refused", "failing", "delta_total", "delta_cert", "window $[\\mathrm{lo},\\mathrm{hi}]$"],
            ref_rows, caption="Certificate refusals on non-admissible variants of R2. The LUT window is "
            "\\eqref{eq:codesign_window} in the variant's own form (port space by default; the state-space row uses "
            "$K$, $D_{\\mathrm{ODD}}$ and $M_z$).")

# ---- constants table
rows = [
    ("N_m", "--"), ("N", "--"), ("m", "--"), ("nu", "--"), ("alpha_plant", "--"), ("gamma_plant", "--"), ("kappa", "--"),
    ("c", "--"), ("K_m", "--"), ("K", "--"), ("n_adv_full_weight", "--"), ("E_safe", "eq:H0"), ("eps_adv", "eq:M_perp"), ("open_loop_growth", "--"), ("K_f", "--"), ("adv_cutoff_weight", "--"), ("Pe", "rem:rda_cutoff"),
    ("lambda_min_Pfe", "eq:eps_choice"), ("omega_min", "eq:eps_choice"),
    ("mu_2P_Ar0", "eq:mu_pe"), ("Q_identity_relerr", "eq:mu_pe"), ("lambda_K1", "eq:tail_dissipation"),
    ("H2a_residual", "lem:h2a_char"), ("A_rperp_full_norm", "lem:h2a_char"),
    ("rho_ODD", "eq:H0"), ("D_ODD", "eq:Cr0"), ("L_Fnl", "assump:generic"),
    ("M_perp", "eq:M_perp"), ("eps_basis_P", "eq:eps_basis_P"), ("A_adv_rperp_norm", "lem:delta_cert_classH"),
    ("L_perp_lin", "eq:Lperp_def"), ("L_perp_nl", "eq:Lperp_def"), ("L_perp", "eq:Lperp_def"), ("L_perp_lin_share", "eq:Lperp_def"),
    ("L_h", "eq:reduced_constants"), ("g_h", "eq:reduced_constants"), ("d_h", "eq:reduced_constants"),
    ("M_z", "eq:reduced_constants"), ("L_Lfh", "eq:reduced_constants"), ("L_Lgh", "eq:reduced_constants"),
    ("L_Tz", "eq:LTz"), ("eta_T", "lem:sector"), ("H_sat_lhs", "eq:sat_cond"), ("u_max", "eq:sat_cond"),
    ("eps_CROWN_backward", "eq:twin_crown_bound"), ("E_q", "eq:crown_explicit"), ("eps_CROWN", "prop:twin_crown"),
    ("D_v", "eq:Dv_closed"), ("D_v_sampled", "eq:Dv_closed"), ("M_v", "eq:eps_A_triangle"), ("Delta_qA", "eq:Dv_closed"),
    ("eps_A", "eq:eps_A_triangle"), ("du_lag", "eq:eps_A_triangle"),
    ("codesign_lo", "eq:codesign_window"), ("b_A", "eq:codesign_window"), ("codesign_hi", "eq:codesign_window"),
    ("T_s", "eq:multirate_hold"), ("tau_A", "eq:multirate_hold"), ("tau_B", "eq:multirate_hold"), ("samples_per_period", "eq:multirate_hold"), ("b_B", "prop:twin_crown"), ("gamma0", "eq:teacher_cbf"), ("w1", "eq:exp_weights"), ("w2", "eq:exp_weights"), ("w3", "eq:exp_weights"),
    ("phi", "eq:delta_total_exact"), ("L_hold", "eq:delta_dip"), ("Delta_h_dip", "eq:delta_dip"),
    ("dt_tail", "eq:delta_total_exact"), ("dt_quant", "eq:delta_total_exact"), ("dt_lag", "eq:delta_total_exact"),
    ("dt_adv", "eq:delta_total_exact"), ("dt_hold", "eq:delta_total_exact"), ("delta_total", "eq:delta_total_exact"),
    ("tail_share", "eq:delta_total_exact"), ("hold_share", "eq:delta_total_exact"),
    ("h_act", "eq:shell"), ("r_min", "eq:shell"), ("lambda_min_MP", "lem:delta_cert_classH"), ("s_nl", "lem:delta_cert_classH"),
    ("c_M", "lem:delta_cert_classH"), ("delta_cert", "lem:delta_cert_classH"), ("delta_cert_branch", "eq:delta_cert_classH"),
    ("H6_margin", "eq:control_authority"), ("gamma0_h_act", "cor:stochastic_w2"),
    ("null_space_rate", "eq:null_space_cond"), ("mu_R", "prop:deadzone_stability"),
    ("L_sec", "prop:deadzone_stability"), ("eps_q_bar", "prop:deadzone_stability"), ("r_ss", "eq:rss_ball"),
]
for n in certs:
    cc, latn = certs[n], systems[n]["lat"]
    plc = cfgs[n]["plant"]
    cc.update(N_m=int(plc["N_m"]), nu=float(plc["nu"]), alpha_plant=float(plc["alpha"]), gamma_plant=float(plc["gamma"]),
              kappa=float(plc["kappa"]), c=float(plc["c"]), K_m=int(cfgs[n]["rom"]["K_m"]), n_adv_full_weight=int(plc["n_adv_full_weight"]))
    cc["D_v"] = cc["D_lut"]
    cc["Delta_qA"] = cc["D_lut"] / 2 ** cc["b_A"]
    cc["open_loop_growth"] = float(np.linalg.eigvals(latn.A_lin).real.max())  # max Re eig of the open-loop linearisation
    cc["L_perp_lin_share"] = cc["L_perp_lin"] / cc["L_perp"]
    cc["samples_per_period"] = 2 * np.pi / cc["omega_min"] / cc["T_s"]  # samples per slowest (fundamental) period
    cc["hold_share"] = cc["dt_hold"] / cc["delta_total"]
    cc["gamma0_h_act"] = cc["gamma0"] * cc["h_act"]  # pointwise penalty gamma0 h_act of Corollary cor:stochastic_w2 at sigma = 0
    cc["mu_R_gt_L_sec"] = bool(cc["mu_R"] > cc["L_sec"])
L.check("R1: Prop. 3(i) contraction mu_R > L_sec holds", certs["R1"]["mu_R_gt_L_sec"], mu_R=certs["R1"]["mu_R"], L_sec=certs["R1"]["L_sec"])
L.check("R2: Prop. 3(i) fails (L_sec > mu_R), so Proposition 3 does not apply on R2", not certs["R2"]["mu_R_gt_L_sec"],
        mu_R=certs["R2"]["mu_R"], L_sec=certs["R2"]["L_sec"])
L.check("R2: gamma0 h_act exceeds the (H6) margin delta_cert - delta_total (Corollary 1 infeasible on R2 even at sigma = 0)",
        certs["R2"]["gamma0_h_act"] > certs["R2"]["H6_margin"], gamma0_h_act=certs["R2"]["gamma0_h_act"], margin=certs["R2"]["H6_margin"])
# budget split of eq:delta_total_exact for the main-paper Fig. 2 (quantization term g_h (w2 eps_A + (1-w2) eps_CROWN)
# split into the CROWN part g_h eps_CROWN and the LUT-rounding part g_h w2 (eps_A - eps_CROWN))
BUDGET = ("tail", "LUT", "CROWN", "lag", "hold", "adversary")
brows = []
for n in ("R1", "R2"):
    cc = certs[n]
    parts = [cc["dt_tail"], cc["g_h"] * cc["w2"] * (cc["eps_A"] - cc["eps_CROWN"]), cc["g_h"] * cc["eps_CROWN"], cc["dt_lag"],
             cc["dt_hold"], cc["dt_adv"]]
    L.check(f"{n}: the six budget parts (tail, LUT, CROWN, lag, hold, adversary) sum to delta_total",
            abs(sum(parts) - cc["delta_total"]) <= float(cfgs[n]["tolerances"]["roundoff_rel"]) * cc["delta_total"],
            parts_sum=sum(parts), delta_total=cc["delta_total"])
    brows.append([n] + parts + [cc["delta_total"], cc["delta_cert"]])
write_csv("exp1_budget", ["config"] + list(BUDGET) + ["delta_total", "delta_cert"], brows)
hdr = ["quantity", "equation", "R2 (Pe $= %s$)" % _fmtpe(certs["R2"]["Pe"]), "R1 (Pe $= %s$)" % _fmtpe(certs["R1"]["Pe"])]
tab = [[tex(k), lab, certs["R2"][k], certs["R1"][k]] for k, lab in rows]
half = (len(tab) + 1) // 2
write_table("exp1_constants", hdr, tab[:half], wide=True,
            caption="Certificate constants (part 1 of 2), computed at run time in float64 from the lattice operators. "
                    "$K_f$ is the advection band, $w$ the disclosed cutoff-mode weight (Table~\\ref{tab:exp9_wsweep}); "
                    "$L_{\\perp,P}(\\mathrm{R2})/L_{\\perp,P}(\\mathrm{R1})=%s$." % num(ratio).strip("$"))
write_table("exp1_constants_b", hdr, tab[half:], wide=True, label="tab:exp1_constants_b",
            caption="Certificate constants (part 2 of 2).")
hyps = list(THEOREM_HYPS) + ["H_sat", "timing_ok", "codesign_window", "theorem_valid", "mu_R_gt_L_sec", "Hdz", "prop3_valid"]
write_table("exp1_hypotheses", ["hypothesis", "content", "R2", "R1"],
            [[tex(h), HYP_DESC.get(h, ""), certs["R2"][h], certs["R1"][h]] for h in hyps], wide=True, align="llcc",
            caption="Hypothesis ledger of Theorem~\\ref{thm:main_cbf}, hypotheses (H0)--(H7), and of "
                    "Proposition~\\ref{prop:deadzone_stability}, conditions (i) $\\mu_R>L_{\\mathrm{sec}}$ and (H-dz); \\checkmark{} = satisfied.")
L.data.update(constants={n: {k: certs[n][k] for k, _ in rows} for n in certs}, hypotheses={n: {h: certs[n][h] for h in hyps} for n in certs},
              train_hist={n: (None if systems[n]["train_hist"] is None else systems[n]["train_hist"][[0, -1]]) for n in systems})
L.finish()
