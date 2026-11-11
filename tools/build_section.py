"""V4: build results/experimental_section.tex from results/experimental_section_template.tex.

Every number in the prose is a {{key}} placeholder resolved from the logs of the final run (results/logs/*.json),
the generated CSV tables (results/tables/*.csv) and the configuration files, so Section V cannot drift from the
tables.  The resolved values are written to results/section_numbers.json (raw floats and their LaTeX rendering).
Exit status 1 on any unresolved placeholder, unused key or texlint finding."""
import csv
import json
import math
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "tools"))
from pgcbf.texnames import num  # noqa: E402
from texlint import lint_files  # noqa: E402

LOG, TAB, RES = (os.path.join(ROOT, "results", d) for d in ("logs", "tables", ""))


def log(name):
    with open(os.path.join(LOG, f"{name}.json")) as f:
        return json.load(f)


def cfg(name):
    with open(os.path.join(ROOT, "config", f"{name}.json")) as f:
        return json.load(f)


def table(name):
    with open(os.path.join(TAB, f"{name}.csv")) as f:
        return list(csv.reader(f))


def assertion(d, needle):
    """The unique assertion record of log d whose text contains needle."""
    hits = [a for a in d["assertions"] if needle in a["assertion"]]
    if len(hits) != 1:
        raise KeyError(f"{d['experiment']}: {len(hits)} assertions match {needle!r}")
    return hits[0]


def f(x):
    """The ONE shared number formatter of the prose (review items 7, 21): pgcbf.texnames.num with its default
    significant-digit rule -- the same function, with the same arguments, that typesets every float cell of the
    generated tables -- so prose and tables agree to the digit.  Integers are printed exactly."""
    if isinstance(x, bool):
        raise TypeError("booleans are not prose numbers")
    if isinstance(x, int):
        return str(x)
    return num(float(x)).strip("$")


V, RAW = {}, {}


def put(key, x, raw=None):
    V[key] = x if isinstance(x, str) else f(x)
    RAW[key] = raw if raw is not None else (x if not isinstance(x, str) else None)


def collect():
    R2c, R1c = cfg("default"), cfg("r1_passive")
    pl2, pl1 = R2c["plant"], R1c["plant"]
    for k in ("N_m", "n_adv_modes", "n_adv_full_weight"):
        put(k, int(pl2[k]))
    put("N", 2 * int(pl2["N_m"]))
    for k in ("nu", "kappa", "alpha", "gamma"):
        put(k, float(pl2[k]))
    put("c_R2", float(pl2["c"])); put("c_R1", float(pl1["c"]))
    put("K_m", int(R2c["rom"]["K_m"])); put("K", 2 * int(R2c["rom"]["K_m"]))
    put("w_star", float(pl2["adv_cutoff_weight"]))
    put("w_star_exp", "2^{-11/2}" if abs(pl2["adv_cutoff_weight"] - 2 ** (-11 / 2)) < 1e-15 else f(pl2["adv_cutoff_weight"]))
    put("Ts_R2", float(R2c["timing"]["T_s"])); put("Ts_R1", float(R1c["timing"]["T_s"])); put("Ts_wave2", float(R2c["experiments"]["wave2_T_s"]))
    put("tauA_R2", float(R2c["timing"]["tau_A"])); put("tauB_R2", float(R2c["timing"]["tau_B"]))
    # main paper, benchmark paragraph (reviewer request): configuration constants from config/*.json
    put("tauA_R1", float(R1c["timing"]["tau_A"])); put("tauB_R1", float(R1c["timing"]["tau_B"]))
    put("cfg_Esafe", float(R2c["safety"]["E_safe"])); put("cfg_gamma0_R2", float(R2c["safety"]["gamma0"])); put("cfg_gamma0_R1", float(R1c["safety"]["gamma0"]))
    put("cfg_umax_R2", float(R2c["controller"]["u_max"])); put("cfg_umax_R1", float(R1c["controller"]["u_max"]))
    put("cfg_LT_R2", float(R2c["controller"]["L_T"])); put("cfg_LT_R1", float(R1c["controller"]["L_T"]))
    put("cfg_thetaGS", float(R2c["controller"]["theta"])); put("cfg_H", int(R2c["controller"]["H"])); put("cfg_L", int(R2c["controller"]["L"]))
    put("cfg_Bmax", int(R2c["hardware"]["B_max"])); put("cfg_Sbits_log2", int(R2c["hardware"]["S_bits_log2"])); put("cfg_bout", int(R2c["hardware"]["b_out"]))
    put("cfg_horizon", float(R2c["sim"]["horizon"])); put("cfg_ablation_horizon", float(R2c["experiments"]["ablation_horizon"]))
    put("eps_adv", float(R2c["safety"]["eps_adv"]))
    put("n_seeds", int(R2c["sim"]["n_seeds"]))
    put("r_W2", float(R1c["stochastic"]["r_W2"]))
    # ---- exp0
    e0 = log("exp0_openloop")["data"]
    put("Pe_R2", e0["Pe"]); put("e0_exits", sum(1 for m in e0["min_h"] if m < 0)); put("e0_minh", min(e0["min_h"]))
    put("e0_lastexit", max(e0["first_exit"]))
    # ---- exp1
    d1 = log("exp1_constants")
    C2, C1 = d1["data"]["constants"]["R2"], d1["data"]["constants"]["R1"]
    put("Pe_R1", C1["Pe"]); put("m", int(C2["m"]))
    assert abs(C2["open_loop_growth"] - e0["growth"]) <= 1e-12 * abs(e0["growth"]), "exp0 / exp1 growth rates differ"
    for k in ("eps_basis_P", "A_adv_rperp_norm", "L_perp", "L_perp_lin", "L_perp_nl", "L_perp_lin_share", "tail_share",
              "delta_total", "delta_cert", "H6_margin", "c_M", "lambda_min_MP", "s_nl", "dt_hold", "M_perp", "L_sec", "mu_R",
              "r_min", "r_ss", "Delta_h_dip", "Q_identity_relerr", "H2a_residual", "A_rperp_full_norm", "lambda_K1", "D_v",
              "D_v_sampled", "omega_min", "samples_per_period", "hold_share", "gamma0_h_act", "open_loop_growth", "L_Lfh",
              "L_hold"):
        for tag, C in (("R2", C2), ("R1", C1)):
            v = C[k]
            put(f"{k}_{tag}", float(v) if not isinstance(v, str) else math.inf)
    put("codesign_lo", int(C2["codesign_lo"])); put("codesign_hi", int(C2["codesign_hi"])); put("b_A", int(C2["b_A"]))
    put("mu2P", C2["mu_2P_Ar0"])
    # story bible Section 9: budget utilisation delta_total / delta_cert, three significant digits
    for tag, Cx in (("R2", C2), ("R1", C1)):
        u = Cx["delta_total"] / Cx["delta_cert"]
        put(f"util_{tag}", num(u, 3).strip("$"), raw=u)
    put("Lperp_ratio", d1["data"]["L_perp_ratio"])
    rr = {r[0]: r for r in table("exp1_refusals")[1:]}
    w2 = [r for k, r in rr.items() if k.startswith("previous sampling period")][0]
    put("wave2_delta_total", float(w2[3])); put("wave2_delta_cert", float(w2[4]))
    put("n_refusals", len(rr)); put("Ts_refusal_factor", int(R2c["experiments"]["Ts_refusal_factor"]))
    ss = [r for k, r in rr.items() if k.startswith("state-space")][0]
    put("ss_window", ss[5].replace("[", "").replace("]", "").replace(", ", ","))
    # ---- exp2
    e2 = log("exp2_crown")["data"]["sweep"]
    ratios = [s["eps_crown"] / s["max_gap"] for s in e2]
    put("e2_ratio_lo", min(ratios)); put("e2_ratio_hi", max(ratios))
    put("e2_bB_lo", min(s["b_B"] for s in e2)); put("e2_bB_hi", max(s["b_B"] for s in e2))
    bB = int(R2c["controller"]["b_B"])
    put("e2_bB", bB)
    put("e2_near", [s for s in e2 if s["b_B"] == bB and s["rounding"] == "nearest"][0]["eps_crown"])
    put("e2_learn", [s for s in e2 if s["b_B"] == bB and s["rounding"] == "learned"][0]["eps_crown"])
    nsamp = re.search(r"\((\d+) samples\)", log("exp2_crown")["assertions"][0]["assertion"]).group(1)
    put("e2_nsamp", int(nsamp))
    # ---- exp3
    e3 = log("exp3_delta_cert")["data"]
    for tag in ("R2", "R1"):
        put(f"e3_bound_{tag}", e3[tag]["delta_cert"]); put(f"e3_brute_{tag}", e3[tag]["brute"])
    # ---- exp4
    e4 = log("exp4_theorem1")["data"]
    for tag in ("R2", "R1"):
        b, i = e4[f"{tag} boundary"], e4[f"{tag} interior"]
        put(f"e4_bmin_{tag}", b["min_h"]); put(f"e4_imin_{tag}", i["min_h"])
        Cx = C2 if tag == "R2" else C1
        put(f"e4_btmin_{tag}", b["t_min"]); put(f"e4_itmin_{tag}", i["t_min"])
        put(f"e4_tailratio_{tag}", max(b["max_tail"], i["max_tail"]) / Cx["eps_basis_P"])
    put("interior_level", float(R2c["experiments"]["interior_level"]))
    # ---- exp4b (supplement, Experiment 4b): realised fraction of each budget term of eq:ode_diff_ineq
    e4b = log("exp4b_budget")["data"]
    for tag in ("R2", "R1"):
        x = e4b[tag]
        st, sa = x["stats"], x["stats_active"]
        put(f"e4b_nodes_{tag}", int(x["nodes"])); put(f"e4b_act_{tag}", int(x["active_nodes"]))
        put(f"e4b_intfrac_{tag}", x["interior_nodes"] / x["nodes"])
        for t in ("tail", "adv", "quant", "hold"):
            put(f"e4b_{t}_max_{tag}", st[t]["max"]); put(f"e4b_{t}_med_{tag}", sa[t]["median"])
        put(f"e4b_eu_max_{tag}", st["f_eu"]["max"]); put(f"e4b_mz_max_{tag}", st["f_mz"]["max"])
        put(f"e4b_lipH_max_{tag}", st["f_lipH"]["max"]); put(f"e4b_lip_med_{tag}", sa["f_lip"]["median"])
        put(f"e4b_cos_med_{tag}", sa["f_cos"]["median"])
        put(f"e4b_minslack_{tag}", x["min_slack"] / x["delta_total"])
        decs = [e4b[f"{tag} {st_}"] for st_ in ("boundary", "interior")]
        worst = min(decs, key=lambda d: d["min_slack"])
        put(f"e4b_teachshare_{tag}", worst["min_slack_decomposition"]["teacher"] / worst["min_slack"])
        # prose claims of the Experiment 4b paragraph, checked against the data (fail rather than print a false sentence)
        assert x["min_slack"] > 0, f"exp4b {tag}: soundness slack not positive"
        assert st["hold"]["max"] == min(st[t]["max"] for t in ("tail", "adv", "quant", "hold")), f"exp4b {tag}: hold not smallest"
        assert worst["min_slack_where"]["s"] == 0, f"exp4b {tag}: minimum slack not at a sampling instant"
        assert max(sa["f_lip"]["median"], sa["f_cos"]["median"]) < min(sa["f_grad"]["median"], sa["f_env"]["median"]), \
            f"exp4b {tag}: L_perp and alignment are not the two least used tail factors"
    # ---- exp5
    d5 = log("exp5_deadzone")
    e5 = d5["data"]
    sw = e5["midrise_switches_per_time"]
    bAs = sorted(int(k) for k in sw)
    put("e5_bAs", ", ".join(str(b) for b in bAs))
    put("e5_bA_lo", bAs[0]); put("e5_bA_hi", bAs[-1])
    put("e5_sw_lo", sw[str(bAs[0])]); put("e5_sw_hi", sw[str(bAs[-1])])
    lag9 = assertion(d5, f"(d) b_A={bAs[-1]}: mid-rise lag error")["lag"]
    put("e5_lag_dep", lag9)
    dz1 = assertion(d5, "b_A = 1 violates (H-dz)")
    put("e5_rss1", dz1["r_ss"]); put("e5_rmin1", dz1["r_min"])
    put("e5_bA_dep", int(e5["deployed_b_A"]))
    # ---- exp6
    d6 = log("exp6_ablations")
    put("e6_a_exits", int(assertion(d6, "(a) K_m too small: full-order")["n_exit"]))
    put("e6_b_exits", int(assertion(d6, "(b) coarse b_B: full-order")["n_exit"]))
    put("e6_e_exits", int(assertion(d6, "(e) adversary eps_adv=")["n_exit"]))
    ca = assertion(d6, "(c) single-rate at the deployed T_s")
    put("e6_c_single", ca["delta_total_single"]); put("e6_c_multi", ca["delta_total_multi"])
    put("e6_c_delta", ca["delta_total_single"])  # story bible Section 9: single-rate delta_total at the deployed T_s
    # ---- main paper, Table 1: the six budget parts per regime from exp1_budget.csv (same CSV as Fig. 2 left)
    bt = table("exp1_budget")
    hdr = bt[0]
    for row in bt[1:]:
        reg = row[hdr.index("config")]
        for term in ("tail", "LUT", "CROWN", "lag", "hold", "adversary"):
            put(f"budget_{term}_{reg}", float(row[hdr.index(term)]))
            put(f"budget_{term}_share_{reg}", float(row[hdr.index(term)]) / float(row[hdr.index("delta_total")]))
    cv = [a for a in d6["assertions"] if a["assertion"].startswith("(c) single-rate at T_s=") and "violated" in a["assertion"]][0]
    put("e6_c_exits", int(cv["n_exit"])); put("e6_TsV", float(R2c["experiments"]["single_rate_Ts_violation"]))
    put("e6_eps_star", d6["data"]["eps_adv_star"])
    put("abl_K_m", int(R2c["experiments"]["ablation_K_m"])); put("abl_b_B", int(R2c["experiments"]["ablation_b_B"]))
    put("abl_eps_adv", float(R2c["experiments"]["ablation_eps_adv"]))
    # ---- exp7
    e7 = log("exp7_stochastic")["data"]
    s = e7["sigma_star"]
    gs = [r for r in e7["sigma_grid"] if r["sigma"] == s["sigma"]][0]  # screening record at sigma* (Table exp7_sigma_grid)
    put("e7_sigma", s["sigma"]); put("e7_exits", int(s["exits"])); put("e7_paths", int(s["paths"]))
    put("e7_cp_hi", s["cp"][1]); put("e7_alpha", s["alpha"]); put("e7_T", s["T"]); put("e7_TTL", s["T_TTL"])
    put("e7_TTL_ratio", s["T_TTL"] / s["T"]); put("e7_h0", s["h0"]); put("e7_meddev", gs["med_dev"])
    put("e7_noisedev", gs["med_ndev"])
    put("e7_z0", float(R1c["stochastic"]["z0_level"])); put("e7_ffrac", float(R1c["stochastic"]["fluct_fraction_of_h0"]))
    put("e7_cp_pct", 100 * float(R1c["stochastic"]["confidence"]))
    grid = [r["sigma"] for r in e7["sigma_grid"]]
    put("e7_grid", ", ".join(f(x) for x in grid))
    fr = min((r for r in e7["sigma_grid"] if not r["valid"]), key=lambda r: r["sigma"])
    put("e7_first_refused", fr["sigma"])
    put("e7_refused_cond", ", ".join({"margin": "margin", "beta0_gt_h0": r"\(\beta_0>h_0\)", "T_le_TTL": r"\(T\le T_{\mathrm{TTL}}\)"}
                                       .get(k, k.replace("_", " ")) for k in fr["fail"]))
    t = e7["noise_tail"]
    put("e7f_Dsig", t["full"]["Delta_sigma"]); put("e7f_sqrtE", math.sqrt(R1c["safety"]["E_safe"]))
    put("e7f_exits", int(t["exits"])); put("e7f_cp_lo", t["cp"][0])
    r2s = e7["R2_sigma0"]
    assert (not r2s["valid"]) and "margin" in r2s["failing"], r2s
    put("e7_R2_gamma0_h_act", r2s["gamma0_h_act"]); put("e7_R2_margin", r2s["margin"])
    # ---- exp8
    n8 = {r[0]: r[1] for r in table("exp8_numerics")[1:]}
    put("e8_amp", float(n8["explicit Euler amplification at conv. step"])); put("e8_rad", float([v for k, v in n8.items() if k.startswith("spectral radius")][0]))
    put("e8_stiff", float(n8["stiffness ratio"])); put("e8_latch", float(n8["latch: rel. error exact"]))
    put("e8_grid", float(n8["latch: rel. error grid-rounded"])); put("e8_modal", float(n8["modal vs nodal rel. error"]))
    e8 = log("exp8_numerics")["data"]
    r8 = [a / b for a, b in zip(e8["errs"], e8["errs"][1:])]
    put("e8_r_lo", min(r8)); put("e8_r_hi", max(r8))
    # ---- exp9
    e9 = log("exp9_wsweep")["data"]
    rec = e9["records"]
    put("e9_n", len(rec)); put("e9_wmin", min(r["w"] for r in rec)); put("e9_wstar", e9["w_star"])
    put("e9_wmax", max(r["w"] for r in rec))
    put("e9_first_ref", e9["w_first_refused_above"])
    big = [r["w"] for r in rec if set(r["failing"]) >= {"H4"}]
    put("e9_wH4", min(big))
    mid = [r["w"] for r in rec if r["failing"] == ["H6"]]
    put("e9_H6_lo", min(mid)); put("e9_H6_hi", max(mid))
    cacc = [r["w"] for r in rec if r["valid"] and not r["abl"]["c"]]
    put("e9_c_acc", max(cacc) if cacc else "none")
    put("e9_share_w1", rec[0]["share"]); put("e9_share_star", e9["w_star_share"])
    put("e9_plain_cert", e9["plain_max_certifiable_c"]); put("e9_plain_unst", e9["plain_min_unstable_c"])
    put("e9_plain_cert_Pe", e9["plain_max_certifiable_c"] / R2c["plant"]["nu"]); put("e9_plain_unst_Pe", e9["plain_min_unstable_c"] / R2c["plant"]["nu"])
    # ---- runtimes / assertion counts
    rt = table("runtimes")[1:]
    tot = [r for r in rt if r[0] == "total"][0]
    put("rt_total", float(tot[1])); put("rt_peak", float(tot[2]))
    n_ass = sum(len(log(os.path.splitext(x)[0])["assertions"]) for x in os.listdir(LOG) if x.startswith("exp") and x.endswith(".json"))
    put("n_assert", n_ass)
    m = re.search(r"(\d+) passed", open(os.path.join(LOG, "pytest.log")).read())
    put("n_tests", int(m.group(1)) if m else 0)


def main():
    collect()
    with open(os.path.join(RES, "experimental_section_template.tex")) as fh:
        tpl = fh.read()
    keys = set(re.findall(r"\{\{([A-Za-z0-9_]+)\}\}", tpl))
    missing = sorted(keys - set(V))
    unused = sorted(set(V) - keys)
    out = re.sub(r"\{\{([A-Za-z0-9_]+)\}\}", lambda m: V.get(m.group(1), m.group(0)), tpl)
    with open(os.path.join(RES, "experimental_section.tex"), "w") as fh:
        fh.write(out)
    with open(os.path.join(RES, "section_numbers.json"), "w") as fh:
        json.dump({k: {"tex": V[k], "raw": RAW[k] if not (isinstance(RAW[k], float) and not math.isfinite(RAW[k])) else str(RAW[k])}
                   for k in sorted(V)}, fh, indent=1)
    issues = lint_files([os.path.join(RES, "experimental_section.tex")])
    from texlint import unsourced_numerals
    issues += [f"numeral {t!r} in the prose does not occur in section_numbers.json" for t in unsourced_numerals(out, list(V.values()))]
    for x in issues:
        print(x)
    print(f"build_section: {len(keys)} placeholders, {len(missing)} unresolved {missing}, {len(unused)} unused keys, "
          f"{len(issues)} lint issues")
    sys.exit(1 if (missing or issues) else 0)


if __name__ == "__main__":
    main()
