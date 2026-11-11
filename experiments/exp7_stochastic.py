"""Experiment 7 - Corollary cor:stochastic_w2 (Ito sampled-data safety, T_TTL) and ablation (f).

R1 benchmark, worst-case adversary eps_adv plus additive Brownian forcing Sigma dW, z(0) on the
level z0_level of C~_r (x_perp(0) = 0).
(i)  S5 positive case: sigma grid search.  sigma* is the largest grid value for which the full
     Corollary accepts at T = horizon AND the fluctuation criterion
         max_t |h(x_t) - h(x_0)| >= fluct_fraction_of_h0 * h(x_0)   on >= fluct_path_share of the paths
     holds (screened on n_screen_paths).  At sigma* a full Monte Carlo (n_paths) gives the exit
     frequency and an exact Clopper-Pearson interval, compared with the bound alpha.  The
     noise-induced part of the fluctuation (deviation from the noise-free trajectory) is reported
     separately, so the reader can see how much of the criterion is drift and how much is noise.
(ii) ablation (f): spatially white small-scale noise orthogonal to the ROM (Sigma_r = 0).  The
     full Corollary refuses (OU tail envelope Delta_sigma >> sqrt(E_safe)); a cruder bound that
     ignores the tail diffusion (sigma_perp := 0) accepts, and Monte Carlo violates its
     probability bound (Clopper-Pearson lower limit > alpha)."""
import numpy as np
from _common import Ledger, load_config, write_table, FIG
from pgcbf.certificates import build, certificate
from pgcbf.simulate import Simulator
from pgcbf.stochastic import noise_matrix, corollary, clopper_pearson
from pgcbf import plotting as pl
from pgcbf.texnames import num as _num


def num(x):
    """math-mode number without the surrounding dollars (for use inside $...$)."""
    return _num(x).strip("$")



L = Ledger("exp7_stochastic")
cfg = load_config("r1_passive")
st = cfg["stochastic"]
sysd = build(cfg, train=True)
c = certificate(cfg, sysd)
lat, rom, pol = sysd["lat"], sysd["rom"], sysd["policy"]
alpha, T, n, conf = float(st["alpha"]), float(st["T_max"]), int(st["n_paths"]), float(st["confidence"])
n_scr, ffrac, fshare = int(st["n_screen_paths"]), float(st["fluct_fraction_of_h0"]), float(st["fluct_path_share"])
COND = {"margin": "margin", "beta0_gt_h0": r"$\beta_0>h_0$", "h_act_minus_pos": r"$h_{\rm act}^{-}>0$",
        "r_min2_pos": r"$r_{\min}^2>0$", "T_le_TTL": r"$T\le T_{\rm TTL}$", "theorem": "Theorem 1"}

# deterministic z(0) on the level z0_level of C~_r (first ROM mode), x_perp(0) = 0
z0 = np.zeros(rom.K)
z0[0] = np.sqrt(float(st["z0_level"]) * (c["E_safe"] - c["eps_basis_P"] ** 2 - c["Delta_h_dip"]) / rom.P[0, 0])
x0 = rom.Phi @ z0
NS = int(round(T / Simulator(lat, rom, pol, cfg, adversary="worst").Ts))
det = Simulator(lat, rom, pol, cfg, adversary="worst").run(x0[:, None], NS, record="substep")
h_det = det["h"][:, 0]
h0 = float(h_det[0])
L.check("z(0) lies strictly inside C~_r and h(x_0) > 0", 0 < h0 and float(z0 @ rom.P @ z0) < c["E_safe"] - c["eps_basis_P"] ** 2,
        h0=h0, z0_level=st["z0_level"])


def mc(Sigma, npaths, seed_offset=0):
    sim = Simulator(lat, rom, pol, cfg, adversary="worst", noise=Sigma)
    rng = np.random.default_rng(int(cfg["seeds"]["mc"]) + seed_offset)
    return sim.run(np.repeat(x0[:, None], npaths, axis=1), NS, rng=rng, record="substep", stop_on_exit=True)


def fluct(r):
    """Fluctuation criterion: share of paths with max_t |h - h(x_0)| >= ffrac h(x_0); noise-induced deviation."""
    H = r["h"]
    dev = np.nanmax(np.abs(H - H[0]), axis=0)
    ndev = np.nanmax(np.abs(H - h_det[: H.shape[0], None]), axis=0)
    return float(np.mean(dev >= ffrac * h0)), float(np.median(dev / h0)), float(np.median(ndev / h0)), float(np.nanstd(H[-1]) / h0)


# ------------------------------------------------------------------ (i) sigma grid search
grid = sorted(float(s) for s in st["sigma_grid"])
grows, recs = [], []
for s in grid:
    Sig = noise_matrix(lat, rom, {"type": st["noise_small"]["type"], "sigma": s})
    o = corollary(cfg, c, sysd, Sig, T, alpha)
    fail = [k for k, v in o["conditions"].items() if not v]
    if o["valid"]:
        share, med_dev, med_ndev, end_std = fluct(mc(Sig, n_scr, seed_offset=1))
    else:
        share = med_dev = med_ndev = end_std = None
    lhs = o["delta_total_sigma"] + o["delta_stoch"] + o["delta_drift"]
    recs.append(dict(sigma=s, valid=o["valid"], fail=fail, share=share, med_dev=med_dev, med_ndev=med_ndev, end_std=end_std,
                     T_TTL=o["T_TTL"], lhs=lhs, dcert=o["delta_cert_sigma"]))
    grows.append([s, o["valid"], ", ".join(COND[k] for k in fail) if fail else "none", lhs, o["delta_cert_sigma"], o["T_TTL"],
                  share if o["valid"] else None, med_dev if o["valid"] else None, med_ndev if o["valid"] else None])
    print(f"  sigma={s:g} valid={o['valid']} fail={fail} share={share} med dev/h0={med_dev} "
          f"noise dev/h0={med_ndev} T_TTL={o['T_TTL']:.4g}", flush=True)
cand = [r for r in recs if r["valid"] and r["share"] is not None and r["share"] >= fshare]
L.check("sigma grid: some grid sigma is accepted by the Corollary at T and meets the fluctuation criterion", len(cand) > 0)
rs = max(cand, key=lambda r: r["sigma"])
sig_star = rs["sigma"]
above = [r for r in recs if r["sigma"] > sig_star]
L.check("sigma* is the largest grid value with Corollary acceptance and the fluctuation criterion (every larger grid value fails one)",
        all((not r["valid"]) or r["share"] is None or r["share"] < fshare for r in above), sigma_star=sig_star,
        first_refused=min((r["sigma"] for r in above), default=None))
L.check("the Corollary's acceptance region on the grid is a down-set (monotone in sigma)",
        all(r["valid"] == (r["sigma"] <= max(q["sigma"] for q in recs if q["valid"])) for r in recs))
write_table("exp7_sigma_grid", ["sigma", "Corollary", "failing condition", r"$\delta_{\rm total,\sigma}+\delta_{\rm stoch}+\delta_{\rm drift}$",
                                "delta_cert_sigma", "T_TTL", "fluct. share", "total dev.", "noise dev."], grows,
            caption=r"Grid search over the profile-noise intensity $\sigma$ (R1, $T=%s$, $\alpha=%s$). Fluct.\ share: fraction of "
                    r"%d screening paths with $\max_t|h(x_t)-h(x_0)|\ge %s\,h(x_0)$; noise dev.: median of "
                    r"$\max_t|h(x_t)-h_{\rm det}(t)|/h(x_0)$ (deviation from the noise-free run). $\sigma^\star=%s$."
                    % (num(T), num(alpha), n_scr, num(ffrac), num(sig_star)), wide=True)

# ------------------------------------------------------------------ full MC at sigma* and ablation (f)
rows, crow_d, crow_s, runs = [], [], [], {}
for key, ncfg in (("sigma_star", {"type": st["noise_small"]["type"], "sigma": sig_star}), ("noise_tail", st["noise_tail"])):
    Sigma = noise_matrix(lat, rom, ncfg)
    full = corollary(cfg, c, sysd, Sigma, T, alpha)
    crude = corollary(cfg, c, sysd, Sigma, T, alpha, ignore_tail_diffusion=True)
    r = mc(Sigma, n)
    runs[key] = (r, ncfg)
    exited = np.nanmin(r["h"], axis=0) < 0
    k = int(exited.sum())
    lo, hi = clopper_pearson(k, n, conf)
    nm = "profile, $\\sigma^\\star$" if key == "sigma_star" else "tail-white (f)"
    rows.append([nm, float(ncfg["sigma"]), h0, full["valid"], crude["valid"], full["T_TTL"], full["T_TTL"] / T, "$%d/%d$" % (k, n),
                 "$[%s,\\,%s]$" % (num(lo), num(hi))])
    for bn, o in (("full", full), ("crude $\\sigma_\\perp:=0$", crude)):
        lab = nm + ", " + bn
        crow_d.append([lab, o["sigma_r"], o["sigma_perp2"], o["Delta_sigma"], o["eb_sigma"], o["delta_total_sigma"], o["T_perp"]])
        crow_s.append([lab, o["r_BM"], o["delta_Ito"], o["sigma_h"], o["Delta_delta_pw"], o["delta_stoch"], o["c_int"]])
    if key == "sigma_star":
        share, med_dev, med_ndev, end_std = fluct(r)
        L.check("sigma*: Corollary conditions hold (margin, beta0 > h0, h_act^- > 0, T <= T_TTL, Theorem 1)", full["valid"],
                **{kk: bool(v) for kk, v in full["conditions"].items()})
        L.check("sigma*: fluctuation criterion holds on the full Monte Carlo (>= fluct_path_share of the paths)", share >= fshare,
                share=share, median_dev_over_h0=med_dev)
        L.check("sigma*: the noise itself is non-trivial: median max_t |h(x_t) - h_det(t)| >= fluct_fraction_of_h0 * h(x_0)",
                med_ndev >= ffrac, median_noise_dev_over_h0=med_ndev, end_std_over_h0=end_std)
        L.check("sigma*: MC exit frequency <= alpha (Clopper-Pearson upper limit <= alpha)", hi <= alpha,
                exits=k, paths=n, freq=k / n, cp_upper=hi, alpha=alpha)
        L.check("sigma*: the probability bound is conservative (T_TTL >> T and zero observed exits)", full["T_TTL"] > T and k == 0,
                T_TTL=full["T_TTL"], T=T, ratio=full["T_TTL"] / T)
        L.data["sigma_star"] = dict(sigma=sig_star, exits=k, paths=n, freq=k / n, cp=[lo, hi], alpha=alpha, T=T, T_TTL=full["T_TTL"],
                                    h0=h0, fluct_share=share, median_dev_over_h0=med_dev, median_noise_dev_over_h0=med_ndev,
                                    end_std_over_h0=end_std, min_h=float(np.nanmin(r["h"])), full=full)
    else:
        L.check("(f) tail-white noise: full Corollary refuses (OU tail envelope exceeds the safe set)", not full["valid"],
                Delta_sigma=full["Delta_sigma"], sqrt_E=np.sqrt(c["E_safe"]), T_perp=full["T_perp"])
        L.check("(f) crude bound ignoring tail diffusion would accept", crude["valid"], T_TTL_crude=crude["T_TTL"])
        L.check("(f) MC violates the crude probability bound (Clopper-Pearson lower limit > alpha)", lo > alpha,
                exits=k, paths=n, cp_lower=lo, alpha=alpha)
        L.data["noise_tail"] = {"exits": k, "paths": n, "cp": [lo, hi], "full": full, "crude": crude}
L.data["sigma_grid"] = recs
# ------------------------------------------------------------------ R2 at sigma = 0 (review item 12)
cfg2 = load_config("default")
sys2 = build(cfg2, train=True)
c2 = certificate(cfg2, sys2)
o2 = corollary(cfg2, c2, sys2, 0 * noise_matrix(sys2["lat"], sys2["rom"], {"type": st["noise_small"]["type"], "sigma": sig_star}),
               float(cfg2["stochastic"]["T_max"]), float(cfg2["stochastic"]["alpha"]))
g0h2, marg2 = c2["gamma0"] * c2["h_act"], c2["H6_margin"]
L.check("R2, sigma = 0: gamma0 h_act exceeds delta_cert - delta_total, so the Corollary margin fails and the Corollary refuses",
        g0h2 > marg2 and (not o2["conditions"]["margin"]) and (not o2["valid"]), gamma0_h_act=g0h2, margin=marg2,
        failing=[kk for kk, v in o2["conditions"].items() if not v])
L.check("stochastic runs use r_W2 = 0 (no distributional drift of the adversary)", float(st["r_W2"]) == 0, r_W2=st["r_W2"])
L.data["R2_sigma0"] = dict(gamma0_h_act=g0h2, margin=marg2, valid=o2["valid"], r_W2=float(st["r_W2"]),
                           failing=[kk for kk, v in o2["conditions"].items() if not v])
write_table("exp7_stochastic", ["noise", "sigma", "$h(x_0)$", "Corollary", "crude", "T_TTL", "$T_{\\mathrm{TTL}}/T$", "exits",
                                "CP $%s\\%%$ interval" % num(100 * conf).strip("$")], rows,
            caption=r"Corollary~\ref{cor:stochastic_w2} (R1, $T=%s$, $\alpha=%s$, $r_{W_{\scriptstyle 2}}=%s$, worst-case adversary, $z(0)$ on the "
                    r"level $%s$ of $\tilde{\mathcal{C}}_r$): Monte Carlo exit frequency vs.\ the bound $\alpha$ with the exact "
                    r"Clopper--Pearson interval; crude: bound with $\sigma_\perp:=0$. On R2 the Corollary refuses already at "
                    r"$\sigma=0$: $\gamma_0h_{\mathrm{act}}=%s$ exceeds the margin $\delta_{\mathrm{cert}}-\delta_{\mathrm{total}}=%s$."
                    % (num(T), num(alpha), num(float(st["r_W2"])).strip("$"), num(float(st["z0_level"])).strip("$"), num(g0h2), num(marg2)),
            wide=True)
write_table("exp7_corollary_det", ["bound", "sigma_r", "sigma_perp2", "Delta_sigma", "eb_sigma", "delta_total_sigma", "T_perp"], crow_d,
            caption=r"Corollary~\ref{cor:stochastic_w2}, diffusion-enlarged deterministic constants (full vs.\ crude bound); "
                    r"the tail envelope must stay below $\sqrt{E_{\mathrm{safe}}}=%s$." % num(np.sqrt(c["E_safe"])), wide=True)
write_table("exp7_corollary_stoch", ["bound", "r_BM", "delta_Ito", "sigma_h", "Delta_delta_pw", "delta_stoch", "c_int"], crow_s,
            caption=r"Corollary~\ref{cor:stochastic_w2}, stochastic constants of \eqref{eq:delta_stoch} and the Kushner bound "
                    r"(full vs.\ crude bound).", wide=True)

fig, ax = pl.figure(2, 1, 3.58, 2.3)
for col, key in enumerate(("sigma_star", "noise_tail")):
    r, ncfg = runs[key]
    a = ax[0, col]
    for j in range(min(n, int(cfg["sim"]["n_seeds"]) * int(cfg["sim"]["n_seeds"]))):
        a.plot(r["t"], r["h"][:, j], lw=0.3, rasterized=True)
    if key == "sigma_star":
        a.plot(det["t"], h_det, "k--", lw=0.9, label="noise-free")
        a.legend(loc="lower right")
        a.set_title(r"profile noise, $\sigma^\star=%g$: raw $h$" % ncfg["sigma"])
    else:
        a.set_title(r"(f) tail-white, $\sigma=%g$: raw $h$" % ncfg["sigma"])
    pl.boundary(a)
    a.set_xlabel("$t$"); a.set_ylabel("$h(x(t))$")
pl.save(fig, FIG, "exp7_stochastic")
L.finish()
