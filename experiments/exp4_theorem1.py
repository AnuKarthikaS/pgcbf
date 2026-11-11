"""Experiment 4 - Theorem thm:main_cbf, positive case (full method).

Full-order sampled-data closed loop with exact multi-rate latch timing, quantized Path A
(LUT) / Path B (QNN), worst-case (HJI-aligned) adversary w* of norm eps_adv, from (H7)
initial states with the maximal admissible tail, in two initial sets per regime:
  boundary: z(0) on the boundary of C~_r (h~_r(z(0)) = 0);
  interior (S6): z(0) on the level set h_r(z) = E_safe - ||z||_P^2 = interior_level * h_act, strictly inside C~_r,
                with tail fraction interior_tail_fraction (0: the tail is excited only by the closed loop, so that the
                minimum of h is attained along the trajectory rather than at t = 0 by a decaying initial tail).
Asserts, for every seed and every integrator node: (a) ||x_perp||_Pfe <= eps_basis_P; (b) h~_r(z(t_k)) >= 0 and
h~_r(z(t)) >= -Delta h_dip; (c) h(x(t)) >= h~_r(z(t)) + Delta h_dip >= 0.  Reports min_t h(x_t) and its time
(t = 0 for both starts: the inter-sample dip floor is not exercised dynamically; budget utilisation along these
trajectories is measured by Experiment 4b, exp4b_budget.py)."""
import numpy as np
from _common import Ledger, load_config, write_table, FIG
from pgcbf.texnames import num
from pgcbf.certificates import build, certificate, refuse_if_invalid
from pgcbf.simulate import Simulator, initial_conditions
from pgcbf import plotting as pl

L = Ledger("exp4_theorem1")
rows = []
fig, ax = pl.figure(2, 3, 3.5, 1.9)
for col, (name, cfgname) in enumerate((("R2", "default"), ("R1", "r1_passive"))):
    cfg = load_config(cfgname)
    sysd = build(cfg, train=True)
    c = certificate(cfg, sysd)
    refuse_if_invalid(c)
    lat, rom, pol = sysd["lat"], sysd["rom"], sysd["policy"]
    sc = cfg["sim"]
    eb, dh = c["eps_basis_P"], c["Delta_h_dip"]
    lev_int = float(cfg["experiments"]["interior_level"]) * c["h_act"]
    L.check(f"{name}: interior level h_r = {cfg['experiments']['interior_level']} h_act lies strictly inside C~_r",
            lev_int > eb**2 + dh, h_r0=lev_int, eps_b2_plus_dh=eb**2 + dh)
    sim = Simulator(lat, rom, pol, cfg, adversary="worst")
    n = int(round(float(sc["horizon"]) / sim.Ts))
    tol = float(cfg["tolerances"]["roundoff_rel"]) * c["E_safe"]
    for start, z2, tf in (("boundary", None, float(sc["ic_tail_fraction"])),
                          ("interior", c["E_safe"] - lev_int, float(cfg["experiments"]["interior_tail_fraction"]))):
        X0 = initial_conditions(lat, rom, c, int(sc["n_seeds"]), np.random.default_rng(int(cfg["seeds"]["ic"])),
                                float(sc["ic_level_fraction"]), tf, int(sc["ic_n_tail_modes"]), z_norm2=z2)
        r = sim.run(X0, n, record="substep")
        htil = r["hPhi"] - eb**2 - dh
        k_idx = np.arange(0, r["h"].shape[0], sum(sim.sub))
        tag = f"{name} {start}"
        L.check(f"{tag}: (a) tail ||x_perp(t)||_Pfe <= eps_basis_P (all seeds, all nodes)", np.all(r["tail"] <= eb * (1 + tol)),
                max_tail_over_eps_b=r["tail"].max() / eb)
        L.check(f"{tag}: (b) h~_r(z(t_k)) >= 0 at every sampling instant", np.all(htil[k_idx] >= -tol), min_htil_tk=htil[k_idx].min())
        L.check(f"{tag}: (b) h~_r(z(t)) >= -Delta h_dip at every node", np.all(htil >= -dh - tol), min_htil=htil.min(), Delta_h_dip=dh)
        L.check(f"{tag}: (c) h(x(t)) >= h~_r(z(t)) + Delta h_dip", np.all(r["h"] >= htil + dh - tol),
                min_gap=(r["h"] - htil - dh).min())
        L.check(f"{tag}: (c) h(x(t)) >= 0 for all t, all seeds (full-order safety)", np.all(r["h"] >= 0), min_h=r["h"].min(),
                delta_total=c["delta_total"], seeds=X0.shape[1])
        L.check(f"{tag}: z(t_k) stays in Z_ODD", not np.any(r["out_of_odd"]))
        i_min, j_min = np.unravel_index(np.argmin(r["h"]), r["h"].shape)
        t_min, h_min = float(r["t"][i_min]), float(r["h"][i_min, j_min])
        if start == "interior":
            # observed (S6): the certified closed loop raises h from every start in C_act, so the minimum is h(x_0)
            L.check(f"{tag}: min_t h(x_t) = min_j h(x_j(0)) (attained at t = 0; h rises along every closed-loop trajectory)",
                    t_min == 0 and h_min >= r["h"][0].min() - tol, t_min=t_min, min_h=h_min)
        else:
            L.check(f"{tag}: min_t h = Delta h_dip at t = 0 (initial state on the certified floor)",
                    t_min == 0 and abs(h_min - dh) <= tol + float(cfg["tolerances"]["roundoff_rel"]) * dh * lat.N, min_h=h_min, Delta_h_dip=dh)
        rows.append([name, start, X0.shape[1], h_min, t_min, htil.min(), r["tail"].max() / eb, c["delta_cert"] - c["delta_total"]])
        L.data[tag] = {"min_h": h_min, "t_min": t_min, "h0_min": float(r["h"][0].min()),
                       "min_h_per_seed": r["h"].min(axis=0), "min_htil": htil.min(), "max_tail": r["tail"].max()}
        a = ax[0 if start == "boundary" else 1, col]
        for j in range(X0.shape[1]):
            a.plot(r["t"], r["h"][:, j], lw=0.5)
        pl.boundary(a)
        a.axhline(dh, color="r", ls=":", lw=0.8, label=r"$\Delta h_{\mathrm{dip}}$")
        a.set_title(f"{name}, {start} start: raw $h$"); a.set_xlabel(r"$t$"); a.set_ylabel(r"$h(x(t))$")
        if start == "boundary":
            a = ax[2, col]
            for j in range(X0.shape[1]):
                a.plot(r["t"], r["tail"][:, j] / eb, lw=0.5)
            a.axhline(1, color="k", ls="--", lw=0.8, label=r"$\varepsilon_{\mathrm{basis},P}$")
            a.set_title(f"{name}: tail ratio, boundary start"); a.set_xlabel(r"$t$")
            a.set_ylabel(r"tail ratio $\|x_\perp\|/\varepsilon_{\mathrm{basis},P}$"); a.legend(loc="center right")
# one shared legend below the panels for the h rows (an in-axes legend would cover the fast initial rise)
fig.legend(*ax[0, 0].get_legend_handles_labels(), loc="outside lower center", ncol=2)
pl.save(fig, FIG, "exp4_theorem1")
write_table("exp4_theorem1", ["config", "start", "seeds", "$\\min_t h(x_t)$", "time of min",
                              "$\\min \\tilde h_r$", "max tail\\,/\\,$\\varepsilon_{\\mathrm{basis},P}$", "H6_margin"], rows,
            caption="Theorem~\\ref{thm:main_cbf}: measured minima over all seeds and integrator nodes for the boundary start "
                    "($\\tilde h_r(z(0))=0$) and the interior start ($h_r(z(0))=%s\\,h_{\\mathrm{act}}$), worst-case adversary."
                    % num(float(cfg["experiments"]["interior_level"])).strip("$"))
L.finish()
