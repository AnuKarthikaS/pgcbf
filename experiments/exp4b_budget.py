"""Experiment 4b - realised budget utilisation along the certified trajectories of Experiment 4.

Derivation: supplement, Experiment 4b (Step 2 of the proof of Theorem thm:main_cbf, eq:ode_diff_ineq).  Same closed
loops as Experiment 4 (both regimes, both starts, all seeds, worst-case adversary w*, same initial conditions).  At the
left end point of every integrator step (latch phase = simulator segment index, right-continuous input) the exact
identity
    r := d/dt h~_r(z) + gamma0 h~_r(z) = Theta(z_k) + H + Q + A + C
is evaluated term by term (C tail coupling, A adversary, Q quantization/lag, H hold shift) and each term is compared
with the bound the proof uses for it:
    (i)   |grad h~_r(z)^T Delta_perp|            <= L_h L_perp,P eps_basis,P
    (ii)  |L_dr h~_r(z) w|                       <= d_h eps_adv
    (iii) |L_gr h~_r(z_k) (u - T(z_k))|          <= g_h ebar_u(phase)
    (iv)  |H|                                    <= L_hold s
    (v)   r - RHS(s) >= 0, RHS(s) = delta_total - L_h L_perp eps_b - d_h eps_adv - L_hold s - g_h e_u(t)
          on periods with z_k in C_act; periods with z_k in B_int use the interior (Lipschitz) branch of the proof.
"""
import numpy as np
from _common import Ledger, load_config, write_table
from pgcbf.certificates import build, certificate, refuse_if_invalid, bracket
from pgcbf.simulate import Simulator, initial_conditions

L = Ledger("exp4b_budget")
TERMS = ("tail", "adv", "quant", "hold")
FIELDS = TERMS + ("phase", "s", "active", "slack", "slackbar", "slack_int", "rho_int", "f_grad", "f_lip", "f_env",
                  "f_cos", "f_eu", "f_gk", "f_mz", "f_lipH", "f_w", "id_res", "dz_res")


class Recorder:
    """Read-only node hook for Simulator.run (see the derivation note, Sections 2-5)."""

    def __init__(self, lat, rom, pol, sim, c, u_max):
        self.lat, self.rom, self.pol, self.sim, self.c = lat, rom, pol, sim, c
        self.u_max = u_max
        self.eb, self.dh, self.E, self.g0 = c["eps_basis_P"], c["Delta_h_dip"], c["E_safe"], c["gamma0"]
        self.B_tail = c["L_h"] * c["L_perp"] * c["eps_basis_P"]
        self.B_adv = c["d_h"] * c["eps_adv"]
        self.L_sum = c["L_Lfh"] + c["L_Lgh"] * u_max + c["gamma0"] * c["L_h"]
        self.k_cur = -1
        self.rows = {f: [] for f in FIELDS}
        self.theta_margin = []

    def htil(self, Z):
        return self.rom.hPhi(Z, self.E) - self.eb**2 - self.dh

    def lie(self, Z):
        rom = self.rom
        PZ = rom.P @ Z
        return PZ, -2 * np.sum(PZ * rom.f_r(Z), axis=0), -2 * (rom.g_r.T @ PZ)

    def ebar(self, i, k):
        c = self.c
        if i == 0:
            return c["eps_CROWN"] + (c["du_lag"] if k >= 1 else 0 * c["du_lag"])
        return c["eps_A"] if i == 1 else c["eps_CROWN"]

    def __call__(self, k, i, s, Y, U, Zk, W):
        lat, rom, c, g0 = self.lat, self.rom, self.c, self.g0
        if k != self.k_cur:  # new sampling instant: quantities of z_k, once per period
            self.k_cur = k
            _, self.Lfh_k, self.Lgh_k = self.lie(Zk)
            self.h_k = self.htil(Zk)
            self.T_k = self.pol.teacher(Zk)
            self.Theta_k = bracket(rom, self.pol, Zk, self.E, self.eb, self.dh, g0)
            self.Zk = Zk
            self.act_k = self.h_k <= c["h_act"]
            self.theta_margin.append(np.where(self.act_k, self.Theta_k - c["delta_cert"], np.inf))
        X = self.sim.T @ Y
        Z = self.sim.Zmap @ Y  # the map the simulator uses for z_k, so H(t_k) = 0 and slack_int(t_k) = 0 identically
        Xp = self.sim.T @ (self.sim.Pperp @ Y)  # modal tail selection: no roundoff leak of Phi z into x_perp
        PhiZ = rom.Phi @ Z
        Xdot = lat.Af @ X + lat.F_nl(X) + lat.G @ U + lat.D @ W
        zdot = rom.Zmap @ Xdot
        Dp = rom.Zmap @ (lat.Af @ Xp + lat.F_nl(PhiZ + Xp) - lat.F_nl(PhiZ))  # eq:Lperp_def
        fr = rom.f_r(Z)
        dz_res = np.linalg.norm(zdot - fr - rom.g_r @ U - rom.d_r @ W - Dp, axis=0) / (
            np.linalg.norm(zdot, axis=0) + np.linalg.norm(fr, axis=0) + np.linalg.norm(Dp, axis=0))
        PZ, Lfh, Lgh = self.lie(Z)
        ht = self.htil(Z)
        grad = -2 * PZ
        C = np.sum(grad * Dp, axis=0)
        A = -2 * np.sum((rom.d_r.T @ PZ) * W, axis=0)
        Q = np.sum(self.Lgh_k * (U - self.T_k), axis=0)
        H = (Lfh - self.Lfh_k) + np.sum((Lgh - self.Lgh_k) * U, axis=0) + g0 * (ht - self.h_k)
        r = np.sum(grad * zdot, axis=0) + g0 * ht
        parts = np.abs(self.Theta_k) + np.abs(H) + np.abs(Q) + np.abs(A) + np.abs(C) + np.abs(r)
        id_res = np.abs(r - (self.Theta_k + H + Q + A + C)) / parts
        e_u = np.linalg.norm(U - self.T_k, axis=0)
        eb_u = self.ebar(i, k)
        base = c["delta_total"] - self.B_tail - self.B_adv - c["L_hold"] * s
        slack = r - (base - c["g_h"] * e_u)
        slackbar = r - (base - c["g_h"] * eb_u)
        dzn = np.linalg.norm(Z - self.Zk, axis=0)
        tailn = np.sqrt(np.sum(Xp * (lat.Pfe @ Xp), axis=0))
        gn, dpn = np.linalg.norm(grad, axis=0), np.linalg.norm(Dp, axis=0)
        nan = np.full(Z.shape[1], np.nan)
        pos = s > 0
        R = self.rows
        R["tail"].append(np.abs(C) / self.B_tail)
        R["adv"].append(np.abs(A) / self.B_adv)
        R["quant"].append(np.abs(Q) / (c["g_h"] * eb_u))
        R["hold"].append(np.abs(H) / (c["L_hold"] * s) if pos else nan)
        R["phase"].append(np.full(Z.shape[1], i))
        R["s"].append(np.full(Z.shape[1], s))
        R["active"].append(self.act_k.copy())
        R["slack"].append(np.where(self.act_k, slack, np.nan))
        R["slackbar"].append(np.where(self.act_k, slackbar, np.nan))
        R["slack_int"].append(np.where(~self.act_k, ht - (self.h_k - c["L_h"] * c["M_z"] * s), np.nan))
        R["rho_int"].append(np.where(~self.act_k, np.maximum(0 * ht, self.h_k - ht) / (c["L_h"] * c["M_z"] * s), np.nan)
                            if pos else nan)
        R["f_grad"].append(gn / c["L_h"])
        R["f_lip"].append(np.where(tailn > 0, dpn / (c["L_perp"] * np.where(tailn > 0, tailn, 1)), np.nan))
        R["f_env"].append(tailn / self.eb)
        R["f_cos"].append(np.where(gn * dpn > 0, np.abs(C) / np.where(gn * dpn > 0, gn * dpn, 1), np.nan))
        R["f_eu"].append(e_u / eb_u)
        R["f_gk"].append(np.linalg.norm(self.Lgh_k, axis=0) / c["g_h"])
        R["f_mz"].append(dzn / (c["M_z"] * s) if pos else nan)
        R["f_lipH"].append(np.where(dzn > 0, np.abs(H) / (self.L_sum * np.where(dzn > 0, dzn, 1)), np.nan) if pos else nan)
        R["f_w"].append(np.linalg.norm(W, axis=0) / c["eps_adv"])
        R["id_res"].append(id_res)
        R["dz_res"].append(dz_res)
        # decomposition of the slack at this node (kept only for the running minimum)
        dec = {"teacher": self.Theta_k - c["delta_total"], "tail": C + self.B_tail, "adv": A + self.B_adv,
               "hold": H + c["L_hold"] * s, "quant": Q + c["g_h"] * e_u}
        sl = np.where(self.act_k, slack, np.inf)
        j = int(np.argmin(sl))
        if sl[j] < getattr(self, "min_slack", np.inf):
            self.min_slack = float(sl[j])
            self.min_dec = {kk: float(v[j]) for kk, v in dec.items()}
            self.min_where = {"k": int(k), "phase": int(i), "s": float(s), "seed": j}

    def arrays(self):
        return {f: np.concatenate(v) for f, v in self.rows.items()}


def stats(D, mask=None):
    out = {}
    for t in TERMS + ("f_grad", "f_lip", "f_env", "f_cos", "f_eu", "f_gk", "f_mz", "f_lipH", "f_w"):
        x = D[t] if mask is None else D[t][mask]
        x = x[np.isfinite(x)]
        out[t] = {"max": float(x.max()), "median": float(np.median(x)), "n": int(x.size)} if x.size else None
    return out


rows, csv_rows = [], []
for name, cfgname in ((("R2", "default"), ("R1", "r1_passive"))):
    cfg = load_config(cfgname)
    sysd = build(cfg, train=True)
    c = certificate(cfg, sysd)
    refuse_if_invalid(c)
    lat, rom, pol = sysd["lat"], sysd["rom"], sysd["policy"]
    sc = cfg["sim"]
    u_max = float(cfg["controller"]["u_max"])
    rr = float(cfg["tolerances"]["roundoff_rel"])
    lev_int = float(cfg["experiments"]["interior_level"]) * c["h_act"]
    sim = Simulator(lat, rom, pol, cfg, adversary="worst")
    n = int(round(float(sc["horizon"]) / sim.Ts))
    tol = rr * (c["delta_cert"] + c["L_h"] * c["L_perp"] * c["eps_basis_P"] + c["d_h"] * c["eps_adv"]
                + c["L_hold"] * sim.Ts + c["g_h"] * (c["eps_A"] + c["eps_CROWN"] + c["du_lag"]))
    rtol = rr * lat.N
    pooled = []
    for start, z2, tf in (("boundary", None, float(sc["ic_tail_fraction"])),
                          ("interior", c["E_safe"] - lev_int, float(cfg["experiments"]["interior_tail_fraction"]))):
        X0 = initial_conditions(lat, rom, c, int(sc["n_seeds"]), np.random.default_rng(int(cfg["seeds"]["ic"])),
                                float(sc["ic_level_fraction"]), tf, int(sc["ic_n_tail_modes"]), z_norm2=z2)
        rec = Recorder(lat, rom, pol, sim, c, u_max)
        sim.run(X0, n, record="period", node_hook=rec)
        D = rec.arrays()
        pooled.append(D)
        tag = f"{name} {start}"
        act = D["active"]
        n_nodes = int(D["tail"].size)
        n_act, n_int = int(act.sum()), int((~act).sum())
        thm = np.concatenate(rec.theta_margin)
        L.check(f"{tag}: node set = every integrator step of every period and seed", n_nodes == n * sum(sim.sub) * X0.shape[1],
                nodes=n_nodes, active_nodes=n_act, interior_nodes=n_int)
        L.check(f"{tag}: decomposition zdot = f_r + g_r u + d_r w + Delta_perp closes to roundoff (every node)",
                np.all(D["dz_res"] <= rtol), max_rel_residual=D["dz_res"].max())
        L.check(f"{tag}: identity r = Theta(z_k) + H + Q + A + C closes to roundoff (every node)",
                np.all(D["id_res"] <= rtol), max_rel_residual=D["id_res"].max())
        for t, lab in (("tail", "(i) tail coupling"), ("adv", "(ii) adversary"), ("quant", "(iii) quantization/lag"),
                       ("hold", "(iv) hold shift")):
            x = D[t][np.isfinite(D[t])]
            L.check(f"{tag}: {lab} realised/bound <= 1 at every node", np.all(x <= 1 + rtol), max_ratio=x.max(),
                    median_ratio=np.median(x))
        fac = {f: D[f][np.isfinite(D[f])].max() for f in ("f_grad", "f_lip", "f_eu", "f_gk", "f_w")}
        L.check(f"{tag}: factor bounds ||2Pz|| <= L_h, ||Delta_perp|| <= L_perp ||x_perp||, e_u <= ebar_u, "
                "||L_gr h~_r(z_k)|| <= g_h, ||w|| <= eps_adv at every node", all(v <= 1 + rtol for v in fac.values()), **fac)
        mz = D["f_mz"][np.isfinite(D["f_mz"])]
        L.check(f"{tag}: ||z - z_k|| <= M_z s at every node", np.all(mz <= 1 + rtol), max_ratio=mz.max())
        sl = D["slack"][act]
        L.check(f"{tag}: (v) soundness: realised rate - RHS(s) of eq:ode_diff_ineq >= 0 at every active-shell node",
                np.all(sl >= -tol), min_slack=sl.min() if sl.size else np.nan, tol=tol, nodes=int(sl.size))
        L.check(f"{tag}: Teacher condition Theta(z_k) >= delta_cert at every active-shell sample",
                np.all(thm[np.isfinite(thm)] >= -tol), min_margin=thm[np.isfinite(thm)].min() if np.isfinite(thm).any() else np.nan)
        si = D["slack_int"][~act]
        L.check(f"{tag}: interior branch h~_r(z(t)) >= h~_r(z_k) - L_h M_z s at every B_int node", np.all(si >= -tol),
                min_slack=si.min() if si.size else np.nan, nodes=int(si.size))
        L.data[tag] = {"nodes": n_nodes, "active_nodes": n_act, "interior_nodes": n_int, "stats": stats(D),
                       "stats_active": stats(D, act),
                       "quant_phase": {str(p): stats(D, D["phase"] == p)["quant"] for p in range(3)},
                       "min_slack": float(sl.min()) if sl.size else None,
                       "median_slack": float(np.median(sl)) if sl.size else None,
                       "min_slackbar": float(D["slackbar"][act].min()) if sl.size else None,
                       "min_slack_decomposition": getattr(rec, "min_dec", None), "min_slack_where": getattr(rec, "min_where", None),
                       "min_teacher_margin": float(thm[np.isfinite(thm)].min()) if np.isfinite(thm).any() else None,
                       "min_slack_int": float(si.min()) if si.size else None,
                       "max_rho_int": float(np.nanmax(D["rho_int"][~act])) if si.size else None}
    P = {f: np.concatenate([d[f] for d in pooled]) for f in FIELDS}
    st = stats(P)
    act = P["active"]
    sl = P["slack"][act]
    sa = stats(P, act)
    qp = {p: stats(P, P["phase"] == p)["quant"] for p in range(3)}
    qa = {p: stats(P, act & (P["phase"] == p))["quant"] for p in range(3)}
    L.data[name] = {"stats": st, "stats_active": sa, "quant_phase_active": {str(p): v for p, v in qa.items()}, "quant_phase": {str(p): v for p, v in qp.items()}, "min_slack": float(sl.min()),
                    "median_slack": float(np.median(sl)), "min_slackbar": float(P["slackbar"][act].min()),
                    "delta_total": c["delta_total"], "delta_cert": c["delta_cert"], "H6_margin": c["H6_margin"],
                    "nodes": int(P["tail"].size), "active_nodes": int(act.sum()), "interior_nodes": int((~act).sum())}
    for t, lab in (("tail", "(i) tail coupling"), ("adv", "(ii) adversary"), ("quant", "(iii) quantization/lag"),
                   ("hold", "(iv) hold shift")):
        rows.append([name, lab, st[t]["max"], st[t]["median"], sa[t]["max"], sa[t]["median"]])
    for p, lab in ((0, "(iii) phase $[0,\\tau_A)$"), (1, "(iii) phase $[\\tau_A,\\tau_B)$"), (2, "(iii) phase $[\\tau_B,T_s)$")):
        rows.append([name, lab, qp[p]["max"], qp[p]["median"], qa[p]["max"], qa[p]["median"]])
    rows.append([name, "(v) slack$\\,/\\,\\delta_{\\mathrm{total}}$: min, median", np.nan, np.nan,
                 float(sl.min()) / c["delta_total"], float(np.median(sl)) / c["delta_total"]])
write_table("exp4b_budget", ["config", "term", "max, all nodes", "median, all nodes", "max, $\\mathcal{C}_{r,\\mathrm{act}}$",
                            "median, $\\mathcal{C}_{r,\\mathrm{act}}$"], rows,
            caption="Experiment 4b: realised fraction of each budget term of \\eqref{eq:ode_diff_ineq} along the certified "
                    "trajectories of Experiment 4 (both starts, all seeds, every integrator node, worst-case adversary); "
                    "row (v): minimum and median of the realised rate minus the right-hand side of \\eqref{eq:ode_diff_ineq}, "
                    "in units of $\\delta_{\\mathrm{total}}$. Columns $\\mathcal{C}_{r,\\mathrm{act}}$: nodes of periods whose "
                    "sample $z_k$ lies in the active shell, where the budget enters the proof.", wide=True)
L.finish()
