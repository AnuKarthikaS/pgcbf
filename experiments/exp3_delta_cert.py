"""Experiment 3 - the certified CBF margin delta_cert (eq:teacher_cbf).

(i) Class (H) (R2, R1), Lemma lem:delta_cert_classH (c_M = lambda_min(M,P) - 2 s_nl, two-branch delta_cert):
    the sound bound used by the certificate is <= the brute-force
    infimum of the exact bracket L_f h~ + L_g h~ T + gamma0 h~ over >= 2e5 samples of C_act
    (both shell boundaries included).
(ii) Lemma lem:delta_cert_closed (class (D): P = I, g_r = I, m = K): the closed form
    eq:delta_cert_closed equals the brute-force infimum of its affine lower bound over rho, is
    <= the brute-force infimum of the exact bracket, and selects branch 1 iff c >= gamma0."""
import math
import numpy as np
from _common import Ledger, load_config, write_table, FIG
from pgcbf.certificates import build, certificate, bracket, sample_shell, classD_closed_form, two_branch
from pgcbf.controller import Teacher, random_orthogonal, project_U
from pgcbf import plotting as pl

L = Ledger("exp3_delta_cert")
rows = []
fig, ax = pl.figure(1, 2, 3.5, 3.0)
fig.get_layout_engine().set(h_pad=6 / 72)
for (name, cfgname), col in zip((("R2", "default"), ("R1", "r1_passive")), ("C0", "C1")):
    cfg = load_config(cfgname)
    sysd = build(cfg, train=True)
    c = certificate(cfg, sysd)
    rng = np.random.default_rng(int(cfg["seeds"]["mc"]))
    n = int(cfg["experiments"]["shell_samples"])
    Z = sample_shell(sysd["rom"], c["r_min2"], c["h_act"], c["Delta_h_dip"], n, rng)
    chunk = int(cfg["experiments"]["crown_samples"]) // int(cfg["sim"]["n_seeds"])
    br = np.concatenate([bracket(sysd["rom"], sysd["policy"], Z[:, i:i + chunk], c["E_safe"], c["eps_basis_P"],
                                 c["Delta_h_dip"], c["gamma0"]) for i in range(0, n, chunk)])
    L.check(f"{name}: sound two-branch delta_cert <= brute-force inf over C_act", c["delta_cert"] <= br.min(),
            delta_cert=c["delta_cert"], brute_min=br.min(), branch=c["delta_cert_branch"], samples=n)
    L.check(f"{name}: (H6) delta_total <= delta_cert <= brute-force inf", c["delta_total"] <= c["delta_cert"] <= br.min(),
            delta_total=c["delta_total"])
    L.check(f"{name}: Lemma delta_cert_classH: c_M = lambda_min(M,P) - 2 s_nl",
            abs(c["c_M"] - (c["lambda_min_MP"] - 2 * c["s_nl"])) <= float(cfg["tolerances"]["roundoff_rel"]) * abs(c["c_M"]),
            c_M=c["c_M"], lambda_min_MP=c["lambda_min_MP"], s_nl=c["s_nl"])
    rows.append([name, "Lemma~\\ref{lem:delta_cert_classH}", c["c_M"], c["lambda_min_MP"], c["s_nl"], c["gamma0"], c["delta_cert_branch"],
                 c["delta_cert"], br.min(), c["delta_total"]])
    L.data[name] = {"c_M": c["c_M"], "lambda_min_MP": c["lambda_min_MP"], "s_nl": c["s_nl"], "delta_cert": c["delta_cert"],
                    "brute": br.min(), "ratio_brute_over_bound": br.min() / c["delta_cert"]}
    rho = sysd["rom"].Pnorm(Z) ** 2
    ax[0, 0].plot(rho, br, ",", color=col, alpha=0.3, rasterized=True)
    ax[0, 0].plot([], [], "s", color=col, label=f"{name}: samples")
    ax[0, 0].axhline(c["delta_cert"], color=col, ls="--", lw=0.8, label=f"{name}: " + r"$\delta_{\mathrm{cert}}$")
    ax[0, 0].axhline(c["delta_total"], color=col, ls=":", lw=0.8, label=f"{name}: " + r"$\delta_{\mathrm{total}}$")
ax[0, 0].set_xlabel(r"$\|z\|_P^2$ on $\mathcal{C}_{r,\mathrm{act}}$"); ax[0, 0].set_ylabel("exact bracket")
ax[0, 0].legend(loc="upper center", bbox_to_anchor=(1 / 2, -0.3), ncol=2, columnspacing=1)
ax[0, 0].set_title("Class (H): exact bracket samples vs. certified bound")

# ---------------- class (D) instance of Lemma lem:delta_cert_closed
cfg = load_config("default")
cd = cfg["experiments"]["classD"]
K = int(cd["K"])
rng = np.random.default_rng(int(cd["seed"]))
R0 = np.diag(np.array(cd["R0_diag"], float))
Ar = rng.standard_normal((K, K))
J0 = (Ar - Ar.T) / 2
W = random_orthogonal(K, rng)
L_F = float(cd["L_F"])
teacher = Teacher(K, int(cd["H"]), int(cd["L"]), float(cd["theta"]), float(cd["L_T"]), int(cd["seed"]))
r2, h_act, dh = float(cd["r_min2"]), float(cd["h_act"]), float(cd["dh"])
radius = np.sqrt(r2 + h_act + dh)
u_max = teacher.L_T * radius  # (H-sat) with equality on the shell: no clipping
c_closed = 2 * (np.linalg.eigvalsh(R0).min() - L_F + teacher.eta)
n = int(cfg["experiments"]["shell_samples"])
G = rng.standard_normal((K, n))
U = G / np.linalg.norm(G, axis=0)
lev = r2 + (h_act + dh) * rng.random(n)
lev[: n // 4] = r2
lev[n // 4: n // 2] = r2 + h_act + dh
Z = U * np.sqrt(lev)
f = (J0 - R0) @ Z + L_F * np.sin(W @ Z)
T = project_U(teacher.forward(Z), u_max)
htil = r2 + h_act - lev
rho_grid = np.linspace(r2, r2 + h_act + dh, int(cfg["sim"]["n_seeds"]) ** 2 + 1)
for g0, col in zip(cd["gamma0_values"], ("C2", "C3")):
    g0 = float(g0)
    exact = -2 * np.sum(Z * f, axis=0) - 2 * np.sum(Z * T, axis=0) + g0 * htil
    cf = classD_closed_form(c_closed, r2, h_act, dh, g0)
    affine = (c_closed - g0) * rho_grid + g0 * (r2 + h_act)
    tol = float(cfg["tolerances"]["roundoff_rel"]) * abs(cf["delta_cert"])
    L.check(f"class D gamma0={g0}: closed form = brute-force inf of affine bound over rho", abs(affine.min() - cf["delta_cert"]) <= tol,
            closed=cf["delta_cert"], brute_affine=affine.min())
    L.check(f"class D gamma0={g0}: closed form <= brute-force inf of exact bracket", cf["delta_cert"] <= exact.min(),
            closed=cf["delta_cert"], brute_exact=exact.min())
    L.check(f"class D gamma0={g0}: branch 1 iff c >= gamma0", (cf["branch"] == 1) == (c_closed >= g0), c=c_closed, branch=cf["branch"])
    rows.append(["class (D)", "\\eqref{eq:delta_cert_closed}", c_closed, None, None, g0, cf["branch"], cf["delta_cert"], exact.min(),
                 "N/A (standalone)"])
    ax[1, 0].plot(lev, exact, ",", color=col, alpha=0.3, rasterized=True)
    g0s = r"$\gamma_0=%g$" % g0
    ax[1, 0].plot([], [], "s", color=col, label="exact, " + g0s)
    ax[1, 0].plot(rho_grid, affine, "-", color="k" if col == "C2" else "C1", lw=0.8, label="affine, " + g0s)
    ax[1, 0].axhline(cf["delta_cert"], color="k" if col == "C2" else "C1", ls="--", lw=0.8, label="closed form, " + g0s)
ax[1, 0].set_xlabel(r"$\rho=\|z\|^2$"); ax[1, 0].set_ylabel("bracket")
ax[1, 0].legend(loc="upper center", bbox_to_anchor=(1 / 2, -0.3), ncol=2, columnspacing=1)
ax[1, 0].set_title("Class (D): two-branch closed form vs. exact bracket")
pl.save(fig, FIG, "exp3_delta_cert")
write_table("exp3_delta_cert", ["case", "bound", "c_M", "lambda_min_MP", "s_nl", "gamma0", "branch", "$\\delta_{\\mathrm{cert}}$ (bound)",
                                "brute-force inf", "delta_total"], rows,
            caption="Certified margin $\\delta_{\\mathrm{cert}}$: class-(H) bound of Lemma~\\ref{lem:delta_cert_classH} "
                    "($c_M=\\lambda_{\\min}(M,P)-2s_{\\mathrm{nl}}$, two-branch infimum) and the class-(D) closed form "
                    "\\eqref{eq:delta_cert_closed} vs.\\ the brute-force infimum of the exact bracket over "
                    "$\\mathcal{C}_{r,\\mathrm{act}}$. For class (D) the column $c_M$ holds the closed-form slope $c$.")
L.finish()
