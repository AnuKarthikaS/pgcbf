"""Typesetting dictionary (V1): every quantity key, hypothesis key and equation label used in the generated
tables and figures maps to LaTeX.  The generated .tex never contains a raw Python identifier, boolean or nan."""
from __future__ import annotations

import math

NAMES = {
    # dimensions / plant
    "N": r"$N$", "K": r"$K$", "m": r"$m$", "K_m": r"$K_m$", "K_f": r"$K_f$", "Pe": r"$\mathrm{Pe}$", "c": r"$c$",
    "N_m": r"$N_m$", "nu": r"$\nu$", "alpha_plant": r"$\alpha$ (plant)", "gamma_plant": r"$\gamma$", "kappa": r"$\kappa$",
    "n_adv_full_weight": r"full-weight advection modes",
    "adv_cutoff_weight": r"$w$ (cutoff weight)", "w": r"$w$", "E_safe": r"$E_{\mathrm{safe}}$",
    "lambda_min_Pfe": r"$\lambda_{\min}(P_{\mathrm{full},\epsilon})$", "omega_min": r"$\omega_{\min}$",
    "lambda_min_P": r"$\lambda_{\min}(P)$",
    # Prop. 2 / (H2)
    "mu_2P_Ar0": r"$\mu_{2,P}(A_{r,0})$", "Q_identity_relerr": r"$\|Q-\gamma P\|/\|\gamma P\|$",
    "lambda_K1": r"$\lambda_{K+1}$", "H2a_residual": r"(H2a) residual $\|\Phi^\top P_{\mathrm{full},\epsilon}P_\perp\|$",
    "A_rperp_full_norm": r"$\|A_{r,\perp}\|_2$ (eigenmode invariance)",
    # ODD / tail
    "rho_ODD": r"$\rho_{\mathrm{ODD}}$", "D_ODD": r"$D_{\mathrm{ODD}}$", "L_Fnl": r"$\mathrm{Lip}(F_{\mathrm{nl}})$",
    "M_perp": r"$M_\perp$", "M_perp_lin": r"$M_\perp$ (linear part)", "eps_basis_P": r"$\varepsilon_{\mathrm{basis},P}$",
    "L_perp": r"$L_{\perp,P}$", "L_perp_P": r"$L_{\perp,P}$", "L_perp_lin": r"$L_{\perp,P}$ (advective, linear)",
    "L_perp_nl": r"$L_{\perp,P}$ ($\sin$ nonlinearity)", "A_adv_rperp_norm": r"$\|A_{\mathrm{adv},r,\perp}\|_2$",
    "tail_share": r"tail share $L_hL_{\perp,P}\varepsilon_{\mathrm{basis},P}/\delta_{\mathrm{total}}$",
    "tail_ratio": r"$\|x_\perp\|_{P_{\mathrm{full},\epsilon}}/\varepsilon_{\mathrm{basis},P}$",
    "tail_init_bound": r"$M_\perp/\lambda_{K+1}$",
    # reduced constants
    "L_h": r"$L_h$", "g_h": r"$g_h$", "d_h": r"$d_h$", "M_z": r"$M_z$", "L_Lfh": r"$\mathrm{Lip}(L_fh)$", "L_Lgh": r"$\mathrm{Lip}(L_gh)$",
    "s_q": r"$s_q$", "L_T": r"$L_{\mathcal{T}}$", "L_Tz": r"$L_{\mathcal{T},z}$", "eta_T": r"$\eta_{\mathcal{T}}$",
    "H_sat_lhs": r"(H-sat) left side $L_{\mathcal{T},z}\sqrt{K}D_{\mathrm{ODD}}/2$", "u_max": r"$u_{\max}$",
    "G_port_min_eig": r"$\lambda_{\min}(G_{\mathrm{port}})$",
    # CROWN / LUT
    "eps_CROWN_backward": r"$\varepsilon_{\mathrm{CROWN}}$ (backward)", "E_q": r"$E_q$ (explicit)",
    "eps_CROWN": r"$\varepsilon_{\mathrm{CROWN}}$", "eps_CROWN_portbox": r"$\varepsilon_{\mathrm{CROWN}}$ (port box)",
    "eps_A": r"$\varepsilon_A$", "du_lag": r"$\Delta u_{\mathrm{lag}}$", "b_A": r"$b_A$", "b_B": r"$b_B$",
    "D_lut": r"$D_v$", "D_v": r"$D_v$", "M_v": r"$M_v$", "Delta_qA": r"$\Delta_{q,A}$ (LUT step)",
    "codesign_lo": r"window lower bound", "codesign_hi": r"window upper bound",
    "lut_address_bits": r"LUT address bits $m(b_A+1)$", "lut_storage_bits": r"LUT storage bits (log$_2$)",
    # timing / margins
    "T_s": r"$T_s$", "tau_A": r"$\tau_A$", "tau_B": r"$\tau_B$", "gamma0": r"$\gamma_0$",
    "w1": r"$w_1$", "w2": r"$w_2$", "w3": r"$w_3$", "phi": r"$\phi(\gamma_0T_s)$", "L_hold": r"$L_{\mathrm{hold}}$",
    "Delta_h_dip": r"$\Delta h_{\mathrm{dip}}$",
    "dt_tail": r"$\delta_{\mathrm{total}}$: tail term", "dt_quant": r"$\delta_{\mathrm{total}}$: quantisation term",
    "dt_lag": r"$\delta_{\mathrm{total}}$: lag term", "dt_adv": r"$\delta_{\mathrm{total}}$: adversary term",
    "dt_hold": r"$\delta_{\mathrm{total}}$: hold term", "delta_total": r"$\delta_{\mathrm{total}}$",
    "h_act": r"$h_{\mathrm{act}}$", "r_min": r"$r_{\min}$", "c_M": r"$c_M$", "lambda_min_MP": r"$\lambda_{\min}(M,P)$",
    "s_nl": r"$s_{\mathrm{nl}}$", "delta_cert": r"$\delta_{\mathrm{cert}}$", "delta_cert_branch": r"$\delta_{\mathrm{cert}}$ branch",
    "H6_margin": r"margin $\delta_{\mathrm{cert}}-\delta_{\mathrm{total}}$",
    "null_space_rate": r"null-space rate", "mu_R": r"$\mu_R$", "L_sec": r"$L_{\mathrm{sec}}$",
    "eps_q_bar": r"$\bar\varepsilon_q$", "mu_R_gt_L_sec": r"$\mu_R>L_{\mathrm{sec}}$", "r_ss": r"$r_{\mathrm{ss}}$", "eps_adv": r"$\varepsilon_{\mathrm{adv}}$",
    "open_loop_growth": r"open-loop growth $\max\operatorname{Re}\lambda(A_{\mathrm{lin}})$",
    "L_perp_lin_share": r"linear share of $L_{\perp,P}$", "D_v_sampled": r"$D_v$ (sampled sup)",
    "samples_per_period": r"samples per period $2\pi/(\omega_{\min}T_s)$",
    "hold_share": r"hold share of $\delta_{\mathrm{total}}$", "gamma0_h_act": r"$\gamma_0h_{\mathrm{act}}$ (Cor.~\ref{cor:stochastic_w2})",
    # stochastic
    "sigma": r"$\sigma$", "sigma_r": r"$\sigma_r$", "sigma_perp2": r"$\sigma_\perp^2$", "Delta_sigma": r"$\Delta_\sigma$",
    "eb_sigma": r"$\varepsilon_{\mathrm{basis},\sigma}$", "r_BM": r"$r_{\mathrm{BM}}$", "delta_Ito": r"$\delta_{\mathrm{It\hat o}}$",
    "sigma_h": r"$\sigma_h$", "Delta_delta_pw": r"$\Delta\delta_{\mathrm{pw}}$", "delta_stoch": r"$\delta_{\mathrm{stoch}}$",
    "delta_total_sigma": r"$\delta_{\mathrm{total},\sigma}$", "delta_drift": r"$\delta_{\mathrm{drift}}$",
    "delta_cert_sigma": r"$\delta_{\mathrm{cert},\sigma}$", "c_int": r"$c_{\mathrm{int}}$", "T_perp": r"$T_\perp$",
    "T_TTL": r"$T_{\mathrm{TTL}}$", "h0": r"$h_0$", "alpha": r"$\alpha$",
    # hypotheses
    "H0": "(H0)", "H1": "(H1)", "H2": "(H2)", "H3": "(H3)", "H_sat": "(H-sat)", "H4": "(H4)", "H5": "(H5)",
    "H6": "(H6)", "H7": "(H7)", "Hdz": "(H-dz)", "timing_ok": r"timing $0<\tau_A<\tau_B<T_s$",
    "codesign_window": r"codesign window \eqref{eq:codesign_window}",
    "theorem_valid": r"Thm.~\ref{thm:main_cbf} hypotheses", "prop3_valid": r"Prop.~\ref{prop:deadzone_stability} hypotheses",
}

HYP_DESC = {
    "H0": "ODD containment and $E_{\\mathrm{safe}}$", "H1": "Assumption~\\ref{assump:generic} constants on $\\mathcal{X}_{\\mathrm{ODD}}$",
    "H2": "(H2a)--(H2c), $\\gamma<2\\omega_{\\min}$, $\\epsilon=\\gamma/2$",
    "H3": "Teacher: $H/2\\ge m$, $G_{\\mathrm{port}}\\succ0$, (H-sat)", "H4": "tail tolerance, $h_{\\mathrm{act}}$ interval",
    "H5": "$b_A\\ge1$ satisfies \\eqref{eq:lut_capacity}", "H6": "authority \\eqref{eq:control_authority}",
    "H7": "initial conditions", "H_sat": "$L_{\\mathcal{T},z}\\sqrt{K}D_{\\mathrm{ODD}}/2\\le u_{\\max}$ \\eqref{eq:sat_cond}",
    "timing_ok": "latch order of \\eqref{eq:multirate_hold}", "codesign_window": "$b_A$ inside the window",
    "theorem_valid": "all of the above", "prop3_valid": "(i) $\\mu_R>L_{\\mathrm{sec}}$, (H-dz) and (H3)", "mu_R_gt_L_sec": "Prop.~\\ref{prop:deadzone_stability}(i) contraction rate", "Hdz": "$r_{\\mathrm{ss}}<r_{\\min}$ (Prop.~\\ref{prop:deadzone_stability} only)",
}

_REF = {"eq": r"\eqref{%s}", "lem": r"Lemma~\ref{%s}", "prop": r"Prop.~\ref{%s}", "thm": r"Thm.~\ref{%s}",
        "rem": r"Remark~\ref{%s}", "assump": r"Assumption~\ref{%s}", "cor": r"Cor.~\ref{%s}", "sec": r"Sec.~\ref{%s}",
        "def": r"Def.~\ref{%s}", "app": r"App.~\ref{%s}"}


def tex(key: str) -> str:
    """LaTeX name of a quantity / hypothesis key (KeyError if the key is unknown: no silent raw identifier)."""
    return NAMES[key]


def ref(label: str) -> str:
    """Equation / lemma label -> \\eqref{...} / Lemma~\\ref{...}; '--' for none."""
    if not label or label == "--":
        return "--"
    kind = label.split(":", 1)[0]
    return _REF[kind] % label


def hyp_list(keys) -> str:
    """Failing-hypothesis list typeset as (H4), (H6), ... ('none' if empty)."""
    out = []
    for k in keys:
        out.append(NAMES.get(k, k.replace("_", " ")))
    return ", ".join(out) if out else "none"


def num(x, digits: int = 4) -> str:
    """Float -> math-mode LaTeX number (scientific as a\\times10^{b}); nan -> N/A; inf -> $\\infty$."""
    x = float(x)
    if math.isnan(x):
        return "N/A"
    if math.isinf(x):
        return r"$\infty$" if x > 0 else r"$-\infty$"
    if x == 0:
        return "$0$"
    ax = abs(x)
    lo, hi = 10 ** -(digits - 1), 10 ** digits
    if lo <= ax < hi:
        s = f"{x:.{digits}g}"
        return f"${s}$"
    m, e = f"{x:.{digits - 1}e}".split("e")
    m = m.rstrip("0").rstrip(".")
    return f"${m}\\times10^{{{int(e)}}}$"


def boolean(b) -> str:
    return r"\checkmark" if bool(b) else r"$\times$"
