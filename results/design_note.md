# Design note: revision 3 numerical validation of the projection-guided multi-rate CBF certificate

Sources of truth: `paper/theory/sections_II_IV.tex` and `paper/theory/appendix_A_D.tex`, with the
fixed names and numbering of the revision 3 patch specification (hypotheses (H0)-(H7) for Theorem 1,
(H-dz) for Prop. 3 only, Lemma `lem:delta_cert_classH` for the class-(H) bound). Every numerical
constant is either a config literal (`config/default.json` = R2, `config/r1_passive.json` = R1) or is
computed at run time in float64 by a function in `src/pgcbf/` whose docstring names the equation it
implements. `tests/test_no_literals.py` (AST scan plus perturbation test) enforces this.

## 1. Benchmark (class (H) RDA lattice)

Full order (eq:full_system), N = 2 N_m = 320 states, coordinates x = (q, p):

    q' = p
    p' = -(nu L_h + alpha I) q - gamma p - c W_f D_h W_f q - kappa (q - sin q) + B u + B_w w

- `A_full = (J_full - R_full) P_full` is the pure damped-wave operator (Omega^2 = nu L_h + alpha I,
  J canonical, R = blkdiag(0, gamma I), class (H) with R_p = gamma I). It does NOT contain the
  advection. Therefore Q_{gamma/2} = gamma P_{gamma/2} holds to round-off (asserted, rel. error
  7.0e-15 < 1e-12), lambda_{K+1} = gamma/2 is unchanged, and (H2a) holds by eigenmode invariance
  (Lemma lem:h2a_char, residual 6.4e-13). The paper's A_{r,perp} = P_0^{-1} Phi^T P_full A_full P_perp
  is identically zero (asserted 1.3e-12); no positivity claim is made about it.
- The nonlinearity is written explicitly as `F_nl(x) = A_adv x + F_sin(x)` with
  A_adv = [[0, 0], [-c W_f D_h W_f, 0]] and F_sin(x) = [0; -kappa (q - sin q)]. F_nl(0) = 0 and F_nl
  is globally Lipschitz (Assumption assump:generic).
- **Graded advection filter.** W_f = Phi_f diag(w_1, ..., w_{K_f}) Phi_f^T with K_f = `n_adv_modes` = 4
  sine modes, w_k = 1 for k <= `n_adv_full_weight` = 2 and w_k = w = `adv_cutoff_weight` for the
  cutoff modes k = 3, 4. The plain projector Pi_f is the special case w = 1. D_h is the skew central
  difference, so W_f D_h W_f is skew. Because K_f = 4 > K_m = 2, the advection couples the ROM modes to
  tail modes 3 and 4: the advective cross-coupling block
  A_{adv,r,perp} = P_0^{-1} Phi^T P_full A_adv P_perp (Lemma lem:delta_cert_classH; review round: P_full,
  not P_{full,eps}) is nonzero (0.1060 in R2, 0.005301 in R1), and the tail
  coupling is active in L_perp_P, M_perp and delta_total (gate G-D).
- **w is a disclosed, swept benchmark parameter** (S1 decision). It appears in the Section V
  benchmark formula and in Table I. exp9 sweeps w over the geometric grid 2^{-k/2}, k = 0..16
  (1 down to 0.0039, 17 values), and the deployed w* = 2^{-11/2} = 0.0221 is the LARGEST grid value
  for which the full certificate accepts with a positive (H6) margin and all ablations (a)-(c), (e)
  still break (ablation (c) is no longer refused for w <= 0.0156). It was selected by the sweep,
  not pre-picked.
- **Structural finding (reported without softening).** Under the Euclidean-ODD supremum that
  defines M_perp, full-strength cutoff-mode advection (w = 1) is not certifiable at any Péclet number
  at which the open loop is unstable: with the plain projector the largest certifiable c is 0.02
  (Pe 0.4) and the smallest open-loop-unstable c is 0.5 (Pe 10), so the two sets are disjoint
  (asserted in exp9). This is the conservatism anticipated in Remark rem:rda_cutoff; it is a
  limitation of the certificate, not of the plant. At every w on the grid the open loop at c = 1
  leaves C for 16/16 seeds.
- B and B_w are Gaussian profiles with a Butterworth roll-off in the sine-mode index,
  f_k = (1 + (k/k_c)^{2p})^{-1/2}, (k_c, p) = (2, 12), renormalised to unit L2 norm.
- R2 (default): c = 1, Pe = 20, open-loop growth rate 0.292 > 0, 16/16 open-loop exits (exp0).
  R1 (`r1_passive.json`): c = 0.05, Pe = 1, passive; (H-dz) holds there.
- ROM: K_m = 2 mass-normalised eigenmodes, K = 4, m = 2 actuators.

## 2. Hypotheses and how each is checked at run time

| hyp. | statement (paper) | function | flag / check |
|---|---|---|---|
| (H0) | eq:H0 ODD containment, E < D^2 lambda_min(P)/4 | `certificates.odd_containment` | `H0` |
| (H1) | assump:generic constants on X_ODD | `lipschitz_Fnl`, `M_perp_bound`, `L_perp_P` | constants computed; F_nl(0) = 0 by construction |
| (H2) | (H2a) Phi^T Pfe P_perp = 0, (H2b) R_p,perp > eps, (H2c) lambda_{K+1} > 0, gamma < 2 omega_min, eps = gamma/2 | `ROM.h2a_residual`, `ROM.h2a_invariance_residual`, `ROM.h2b_margin`, `tail_rate`, `pencil_rate`, `eps_choice_check` | `H2` |
| (H3) | teacher eq:teacher_def, H/2 >= m, G_port > 0, (H-sat) eq:sat_cond | `teacher_constants` | `H3` (includes `H_sat`) |
| (H4) | tail tolerance eq:eps_basis_def, h_act in (L_h M_z T_s, E - eps_b^2 - dh) | `shell` | `H4` |
| (H5) | LUT capacity eq:lut_capacity, b_A >= 1; codesign window eq:codesign_window | `codesign_window` | `H5`, `codesign_window` (lower bound required: stricter than the paper) |
| (H6) | control authority eq:control_authority, delta_total <= delta_cert | `delta_total_exact`, `delta_cert_bound_classH`, `two_branch` | `H6`, margin reported |
| (H7) | initial conditions z(0) in C~_r, ||x_perp(0)|| within the tail tolerance | `simulate.initial_conditions` | `H7` |
| timing | 0 < tau_A < tau_B < T_s (eq:multirate_hold) or single rate | `certificate` | `timing_ok` |
| (H-dz) | eq:deadzone_cond r_ss < r_min (Prop. 3 only) | `rss_ball`, `sector_modulus` | `Hdz` |

Two flags summarise this: `theorem_valid` = (H0)-(H7) and timing and codesign window, which is TRUE
in both R1 and R2 (asserted); `prop3_valid` = (H-dz) and (H3), which is true in R1 and false in R2
(L_sec = 2.39 > mu_R = 0.5). (H-dz) is not part of Theorem 1. `certificates.refuse_if_invalid` raises
`CertificateRefused` unless `theorem_valid` holds; every experiment relying on Theorem thm:main_cbf
calls it.

## 3. Soundness choices for the constants

- **M_perp (eq:M_perp).** Triangle inequality: linear advective part ||P_perp A_adv|| rho_ODD, sine
  part via ||q - sin q|| <= ||q||^3/6, input and disturbance parts via operator norms times u_max and
  eps_adv.
- **L_perp_P (eq:Lperp_def).** Reported as two parts: the linear part ||A_{adv,r perp}||-type term
  (Pfe-weighted operator norm of the advective cross-coupling) and the sine part
  kappa (1 - cos rho_ODD). R2: 0.04741 = 0.04654 (linear) + 8.72e-4 (sine); R1: 0.003199
  (0.002327 linear). Ratio R2/R1 = 14.8.
- **eq:reduced_constants.** Suprema over the P-ellipsoid come from exact ellipsoid operator norms.
  D_v has the closed form D_v = 2 sqrt(E_safe) max_j sqrt((g_r^T P g_r)_{jj}) on
  C_{r,0} = {||z||_P^2 <= E_safe}; it is computed in closed form and asserted equal to the sampled
  supremum (0.447210 vs 0.447169).
- **delta_cert (Lemma lem:delta_cert_classH).** Class (H) is not class (D), so eq:delta_cert_closed
  does not apply verbatim. The sound bound uses
  M = gamma P - (P A_adv,r + A_adv,r^T P) + 2 eta_T P g_r g_r^T P, c_M = lambda_min(M, P) - 2 s_nl,
  and the two-branch infimum over the shell. Reported: c_M, lambda_min(M, P), s_nl (R2: 0.3725,
  0.3731, 3.0e-4). exp3 checks it against a brute-force shell minimum; the class-(D) closed form is
  checked on a synthetic instance.
- **eps_CROWN.** Twin-CROWN (eq:twin_crown_bound); the certificate uses min(twin-CROWN, eq:crown_explicit).
- **eps_A (eq:eps_A_triangle).** Port-space LUT is the default (rem:lut_variants(i)), with port
  names v_1, v_2; the state-space LUT (m = K) is reported as the refused variant (codesign window
  empty).

## 4. Chosen parameters and justification

| parameter | R2 | R1 | why |
|---|---|---|---|
| N_m | 160 | 160 | stiff enough for the integrator test, runs stay a few seconds |
| nu, gamma, alpha | 0.05, 1, 1 | same | omega_min = 1.22 so gamma < 2 omega_min; lambda_{K+1} = gamma/2 |
| c (Pe) | 1.0 (20) | 0.05 (1) | R2 open-loop unstable; R1 passive |
| K_f, n_adv_full_weight, w | 4, 2, 0.0221 | same | K_f > K_m activates the tail coupling; w selected by exp9 |
| kappa, E_safe | 0.05, 0.05 | same | sine nonlinearity active on the ODD (rho = 0.28 rad) |
| gamma0 | 0.35 | 0.3 | delta_cert in branch 1 with slack |
| eps_adv | 0.002 | 0.002 | below the largest certifiable eps* = 0.0055 (exp6 e) |
| L_T, u_max | 5, 3.4 | 0.5, 0.35 | u_max just above the (H-sat) left-hand side |
| b_A, b_B | 9, 12 | 9, 12 | b_A inside the codesign window [8, 11] |
| T_s, tau_A, tau_B | 3e-4, 1.5e-5, 3e-5 | 2e-3, 1e-4, 2e-4 | see trade-off below |
| MC paths | 2000 | -- | budget |

**S1 trade-off (honest).** Activating the tail coupling (K_f = 4) raises delta_total in R2. At the
revision 2 sampling period T_s = 5e-4 the certificate refuses ((H6): delta_total 0.01821 > delta_cert
0.01455 after the review-round correction R-L1, 0.01726 before it; asserted as an exp1 refusal row
"previous sampling period"). After R-L1, 4e-4 (delta_total 0.01595) and 3.5e-4 (0.01482 > 0.01456)
also fail (before R-L1, 3.5e-4 passed). The deployed T_s is 3e-4 (latches scaled with it), i.e. a 1.67x
faster loop than revision 2. The (H6) margin at the deployed point is 9.994e-4 (delta_total 0.01356,
delta_cert 0.01456; 0.00157 before R-L1). Pe is unchanged (20).

## 5. Tuning record

1. Advection: revision 2 used Pi_f with n_adv_modes = K_m = 2, so the tail coupling was inactive (S1).
   revision 3 uses K_f = 4 with the graded W_f; w is swept (exp9) and the deployed value is the frontier.
2. Actuator and disturbance profiles are band-limited (Butterworth) to reduce tail injection.
3. Teacher theta = 0.05 with L_T, u_max set just above the (H-sat) left-hand side.
4. The state-space LUT has an empty codesign window (refused); the port-space LUT (K -> m = 2) is used.
5. T_s: 5e-4 -> 3e-4 (see Section 4).
6. R1 T_s = 2e-3, chosen so the stochastic Corollary margin is positive.
7. Stochastic positive case (S5): initial condition at level 0.9 of the shell, sigma chosen by a grid
   search over [0.002, ..., 0.02] as the largest value the Corollary accepts (sigma* = 0.007), with a
   fluctuation criterion asserting that the noise-induced deviation is non-trivial (max_t |h - h0| >= 0.1 h0
   on at least half of the paths; at sigma* the median deviation from the noise-free trajectory is 0.174 h0,
   keys e7_ffrac and e7_noisedev; the earlier "median 0.25 h0" was stale).
8. Ablation (b) uses b_B = 2 (b_B = 3 with learned rounding gave no violation); reported.
9. (H-dz) is not satisfiable in R2 (L_sec = 2.39 > mu_R = 0.5); Prop. 3 is demonstrated in R1.
   This does not affect Theorem 1 (`theorem_valid` is true in R2).
10. The certificate requires the codesign lower bound, which the paper marks optional (stricter).

## 6. Numerical method

- Exponential integrator (ETD2RK / ETD1) in the orthonormal lattice eigenbasis; exp(Ah), phi_1, phi_2
  from one augmented matrix exponential.
- Latch instants tau_A, tau_B are integration nodes (exact latch timing).
- Worst-case adversary w* = eps_adv d_r^T P z / ||.|| (rem:hji_game).
- Stochastic runs sample the additive noise exactly per sub-step from its Van Loan covariance.
- Trajectories are raw and unclipped; seeds come from config; two BLAS threads.


## Review round (correction R-L1 and item 3)

- `rom.ROM.A_adv_rperp` is now P_0^{-1} Phi^T P_full A_adv P_perp (`Zmap @ A_adv @ Pperp`), as in Lemma
  lem:delta_cert_classH: 0.1060 (R2), 0.005301 (R1); previously P_{full,eps} gave 0.1075 / 0.005376.
- `certificates.certificate`: the Lipschitz modulus of L_{f_r} h_r on C_{r,0} restores the factor 2 of the
  gradient of the quadratic form (defect R-L1). L_Lfh 0.8855 -> 1.770 (R2), 0.3980 -> 0.7948 (R1);
  L_hold 40.22 -> 44.04 (R2), 1.161 -> 1.474 (R1); delta_total 0.01299 -> 0.01356 (R2),
  0.002321 -> 0.002634 (R1). (H6) margins remain positive (9.994e-4, 0.02898). New regression test
  `tests/test_soundness.py::test_Lfh_lipschitz_dominates_finite_differences`.
- Downstream: exp6 eps*_adv 0.005501 -> 0.004223; exp7 sigma* 0.01 -> 0.007 (sigma = 0.01 now refused on
  the margin); exp4 Delta h_dip 1.888e-6 -> 2.060e-6 (R2). The w-sweep frontier w* = 0.0221 and all exp9
  verdicts are unchanged. See `results/defect_closure.md`, "Review round".


## revision 4 (Phase 4.1)

- **Experiment 4b (budget utilisation, `experiments/exp4b_budget.py`).** Re-runs the Exp 4 closed loops (both
  regimes, both starts, all seeds, worst-case adversary) with `simulate.Simulator.run(node_hook=...)`, which hands the
  state, the held inputs, the sample z_k and the adversary to a recorder before every integrator step; the integration
  itself is unchanged (the determinism hash of `tests/test_determinism.py` equals that of the revision 3 code on the same
  machine). At each node the recorder evaluates every term of eq:ode_diff_ineq (realised tail coupling, adversary,
  quantization/lag, hold shift) next to the bound the proof uses for it, the soundness slack (v), and diagnostic
  factor ratios (||2Pz||/L_h, ||Delta_perp||/(L_perp ||x_perp||), cosine, e_u/ebar_u, ||z - z_k||/(M_z s), hold shift
  over its Lipschitz bound). Definitions: the Experiment 4b section of the supplement (derivation cross-checked independently). Z(t) and
  z_k use the same map (`Zmap @ X`), so the interior-branch slack is exactly 0 at s = 0; s = 0 nodes are masked
  out of the ratios that divide by s; s is accumulated from the segment ends (0, tau_A, tau_B). The tail is selected
  modally (`T @ (Pperp @ Y)`); a nodal difference gave a round-off artefact at near-zero tails in R1.
- Results (keys `e4b_*`, Table `exp4b_budget`): every realised/bound ratio is <= 1 and the slack (v) is positive at
  every active-shell node in both regimes. The adversary term uses the largest fraction of its bound, the hold
  shift the smallest; the conservative constants are named in the generated Section V prose.
- **P1.** Exp 4 no longer reports min_t h / Delta h_dip ratios (table column, CSV column, keys `e4_iratio_*`, and the
  assertion "floor not tight from the interior" are removed). The prose states only: min of h at t = 0 for both
  starts, (a)-(c) at every node, the dip floor not exercised dynamically, and points to Exp 4b.
- **P2.** Tuning record item 7: sigma* = 0.007 (and the stale "median 0.25 h0" replaced, see item 7).
- **P4.** `LICENSE` (MIT, 2026, The pgcbf authors) and `CITATION.cff` (placeholder author).
- **Main-paper support.** Keys `util_R1`, `util_R2` (delta_total/delta_cert, 3 significant digits) and `e6_c_delta`
  (single-rate delta_total at the deployed T_s); `tools/build_main.py` (any template, same placeholder resolution
  from `results/section_numbers.json`, fails on unresolved placeholders, unused new keys, untraceable numerals and
  texlint findings); `tools/main_figures.py` renders `main_fig1_ablation` and `main_fig2_budget` (6.0 in, jmlr text
  width) from `exp6_ablations.csv`, `exp4_theorem1.csv`, `exp1_budget.csv` (new: the six budget parts per regime,
  asserted to sum to delta_total) and `exp9_wsweep_deltas.csv` (new: delta_total, delta_cert per w).
