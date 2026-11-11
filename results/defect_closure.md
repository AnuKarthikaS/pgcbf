# revision 3 defect closure (repo side)

One entry per defect-register ID owned (fully or jointly) by the rectifier. "Before" is the revision 2
repository (`repo_baseline`, final revision 2 run); "after" is the final revision 3 run
(`results/logs/*.json`). Every "guard" is an assertion in an experiment ledger (fails the run) or a
pytest test.

## S1: tail coupling inactive (gate G-D)

- Change. `lattice.Lattice.__init__`: graded advection filter W_f = Phi_f diag(w) Phi_f^T with
  K_f = `plant.n_adv_modes` = 4 > K_m = 2, full weight on the first `plant.n_adv_full_weight` = 2
  modes and cutoff weight w = `plant.adv_cutoff_weight` on modes 3-4 (config literal; W_f = Pi_f at
  w = 1). A_full stays the pure damped-wave operator; the advection is moved into
  F_nl(x) = A_adv x + F_sin(x) (`Lattice.A_adv`, `Lattice.F_nl`). `rom.ROM`: new block
  A_adv_rperp (review round: corrected to P_0^{-1} Phi^T P_full A_adv P_perp, see item 3 below). `certificates.L_perp_P`: returns linear
  part, sine part, ||A_adv,r perp||_2 and ||A_r perp||. `certificates.certificate`: tail share of
  delta_total. New `experiments/exp9_wsweep.py` (w-sweep, frontier w*, plain-projector Pe sweep).
  `config/*.json`: T_s 5e-4 -> 3e-4 in R2 (latches scaled), `experiments.wave2_T_s` for the refusal row.
- Guards. exp1: "Q_{gamma/2} = gamma P_{gamma/2} (rel < 1e-12)", "(H2a) by invariance",
  "A_{r,perp} = 0 (class (H) eigenmode invariance)", "||A_adv,r perp||_2 > 0 in R2",
  "L_perp_P(R2) materially larger than L_perp_P(R1)", "tail share reported", refusal row
  "revision 2 sampling period" (certificate refuses at T_s = 5e-4). exp9: frontier exists, accept set is a
  down-set of the grid, deployed w equals the frontier, ablation (a) still refused at w*, open loop
  exits at every w, plain-projector disjointness. tests/test_no_literals perturbation of
  `plant.adv_cutoff_weight`.
- Before -> after.

| quantity | revision 2 | revision 3 R2 | revision 3 R1 |
|---|---|---|---|
| K_f (advected modes) | 2 (= K_m) | 4 | 4 |
| ||A_adv,r perp||_2 | 0 | 0.1075 | 0.005376 |
| L_perp_P | 8.719e-4 (both regimes) | 0.04741 (lin 0.04654 + sin 8.72e-4) | 0.003199 (lin 0.002327) |
| L_perp_P(R2)/L_perp_P(R1) | 1.0 | 14.8 | |
| eps_basis_P | 0.03831 | 0.09831 | 0.007293 |
| tail share of delta_total | about 0 | 0.282 | 0.0079 |
| delta_total / delta_cert | 0.01361 / 0.01751 | 0.01299 / 0.01456 (margin 0.00157) | 0.002321 / 0.03161 |
| T_s | 5e-4 | 3e-4 (5e-4 refused: delta_total 0.01726 > 0.01455) | 2e-3 |
| codesign window | [7, 11] | [8, 11] | [8, 11] |

- w-sweep (exp9, 17 values 1 -> 0.0039): w >= 0.0884 fails (H4), (H6), (H7); 0.03125-0.0625 fail
  (H6); every w <= 0.0221 accepts; ablation (c) is no longer refused for w <= 0.0156, so the frontier
  is w* = 2^{-11/2} = 0.0221 (deployed). Tail share 0.997 at w = 1 -> 0.282 at w*. Open loop at c = 1
  exits 16/16 at every w. Plain projector: largest certifiable c = 0.02 (Pe 0.4), smallest unstable
  c = 0.5 (Pe 10): disjoint. Structural finding stated in Section V without softening.

## S2: class-(H) lemma, nan in Table VI

- Change. `certificates.delta_cert_bound_classH` docstring and tables cite Lemma
  `lem:delta_cert_classH`; it reports c_M, lambda_min(M, P) and s_nl. `experiments/exp3_delta_cert.py`
  prints `N/A` (typeset `--`) for non-applicable cells instead of `nan`; `tools/redteam_gates.py`
  checks the lemma label.
- Guards. exp3 "class-(H) bound <= brute-force shell minimum" (R1, R2); tests/test_soundness::
  test_delta_cert_below_brute_force; tests/test_tex_lint::test_generated_table_is_clean (fails on nan).
- Before -> after. Table VI: nan cells -> none. c_M (R2) 0.3725 (unchanged), lambda_min(M, P) 0.3731,
  s_nl 3.0e-4 now printed; R1 c_M 0.9829.

## S3: Path A on the port vector

- Change. `controller.Policy`, `certificates.port_box`: the port-space LUT is the default
  (`controller.lut_space = "port"`), ports named v_1, v_2 in tables; the state-space LUT (m = K)
  appears only as the refused variant in exp1. D_v in closed form
  2 sqrt(E_safe) max_j sqrt((g_r^T P g_r)_jj).
- Guards. exp1 refusal row "state-space LUT (m = K)" (codesign window empty); exp1 "D_v closed form
  equals sampled sup"; exp2 port-box CROWN soundness.
- Before -> after. D_v: sampled only -> closed form 0.447210 vs sampled 0.447169.

## S4: hypothesis renumbering

- Change. `certificates.THEOREM_HYPS` = (H0)-(H7) (H5 = LUT capacity / codesign, H6 = authority,
  H7 = initial conditions), `Hdz` for Prop. 3, flags `theorem_valid` and `prop3_valid`
  (`certificates.certificate`, `refuse_if_invalid`); `texnames.HYP_DESC`; exp1/exp6 tables and prose.
- Guards. exp1 "theorem_valid in R2" and "theorem_valid in R1", "prop3_valid true in R1, false in R2";
  tests/test_refusal::test_deployed_configs_accepted; tools/redteam_gates.py (no `(H8)`, `(H5)` never
  the deadzone); test_tex_lint (check mark / cross, no True/False).
- Before -> after. Table II `theorem_full_valid = False` for R2 -> single `theorem_valid` column,
  check mark on R1 and R2; prop3_valid check mark R1, cross R2.

## S5: trivial stochastic positive case

- Change. `experiments/exp7_stochastic.py`: sigma grid search over `stochastic.sigma_grid`
  [0.002, ..., 0.02], z(0) at level `stochastic.z0_level` = 0.9, fluctuation criterion
  max_t |h(x_t) - h(x_0)| >= 0.1 h(x_0) on >= 50 % of screen paths; noise-induced deviation measured
  against the noise-free path; Corollary constants split into two tables.
- Guards. exp7 "sigma* is the largest accepted grid value", "fluctuation criterion", "noise itself is
  non-trivial (median noise-induced deviation >= 0.1 h0)", "CP upper <= alpha", ablation (f)
  assertions.
- Before -> after. sigma 0.003, z0 = 0, noise deviation 0.3 % of h0 -> sigma* = 0.01 (0.014 and 0.02
  refused), fluctuation share 1.0, median noise-induced deviation 0.246 h0. Exits 0/2000 (CP upper
  0.00184) vs Corollary bound alpha = 0.05; T_TTL 2438 vs horizon 2 (ratio 1219): the bound is very
  conservative, stated in the prose. Ablation (f): 1875/2000 -> 2000/2000 exits (CP lower 0.998).

## S6: Theorem 1 tightness

- Change. `experiments/exp4_theorem1.py`: second initial set on h_r = h_act/2 with zero tail under
  the worst-case adversary; reports min_t h, its time, and min h / Delta h_dip.
- Guards. exp4 "boundary start: min h = Delta h_dip at t = 0", "interior start: min h >= 0",
  "ratio > 1", "tail <= eps_basis_P".
- Before -> after. Only the boundary run (min h = Delta h_dip at t = 0; R2 5.21e-6) -> boundary
  R2 1.888e-6 / R1 2.382e-6 at t = 0; interior min h R2 0.01034 (ratio 5476), R1 0.01280 (ratio 5373),
  attained at t = 0 in both: the controller keeps h non-decreasing from the interior, so the
  certificate is not tight along these trajectories (stated in prose).

## S7: deadzone ablation at the benchmark b_A

- Change. `experiments/exp5_deadzone.py`: ablation (d) swept over b_A in {3, 5, 7, 9} (9 = deployed);
  new figure exp5_deadzone_sweep.
- Guards. exp5 "deadzone gives u_A = 0 late at every b_A", "mid-rise chatters at every b_A including
  the deployed one", "b_A = 1 violates (H-dz)".
- Before -> after. Only b_A = 3 (about 950 switches) -> switches per unit time 953, 860, 469, 122 at
  b_A = 3, 5, 7, 9; lag at 9 = 3.5e-4; deadzone 0 at all four.

## V1: raw identifiers and booleans

- Change. `src/pgcbf/texnames.py` (NAMES, HYP_DESC, `tex`, `ref`, `hyp_list`, `boolean`);
  `experiments/_common.write_table` routes every header and cell through it; `tools/texlint.py`.
- Guards. tests/test_tex_lint (linter not vacuous, every table clean, section clean);
  tools/page_sweep.py on the standalone PDF (logs/page_sweep.log).
- Before -> after. Page sweep pages 18-24 with raw identifiers / True / False / nan -> 0 hits.

## V2: unreadable Table XII, resizebox

- Change. No `\resizebox`; `\footnotesize` minimum; tables <= 7 columns unless typeset as
  `table*`; Corollary constants split into exp7_corollary_det and exp7_corollary_stoch; multirate
  table trimmed.
- Guards. tests/test_tex_lint::test_tables_have_no_resizebox_and_footnotesize_minimum,
  ::test_tables_have_at_most_seven_columns_unless_wide; page_sweep small-font check.
- Before -> after. 97 spans < 5.5 pt -> 0 (page_sweep).

## V3: clipped labels, caption mismatch, full-page floats

- Change. `plotting.save`: constrained_layout, IEEE widths (3.5 / 7.16 in), titles <= 60 chars,
  artist bounding-box check; captions in the template rewritten to match each layout.
- Guards. tests/test_figures (checker detects a clipped artist; refuses long titles / non-IEEE
  widths; every figure renders clean), page_sweep clipped-artist check.
- Before -> after. Fig. 3 x-label clipped, Fig. 8 caption "left/right" on stacked panels -> 0 clipped
  artists in the 10 figures; captions verified against the rendered PNGs.

## V4: prose numbers typed by hand

- Change. `tools/build_section.py` reads the CSV tables and ledgers, writes
  `results/section_numbers.json`, and fills `results/experimental_section_template.tex` into
  `results/experimental_section.tex`; it exits nonzero on any unresolved placeholder or lint hit;
  run by run_all.sh after every run.
- Guards. build_section exit status; tests/test_tex_lint::test_experimental_section_is_clean.
- Before -> after. Hand-typed numbers -> every number in the prose comes from a placeholder.

## V5 / V6

- V5: Section V uses `\eqref{}` / `Lemma~\ref{}` with the revision 1 labels and the new hypothesis
  numbering (checked by tools/redteam_gates.py).
- V6: `tools/redteam_gates.py` and `tools/page_sweep.py` (via `tools/standalone_pdf.py`, tectonic
  compile of Section V alone) are run by run_all.sh; any hit gives a nonzero exit.


## Review round (hostile review, Final Verdict Table items 1, 3-17, 19-22, plus cross-check patches)

"Before" is the delivered revision 3 run (revision-3 release); "after" is the review-round run
(`results/logs/run_all.log`). Items 2 and 18 live in the combined manuscript and were applied by the parent.

| item | review round | change | guard |
|---|---|---|---|
| 1 | DONE | `\pgcfig` macro in the generated section: natural size, capped at 0.62 textwidth, linewidth and 0.72 textheight with keepaspectratio, and it writes the size it used to the log. Figure widths 3.5 in / 4.4 in (4.4 <= 0.62 x 7.16 in). Titles wider than their axes are wrapped (`plotting.wrap_titles`). Legends were moved off the data in exp4, exp5 and exp9 | `tools/standalone_pdf.py` compiles a one-column (primary) and a two-column IEEEtran build and fails on any figure above the caps, any overfull box > 1 pt, overlapping text, objects in the bottom margin, or page-sweep hits. tests/test_figures: size caps, overlapping-axes checker, in-axes-legend checker, title wrapping |
| 3 (BLOCKING) | DONE | `rom.ROM.A_adv_rperp = Zmap @ A_adv @ Pperp` = P_0^{-1} Phi^T P_full A_adv P_perp. Printed as $\|A_{\mathrm{adv},r,\perp}\|_2$ and cited as Lemma lem:delta_cert_classH in exp1, exp9 and the prose | exp1 "||A_adv,r,perp||_2 > 0 in R2" (kept) |
| 4 | DONE | exp3 "N/A (standalone)"; exp6 eps*_adv moved to the caption; `write_table` sets tabcolsep to 3 pt for tables with >= 8 columns | 0 overfull hboxes in both layouts (standalone_pdf) |
| 5, 6, 17, 22 | DONE | D_v and Delta_{q,A} -> eq:Dv_closed; M_v -> eq:eps_A_triangle; delta_cert branch -> eq:delta_cert_classH; D_ODD -> eq:Cr0; mu_{2,P}(A_{r,0}); prose cites Sections sec:problem--sec:theory | redteam_gates, test_tex_lint |
| 7, 9, 21 | DONE | `build_section.py` uses one formatter, `num()` with 4 significant digits, the same function as the tables. The R1 tail share is printed as "a share ... of delta_total = ... in R1" | test_tex_lint::test_every_prose_numeral_occurs_in_section_numbers (fails if a numeral in the prose is not in `results/section_numbers.json`); build_section exits nonzero on unsourced numerals |
| 8 | DONE | Table I/II rows: N_m, N, m, nu, alpha, gamma, kappa, c, K_m, K, full-weight modes, E_safe, eps_adv, open-loop growth, L_perp linear share, D_v (sampled sup), tau_A, tau_B, samples per period, b_B, gamma_0, hold share, gamma_0 h_act. Table XI has a final-lag column; Table XV has h(x_0) | same lint test |
| 10 | DONE | Prose: the plain-projector result is the L^2-ODD conservatism of Remark rem:rda_cutoff. Keeping A_adv inside F_nl preserves (H2a) exactly (Lemma lem:h2a_char) at the price of charging A_{adv,r,perp} in M_perp and L_{perp,P}. The weight w attenuates the cutoff modes the way an H^1-bounded ODD would for a smooth field. Certifying w = 1 at Pe >= 10 would need that H^1 formulation, which is outside the theorem | exp9 disjointness assertion |
| 11 | DONE | Prose: Delta h_dip is attained at t = 0 by construction, h increases from t = 0+ because delta_cert > delta_total, so the runs do not exercise the inter-sample dip dynamically. Stated as a limitation | exp4 assertions (t_min = 0, ratio > 1) |
| 12 | DONE | exp7 computes gamma_0 h_act and delta_cert - delta_total on R2 and runs the Corollary at sigma = 0 on R2. Both are printed, and the prose states that the Corollary is infeasible on R2 for that reason and that r_W2 = 0 | exp7 "R2: Corollary refuses at sigma = 0", "R2: gamma0 h_act > margin"; exp1 "R2 gamma0 h_act > margin" |
| 13 | DONE | samples per period 2 pi/(omega_min T_s) and hold share of delta_total are computed (Table II) and printed | exp1 table rows |
| 14 | DONE | Table III row "mu_R > L_sec" and Prop. 3 row "(i) mu_R > L_sec, (H-dz) and (H3)". Prose sentence for R1 (holds) and R2 (fails) | exp1 "R1 mu_R > L_sec", "R2 mu_R <= L_sec" |
| 15, 16 | DONE | Fig. exp5z cited; \mathcal C_{r,act} | redteam_gates (every figure label cited) |
| 19 | DONE | exp2 ylabel $\|\Psi_B(z)-\hat{\mathcal{T}}(z)\|_2$ | rendered PNG inspected |
| 20 | DONE | exp9 non-evaluated "abl. refused" cells -> "--"; exp6 titles "(c1) single rate ($T_s=0.2$)" / "(c2) multi rate ($T_s=0.2$)" | rendered PNG inspected |
| X1 (cross-check 1) | DONE | same as item 3 | as item 3 |
| X2 (cross-check 2) = R-L1 | DONE | `certificates.certificate`: L_Lfh = 2 ||sym P_ih|| sqrt(E) + ... (missing factor 2 restored) | tests/test_soundness::test_Lfh_lipschitz_dominates_finite_differences (random pairs in C_{r,0} plus finite differences at the boundary maximiser, both regimes; it also asserts the empirical quotient exceeds the pre-fix constant, so it would have caught the bug) |

### R-L1: soundness defect found after revision 2

- Cause. Factor 2 missing in the gradient of a quadratic form: the Lipschitz modulus of
  L_{f_r} h_r on C_{r,0} used ||sym P_ih|| sqrt(E) where the gradient of z^T M z is 2 sym(M) z. The
  certified L_hold, and hence delta_total, were under-estimates (unsound). Found by the parent's
  cross-check (`verify_crosscheck.py`).
- Fix. `src/pgcbf/certificates.py` (L_Lfh, comment "R-L1"). Guard:
  `tests/test_soundness.py::test_Lfh_lipschitz_dominates_finite_differences`.
- Before -> after (items 3 and R-L1 together):

| quantity | R2 before | R2 after | R1 before | R1 after |
|---|---|---|---|---|
| ||A_adv,r,perp||_2 | 0.1075 (P_full,eps) | 0.1060 | 0.005376 | 0.005301 |
| Lip(L_f h) | 0.8855 | 1.770 | 0.3980 | 0.7948 |
| L_hold | 40.22 | 44.04 | 1.161 | 1.474 |
| delta_total | 0.01299 | 0.01356 | 0.002321 | 0.002634 |
| delta_cert | 0.01456 | 0.01456 | 0.03161 | 0.03161 |
| (H6) margin | 0.001573 | 9.994e-4 | 0.02929 | 0.02898 |
| tail share of delta_total | 0.2823 | 0.2704 | 0.007906 | 0.006967 |
| hold share of delta_total | 0.4644 | 0.4872 | 0.5002 | 0.5596 |
| refusal row: delta_total at the revision 2 T_s = 5e-4 | 0.01726 | 0.01821 (still refused) | | |
| exp4 boundary min h = Delta h_dip | 1.888e-6 | 2.060e-6 | 2.382e-6 | 3.008e-6 |
| exp4 interior min h / Delta h_dip | 5476 | 5019 | 5373 | 4255 |
| exp6 largest certifiable eps*_adv | 0.005501 | 0.004223 | | |
| exp7 sigma* (R1) | | | 0.01 | 0.007 |
| exp7 T_TTL (ratio to T) | | | 2438 (1219) | 2022 (1011) |
| exp7 median noise-induced deviation | | | 0.246 h0 | 0.174 h0 |

- Verdicts. No sign change: the (H6) margin stays positive in R1 and R2, theorem_valid holds in both,
  the w-sweep frontier w* = 0.0221 and every exp9 accept/refuse cell are unchanged, the plain-projector
  disjointness (certifiable c <= 0.02, unstable c >= 0.5) is unchanged, and every ablation still breaks.
  The one verdict that changed is the exp7 grid point sigma = 0.01. The Corollary accepted it before and
  now refuses it (failing condition: margin), so sigma* drops from 0.01 to 0.007. The positive stochastic
  case still has 0/2000 exits (CP upper 0.001843 <= alpha = 0.05). The generated Section V states the
  correction and this changed verdict. Blank cells mean the value was not printed in the earlier release.
- Suite. 200/200 experiment assertions (revision 3: 195), pytest 58 passed + 2 skipped, two-layout
  standalone build: 0 page-sweep issues, 0 overfull boxes, 10/10 figures within the caps.


## revision 4 (Phase 4.1): what changed and why

| item | change | why |
|---|---|---|
| Exp 4b | new experiment `exp4b_budget.py` (48 assertions), table `exp4b_budget`, keys `e4b_*`, Section V subsection | story bible C13: report the fraction of each budget term used along the certified trajectories, whatever it is, and check the proof inequality (v) at every node |
| simulator hook | `Simulator.run(..., node_hook=None)` | lets Exp 4b observe every integrator node without a second integrator; default path unchanged (same SHA-256 as the revision 3 code on this machine) |
| Exp 4b tail selection | modal (`T @ (Pperp @ Y)`) instead of a nodal difference | the nodal difference left a round-off residue at near-zero tails in R1 that made the diagnostic ||Delta_perp||/(L_perp ||x_perp||) exceed one; after the fix it is <= 1 at every node and is asserted |
| P1 | removed: Exp 4 ratio column (table + CSV), keys `e4_iratio_R2`, `e4_iratio_R1`, assertion "floor not tight from the interior"; Exp 4 prose rewritten | C12: no min_t h / Delta h_dip ratio, never "the dip is tight"; Exp 4b now carries the utilisation statement. The historical S6 entry above (ratios 5476 / 5373) and the review-round row "exp4 interior min h / Delta h_dip" are superseded |
| P2 | design note item 7: sigma* = 0.007; stale "median 0.25 h0" replaced by the run values (criterion 0.1 h0 on half the paths; median deviation keyed `e7_noisedev`) | consistency with the R-L1 correction |
| P4 | `LICENSE` (MIT), `CITATION.cff` | release metadata |
| Section V notation | Exp 4 and 4b prose avoid second-level subscripts (||x_perp|| with the weighting stated in words, d_r^T grad h~_r, L_hold ||z - z_k|| / M_z) | the first revision 4 run failed the page sweep (spans under 5.5 pt) on these symbols |
| main paper | keys `util_R1`, `util_R2`, `e6_c_delta`; `tools/build_main.py` (+ `tests/test_build_main.py`); `tools/main_figures.py`; CSVs `exp1_budget.csv`, `exp9_wsweep_deltas.csv`; new exp1 assertions "the six budget parts sum to delta_total" (R1, R2) | revision 4 main-paper build |

Theory-side constants: the diff of `results/section_numbers.json` against the revision 3 release shows only the two
removed P1 keys, the new keys above, and the run metadata `n_assert`, `n_tests`, `rt_total`.
