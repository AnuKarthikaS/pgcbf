# Traceability: paper label → function → test → experiment → figure/table

Functions are in `src/pgcbf/`, tests in `tests/`, experiments in `experiments/`. Tables are in
`results/tables/` (CSV + .tex) and figures in `results/figures/` (PDF + PNG). Hypothesis keys follow
the revision 3 numbering: (H0)-(H7) for Theorem 1, (H-dz) (key `Hdz`) for Prop. 3 only; summary flags
`theorem_valid` and `prop3_valid`. All LaTeX names of code identifiers come from `pgcbf.texnames`.

| label (paper) | function(s) | pytest | experiment / assertion | figure / table |
|---|---|---|---|---|
| eq:full_system | `lattice.Lattice` (A_full pure damped wave, A_adv, F_sin, G, D), `simulate.Simulator` | test_soundness::test_theorem1_short_run | exp0, exp4, exp6 | exp0_openloop_h, exp4_theorem1 |
| assump:generic (F_nl = A_adv x + F_sin) | `lattice.Lattice.F_nl`, `certificates.lipschitz_Fnl` | test_no_literals::test_every_certificate_quantity_depends_on_the_inputs | exp1 (L_Fnl) | exp1_constants |
| eq:h_def | `lattice.Lattice.h_barrier` | test_identities::test_quad_split_1e12 | exp0, exp4 | exp0_openloop_h, exp4_theorem1 |
| eq:eps_choice | `certificates.eps_choice_check`, `Lattice.eps` | test_identities | exp1 (H2) | exp1_hypotheses |
| eq:mu_p0, eq:mu_pe, prop:dahlquist | `certificates.mu_2P`, `certificates.dahlquist` | test_identities::test_prop2_mu_equals_minus_gamma_half | exp1 "Q_{gamma/2} = gamma P_{gamma/2} (rel < 1e-12)" | exp1_constants |
| eq:tail_dissipation, (H2c) | `certificates.tail_rate` | test_identities::test_tail_rate_equals_pencil_rate | exp1 "tail rate" | exp1_constants |
| lem:pencil_rate, rem:tail_rates | `certificates.pencil_rate` | test_identities::test_tail_rate_equals_pencil_rate | exp1 | exp1_constants |
| (H2a), lem:h2a_char | `rom.ROM.h2a_residual`, `rom.ROM.h2a_invariance_residual` | test_identities::test_quad_split_1e12 | exp1 "(H2a)", "A_{r,perp} = 0 by invariance" | exp1_constants |
| (H2b) | `rom.ROM.h2b_margin` | test_refusal | exp1 hypothesis (H2) | exp1_hypotheses |
| eq:pg_coords, eq:reduced_dynamics | `rom.ROM` (Zmap, Pperp, P, Ar0, g_r, d_r, F_nl_r) | test_identities | exp1, exp8 (modal vs nodal) | exp8_numerics |
| eq:quad_split, lem:quad_split | `rom.ROM.hPhi`, `lattice.Lattice.h_barrier` | test_identities::test_quad_split_1e12 | exp1 | exp1_constants |
| eq:H0, lem:odd_containment, (H0) | `certificates.odd_containment` | test_refusal::test_refusals | exp1 | exp1_hypotheses |
| eq:M_perp | `certificates.M_perp_bound` | test_no_literals (perturbation), test_soundness::test_theorem1_short_run | exp4 "tail <= eps_basis_P" | exp4_theorem1 |
| eq:eps_basis_P | `certificates.eps_basis_P` | same | exp4, exp6 (a), exp9 | exp4_theorem1, exp9_wsweep |
| eq:Lperp_def (linear + sine parts) | `certificates.L_perp_P` | test_no_literals | exp1 "L_perp_P(R2)/L_perp_P(R1)", "A_{adv,r perp} > 0 in R2", exp9 | exp1_constants, exp9_wsweep |
| lem:delta_cert_classH: A_{adv,r,perp} = P_0^{-1} Phi^T P_full A_adv P_perp (review item 3: P_full, not P_{full,eps}) | `rom.ROM.A_adv_rperp` (= `Zmap @ A_adv @ Pperp`), `certificates.L_perp_P` (key `A_adv_rperp_norm`) | test_no_literals | exp1 "||A_{adv,r,perp}||_2 > 0 in R2", exp9 | exp1_constants, exp9_wsweep |
| Lipschitz modulus of L_{f_r} h_r on C_{r,0} (L_Lfh, enters L_hold of eq:delta_total_exact; defect R-L1, factor 2 of the quadratic-form gradient) | `certificates.certificate` (L_Lfh, L_hold) | test_soundness::test_Lfh_lipschitz_dominates_finite_differences (random pairs + finite differences at the boundary maximiser) | exp1 (L_Lfh, L_hold, delta_total) | exp1_constants |
| eq:reduced_constants, D_v closed form | `certificates.reduced_constants`, `certificates.port_box` | test_no_literals | exp1 "D_v closed form = sampled sup" | exp1_constants |
| eq:teacher_weights, eq:teacher_def, eq:groupsort_relu, (H3) | `controller.Teacher`, `controller.project_U`, `controller.gs2`, `certificates.teacher_constants` | test_soundness | exp2, exp3 | exp2_crown |
| eq:LTz, lem:sector, eq:sat_cond (H-sat) | `certificates.teacher_constants` | test_refusal::test_refusals (u_max case) | exp1 refusal "u_max below the (H-sat) lhs" | exp1_refusals |
| eq:twin_slopes, eq:twin_affine_relax, eq:twin_crown_bound, prop:twin_crown | `crown.twin_slopes`, `crown.TwinCrown.eps_crown` | test_soundness::test_twin_crown_sound_1e5, ::test_twin_crown_monotone_in_bB | exp2 | exp2_crown |
| eq:crown_explicit | `crown.E_q` | test_no_literals | exp2 | exp2_crown |
| eq:eps_A_triangle | `certificates.eps_A_triangle` | test_no_literals | exp1, exp5 | exp1_constants, exp5_deadzone |
| eq:lut_capacity, eq:codesign_window, (H5) | `certificates.codesign_window` | test_refusal (state-space LUT) | exp1 refusal "state-space LUT (m = K)", exp5 b_A sweep | exp1_refusals, exp5_deadzone |
| rem:lut_variants (i) port-space LUT (default) | `certificates.port_box`, `controller.Policy` | test_soundness | exp2 "port box" | exp2_crown |
| eq:exp_weights | `certificates.exp_weights` | test_identities::test_exp_weights_and_phi | exp6 (c) | exp6_multirate |
| eq:multirate_hold, timing | `simulate.Simulator.run` | test_determinism::test_exact_latch_timing, test_refusal::test_timing_order_refused | exp8 | exp8_numerics |
| eq:delta_dip | `certificates.delta_dip` | test_no_literals | exp4 (boundary min h = Delta h_dip) | exp4_theorem1 |
| eq:delta_total_exact, eq:exact_exp_integral | `certificates.delta_total_exact`, `certificates.phi_hold` | test_identities::test_exp_weights_and_phi | exp1, exp4, exp6, exp9 (tail share) | exp1_constants |
| eq:h_tilde_def, eq:eps_basis_def, eq:shell, (H4) | `certificates.shell`, `certificates.sample_shell` | test_soundness::test_delta_cert_below_brute_force | exp3 | exp3_delta_cert |
| eq:teacher_cbf, eq:null_space_cond | `certificates.bracket`, `certificates.null_space_rate` | test_soundness::test_delta_cert_below_brute_force | exp3 | exp3_delta_cert |
| lem:delta_cert_classH (class-(H) bound; c_M, lambda_min(M,P), s_nl) | `certificates.delta_cert_bound_classH` | test_soundness::test_delta_cert_below_brute_force | exp3 (i) "bound <= brute force" | exp3_delta_cert |
| eq:delta_cert_closed, lem:delta_cert_closed | `certificates.two_branch`, `certificates.classD_closed_form` | test_soundness::test_classD_two_branch_matches_affine_brute_force | exp3 (ii) class (D), both branches | exp3_delta_cert |
| eq:control_authority, (H6) | `certificates.certificate` (flag `H6`, margin), `refuse_if_invalid` | test_refusal | exp1 refusals ("revision 2 sampling period"), exp6, exp9 | exp1_refusals, exp6_ablations, exp9_wsweep |
| (H7) initial conditions | `simulate.initial_conditions` | test_soundness::test_theorem1_short_run | exp4 boundary and interior starts | exp4_theorem1 |
| thm:main_cbf (a)-(c), `theorem_valid` | `certificates.certificate`, `simulate.Simulator` | test_refusal::test_deployed_configs_accepted, test_soundness::test_theorem1_short_run | exp1 "theorem_valid in R1 and R2", exp4 (16 seeds, boundary + interior; min of h at t = 0, (a)-(c) at every node) | exp4_theorem1, exp1_hypotheses |
| thm:main_cbf proof chain along trajectories (revision 4, Exp 4b): decomposition H = C + A + Q + r_hold, (i) tail <= L_h L_perp eps_basis_P, (ii) adversary <= d_h eps_adv, (iii) quantization/lag <= g_h ebar_u (phase-wise, eq:multirate_hold), (iv) hold <= L_hold s, (v) slack >= 0 | `simulate.Simulator.run(node_hook=...)`, `experiments/exp4b_budget.Recorder`; derivation in the Experiment 4b section of the supplement | (assertions in exp4b) | exp4b: node count, two decomposition identities to round-off, (i)-(iv) ratios <= 1, factor bounds <= 1, ||z - z_k|| <= M_z s, (v) slack >= 0 at every node, Teacher Theta(z_k) >= delta_cert, interior-branch slack | exp4b_budget |
| rem:hji_game | `simulate.Simulator.w_star` | test_soundness::test_theorem1_short_run | exp4, exp6 (e) | exp6_ablations |
| rem:rda_cutoff, graded filter W_f | `lattice.Lattice` (adv_cutoff_weight), `lattice.Lattice.peclet` | test_refusal | exp0, exp1 refusal "pure lattice advection", exp9 w-sweep and plain-projector disjointness | exp9_wsweep, exp9_plain_peclet |
| eq:rss_ball, eq:deadzone_cond (H-dz), prop:deadzone_stability, `prop3_valid` | `certificates.rss_ball`, `certificates.sector_modulus` | (flag) | exp5 (deadzone vs mid-rise at b_A in {3,5,7,9}) | exp5_deadzone_z, exp5_deadzone_sweep, exp5_deadzone |
| eq:delta_stoch, eq:T_TTL, cor:stochastic_w2 | `stochastic.corollary`, `stochastic.noise_matrix`, `stochastic.clopper_pearson` | -- (MC in exp7) | exp7 sigma grid, fluctuation criterion, 2000 paths, CP interval | exp7_stochastic, exp7_sigma_grid, exp7_corollary_det, exp7_corollary_stoch |
| numerics (not a paper label) | `simulate.phi_matrices`, `simulate.van_loan_cov` | test_determinism::test_two_processes_identical_sha256 | exp8 | exp8_convergence, exp8_numerics |
| main-paper budget split (Fig. 2 left): tail, LUT rounding g_h w_2 (eps_A - eps_CROWN), CROWN g_h eps_CROWN, lag, hold, adversary of eq:delta_total_exact | `certificates.delta_total_exact` (terms), `experiments/exp1_constants.py` | -- | exp1 "the six budget parts sum to delta_total" (R1, R2) | exp1_budget.csv, main_fig2_budget |
| main-paper w-sweep (Fig. 2 right) | `certificates.certificate` | -- | exp9 (delta_total, delta_cert per w) | exp9_wsweep_deltas.csv, main_fig2_budget |
| main-paper ablation panel (Fig. 1) | `plotting.main_fig1_ablation` | test_figures::test_main_figures_rendered_clean | exp6 (a), (b), (c), (e); exp4 (deployed R2, boundary start); exp1 (theorem_valid) | main_fig1_ablation |
| main-paper prose numbers (util_R1, util_R2 = delta_total/delta_cert; e6_c_delta) | `tools/build_section.py` (keys), `tools/build_main.py` (any template) | test_build_main (5 tests) | -- | results/section_numbers.json |
| presentation (V1-V6) | `pgcbf.texnames`, `tools/texlint.py`, `tools/build_section.py`, `tools/standalone_pdf.py`, `tools/redteam_gates.py`, `tools/page_sweep.py` | test_tex_lint (all, incl. test_every_prose_numeral_occurs_in_section_numbers), test_figures (all, incl. size caps and overlapping-axes checker) | run_all.sh steps after the experiments | results/experimental_section.tex, results/section_v_standalone.pdf |

## Ablation map

| ablation | where | what is asserted |
|---|---|---|
| (a) no tail bookkeeping, K_m = 1 | exp6, exp9 | full certificate refuses ((H3), (H4), (H6), (H7)); naive certificate accepts; seeds exit C |
| (b) coarse b_B = 2, no CROWN margin | exp6, exp9 | full certificate refuses ((H6)); naive accepts; seeds exit C |
| (c) single rate tau_A = tau_B = T_s | exp6, exp9, main_fig1_ablation | at deployed T_s single rate fails (H6) while multi-rate passes; at T_s = 0.2 single rate exits, multi-rate does not |
| (d) no deadzone (mid-rise) | exp5 | at every b_A in {3, 5, 7, 9}, including the deployed 9, the mid-rise quantiser chatters while the deadzone gives u_A = 0 late; b_A = 1 violates (H-dz) |
| (e) adversary above eps_adv | exp6, exp9 | eps_adv = 1 refused and 16/16 seeds exit; the largest certifiable eps* has 0 exits |
| (f) noise above the Corollary bound | exp7 | full Corollary refuses; crude bound (tail diffusion ignored) accepts; MC exit frequency far above alpha (CP lower bound > alpha) |
