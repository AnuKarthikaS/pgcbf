# pgcbf: companion code for "One Budget for Every Gap: Certifying Reduced-Order, Quantized, Sampled-Data Neural Barrier Controllers on the Full Plant"

Anu Karthika (Independent Researcher, Chennai, India), iamanukarthika@gmail.com. Submitted to L4DC 2027.
The submitted paper (`main.pdf`) and the supplementary material with all proofs and the full experimental
section (`supplement.pdf`) are attached to the release `v1.0-l4dc-submission` of this repository
(https://github.com/AnuKarthikaS/pgcbf/releases).

Implementation built from the theory sources in `paper/theory/` (the statements and proofs that appear in the
supplement). It builds a class-(H) reaction-diffusion-advection lattice benchmark with an active advective
ROM-tail coupling, computes every certificate quantity at run time from the operators, and checks the paper's
claims on it:

- Theorem thm:main_cbf (hypotheses (H0)-(H7), flag `theorem_valid`, true in R1 and R2); the main paper's Theorem 1;
- Propositions prop:twin_crown, prop:dahlquist and prop:deadzone_stability ((H-dz), flag `prop3_valid`);
- Lemmas quad_split, h2a_char, delta_cert_closed and the class-(H) Lemma lem:delta_cert_classH;
- Corollary cor:stochastic_w2;
- the realised fraction of every budget term along the certified trajectories (Experiment 4b).

Ablations (a)-(f) each produce a certificate refusal and a safety violation. The main-paper figure (Figure 1, two panels) and
Table 1 are rendered from the same CSV tables and logs (`tools/main_figures.py`, `tools/build_main.py`).

## How to run

    pip install -r requirements.txt        # numpy, scipy, matplotlib, pytest (+ pymupdf for the PDF gate)
    ./run_all.sh                           # everything below; exits nonzero on any failure

`run_all.sh`:

1. pins BLAS to 2 threads and the process to cores 0,1 (`taskset`), deletes all old tables, figures
   and logs;
2. runs `experiments/exp0..exp9` (per-script logs and JSON ledgers in `results/logs/`, wall-clock and
   peak RSS via `os.wait4` in `results/logs/runtimes.csv`), then `pytest tests`;
3. `tools/summarize.py`: `results/summary.md`, `results/tables/assertions.csv`, runtime table;
4. `tools/build_section.py`: fills `results/experimental_section_template.tex` from the CSV tables,
   ledgers and config into `results/experimental_section.tex` and `results/section_numbers.json`
   (no number in the prose is typed by hand; fails on an unresolved placeholder or a lint hit);
5. final lint (`tests/test_tex_lint.py`, `tests/test_figures.py`) on the generated tables, Section V
   and figures;
6. `tools/redteam_gates.py` (hypothesis label sets, no `(H8)`, `(H5)` never the deadzone, lemma label,
   clean-room scan);
7. `tools/standalone_pdf.py`: compiles Section V alone in two layouts, IEEEtran one-column (the
   combined-manuscript layout, `results/section_v_standalone.pdf`) and two-column
   (`results/section_v_standalone_twocolumn.pdf`), runs `tools/page_sweep.py` on both, renders every page to
   `results/pdf_pages/{onecolumn,twocolumn}/`, and fails on any page-sweep issue, compile error, undefined
   reference, `??`, overfull box wider than 1 pt, overlapping text, objects in the bottom margin, or a
   figure above the caps (width <= 0.62 textwidth, height <= 0.72 textheight, read from the `\pgcfig` log lines). Set `PDF_GATE=skip` if no TeX engine is
   installed.

To run one experiment: `cd experiments && OPENBLAS_NUM_THREADS=2 python3 exp4_theorem1.py`.

## Layout

    config/default.json        R2 benchmark (Pe = 20, open-loop unstable): all parameters, seeds, tolerances, grids
    config/r1_passive.json     R1 benchmark (Pe = 1): (H-dz) and Corollary regime
    src/pgcbf/lattice.py       class-(H) lattice: pure damped-wave A_full, A_adv with graded filter W_f, F_sin, h
    src/pgcbf/rom.py           Galerkin ROM, P-orthogonal tail projector, A_adv,r perp, (H2a)/(H2b)
    src/pgcbf/controller.py    GroupSort Teacher, b_B-bit Student (Adam learned rounding), port-space Path-A LUT with deadzone
    src/pgcbf/crown.py         twin-CROWN and E_q
    src/pgcbf/certificates.py  one function per certificate quantity, hypothesis ledger, theorem_valid / prop3_valid, refusal
    src/pgcbf/simulate.py      exponential sampled-data simulator, exact latch instants, worst-case adversary, exact noise
    src/pgcbf/stochastic.py    Corollary cor:stochastic_w2, Clopper-Pearson
    src/pgcbf/texnames.py      LaTeX names of every quantity / hypothesis key, booleans -> check mark / cross, refs
    src/pgcbf/plotting.py      raw-trajectory plotting, IEEE widths, constrained layout, artist bounding-box check
    experiments/exp0..exp9     experiments (below)
    tests/                     no-literal AST + perturbation, identities, soundness, refusal, determinism, tex lint, figures
    tools/                     timed_run, summarize, build_section, texlint, redteam_gates, standalone_pdf, page_sweep
    paper/theory/              the theory sources used
    results/                   tables, figures, logs, pdf_pages, experimental_section(.tex, _template.tex),
                               section_v_standalone.pdf, design_note.md, traceability.md, defect_closure.md, summary.md

## What each experiment proves

| script | claim tested |
|---|---|
| exp0_openloop | with Pe = 20 every (H7) initial condition leaves C open loop (16/16) |
| exp1_constants | all constants (Tables I-II) incl. ||A_adv,r perp||_2, L_perp_P split into linear and sine parts, tail share; Q = gamma P to 1e-12; A_r perp = 0 by invariance; D_v closed form = sampled sup; hypothesis ledger, theorem_valid in R1 and R2; refusals: state-space LUT, u_max below (H-sat), unfiltered advection, plain projector, previous T_s = 5e-4, T_s x 10 |
| exp2_crown | Prop. 1: twin-CROWN bound >= sampled gap for b_B = 3..14, both roundings, port box |
| exp3_delta_cert | Lemma lem:delta_cert_classH bound <= brute force (R1, R2), c_M, lambda_min(M,P), s_nl; class-(D) closed form |
| exp4_theorem1 | Theorem 1 on 16 seeds x {R2, R1}, worst-case adversary, boundary start and interior start (h_r = h_act/2): min h, time of min, ratio to Delta h_dip, tail <= eps_basis_P |
| exp4b_budget | each term of the budget (7) evaluated at every integrator node of the certified R2/R1 runs next to the bound the proof uses: realised fractions (adversary, quantization/lag, tail, hold), the two conservative constants of the hold term, the slack of the realised rate over the proof's lower bound |
| exp5_deadzone | Prop. 3 in R1; ablation (d) over b_A in {3, 5, 7, 9} (9 deployed): mid-rise chatters at every b_A, deadzone silent; b_A = 1 violates (H-dz) |
| exp6_ablations | ablations (a) K_m = 1, (b) b_B = 2, (c) single rate, (e) eps_adv = 1: refusal + violation |
| exp7_stochastic | Corollary 1: sigma grid search, fluctuation criterion, 2000 paths, CP interval vs alpha; ablation (f) |
| exp8_numerics | Euler instability, ETD orders, exact latch timing, modal = nodal, SHA-256 determinism |
| exp9_wsweep | w-sweep 1 -> 0.0039 (17 values): accept/refuse, failing hypothesis, eps_basis_P, tail share, ||A_adv,r perp||_2, margin, ablation refusals, open-loop exits; frontier w*; plain-projector Pe scan (disjointness) |

Full label map: `results/traceability.md`. Defect-by-defect changes: `results/defect_closure.md`.

## No hardcoded constants

The only literals are the inputs in `config/*.json` (including the cutoff weight
`plant.adv_cutoff_weight`, the sigma grid and the w-sweep thresholds). `tests/test_no_literals.py`
fails on any float literal in a numerical module, on float literals in experiments outside figure
cosmetics, and if any certificate quantity is insensitive to every config perturbation.

## Tuning record and trade-offs

Details in `results/design_note.md` §4-5.

1. **S1 tail coupling.** Advection band K_f = 4 > K_m = 2 with the graded filter
   W_f = Phi_f diag(1, 1, w, w) Phi_f^T. w is a disclosed benchmark parameter selected by exp9 as the
   largest grid value whose full certificate accepts with a positive (H6) margin while every ablation
   is still refused: w* = 2^{-11/2} = 0.0221.
2. **Structural finding.** Under the Euclidean-ODD supremum defining M_perp, full-strength cutoff-mode
   advection (w = 1) is not certifiable at any Pe where the open loop is unstable (certifiable c <= 0.02,
   unstable c >= 0.5). This is a limitation of the certificate (Remark rem:rda_cutoff), not of the plant.
3. **T_s.** With the tail coupling active, the earlier T_s = 5e-4 is refused ((H6): delta_total 0.01726 >
   delta_cert 0.01455). Deployed R2 T_s = 3e-4 (tau_A 1.5e-5, tau_B 3e-5); Pe stays 20.
4. **(H5) codesign window.** State-space LUT window empty; port-space LUT (default) with b_A = 9 in [8, 11].
5. **(H-sat).** u_max just above the left-hand side (3.386 R2, 0.3386 R1).
6. **Stochastic case.** z0 level 0.9, sigma* = 0.01 = largest accepted grid value; R1 T_s = 2e-3.
7. **Ablation (b)** uses b_B = 2 (b_B = 3 with learned rounding gave no violation).
8. **(H-dz) fails in R2** (L_sec 2.39 > mu_R 0.5); it is not a Theorem-1 hypothesis; Prop. 3 is shown in R1.
9. **Stricter refusal rule**: the codesign lower bound is required.

## Measured runtimes (final run, 2 cores)

| script | exit | wall-clock (s) | peak RSS (MB) |
|---|---|---|---|
| exp0_openloop | 0 | 18.0 | 213.9 |
| exp1_constants | 0 | 4.3 | 111.6 |
| exp2_crown | 0 | 7.2 | 153.1 |
| exp3_delta_cert | 0 | 5.9 | 205.5 |
| exp4_theorem1 | 0 | 52.0 | 342.9 |
| exp4b_budget | 0 | 213.4 | 1009.7 |
| exp5_deadzone | 0 | 53.3 | 209.1 |
| exp6_ablations | 0 | 66.4 | 296.8 |
| exp7_stochastic | 0 | 179.5 | 1050.5 |
| exp8_numerics | 0 | 9.1 | 219.5 |
| exp9_wsweep | 0 | 32.1 | 194.1 |
| pytest | 0 | 24.8 | 1192.9 |
| total suite | 0 | 666.4 | 1193 (max) |

Section build + main-paper figures + lint + red-team gates + standalone PDF compile/page sweep (one- and two-column): 11.4 s (tectonic, cached bundle).
Result (final run): 248/248 experiment assertions; pytest 65 passed + 3 skipped in the main pass (the Section V lint tests skip until the section is built) and 36/36 in the final lint pass; red-team hits 0; page sweep 0 issues on 14 one-column and 13 two-column pages; status PASS. Full log: `results/logs/run_all.log`.

Hardware: cpu: Intel(R) Xeon(R) Processor @ 2.90GHz; cores on machine: 4; used: 2 (BLAS threads=2, taskset -c 0,1); memory: 7957 MB; python 3.14.8 numpy 2.5.4 scipy 1.18.1 matplotlib 3.11.2

## Honest limitations

- exp4: from both the boundary and the interior start the minimum of h is attained at t = 0. The
  certified floor Delta h_dip is attained only from the boundary, by construction of the initial state;
  from the interior, min h exceeds it by a factor of about 5.0e3 (R2) and 4.3e3 (R1). h increases from
  t = 0+ because delta_cert > delta_total, so these runs do not exercise the inter-sample dip dynamically.
- exp7: the Corollary is very conservative: T_TTL = 2022 against a horizon T = 2, and 0/2000 exits
  against alpha = 0.05. The fluctuation criterion is met mostly by the closed-loop drift; the
  noise-induced deviation alone is 0.174 h(x0). On R2 the Corollary is infeasible even at sigma = 0
  (gamma_0 h_act = 0.007237 > margin 9.994e-4); all stochastic runs use r_W2 = 0.
- A correction of a factor 2 in Lip(L_f h) during review: no Theorem 1 verdict changed; exp7 sigma = 0.01 is
  now refused (sigma* 0.01 -> 0.007). See `results/defect_closure.md`.
- exp9: below w* the single-rate ablation (c) is no longer refused, which is what bounds w* from
  below; above w*, (H6) fails. The frontier is a grid point (factor sqrt 2 resolution).
- Ablation (d): R1 is open-loop stable; without the deadzone the limit cycle appears in the actuator
  signal, not as state growth.
- The class-(H) delta_cert is a sound bound (Lemma lem:delta_cert_classH), not the class-(D) closed form,
  which is verified separately on a class-(D) instance.
- The PDF gate needs a TeX engine (tectonic) and pymupdf; these are tooling, not experiment dependencies.
