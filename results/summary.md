# Final run summary

Assertions: 248/248 passed in experiments; script failures: none.

Total wall-clock: 642.9 s (budget 1200 s); peak RSS of any script: 1193 MB.

| script | exit | wall [s] | peak RSS [MB] |
|---|---|---|---|
| exp0_openloop | 0 | 18.40 | 213.7 |
| exp1_constants | 0 | 4.75 | 111.7 |
| exp2_crown | 0 | 6.88 | 149.8 |
| exp3_delta_cert | 0 | 5.91 | 205.5 |
| exp4_theorem1 | 0 | 52.60 | 342.8 |
| exp4b_budget | 0 | 203.64 | 1010.0 |
| exp5_deadzone | 0 | 54.38 | 209.5 |
| exp6_ablations | 0 | 65.25 | 298.6 |
| exp7_stochastic | 0 | 166.43 | 1039.7 |
| exp8_numerics | 0 | 9.20 | 219.8 |
| exp9_wsweep | 0 | 31.42 | 194.0 |
| pytest | 0 | 23.74 | 1193.0 |

| experiment | assertion | result |
|---|---|---|
| exp0_openloop | open-loop linearisation is exponentially unstable (max Re eig > 0) | PASS |
| exp0_openloop | initial states lie in C (h(x0) >= 0) | PASS |
| exp0_openloop | open loop exits C for every seed within the horizon | PASS |
| exp1_constants | R2: Prop2 mu_2P(A_r0) = -gamma/2 | PASS |
| exp1_constants | R2: Prop2 Q = gamma P identity | PASS |
| exp1_constants | R2: Prop2 full-order mu_2,Pfe(A_full) = -gamma/2 | PASS |
| exp1_constants | R2: tail rate lambda_K+1 equals modal pencil rate | PASS |
| exp1_constants | R2: tail energy decays exactly as exp(-gamma t/2) (Q = gamma P) | PASS |
| exp1_constants | R2: Lemma quad_split h(Phi z + x_perp) = h(Phi z) - //x_perp//^2 (rel 1e-12) | PASS |
| exp1_constants | R2: (H2a) Phi^T P_eps P_perp = 0 | PASS |
| exp1_constants | R2: Q = gamma P relative error < 1e-12 with K_f = 4 > K_m = 2 (A_full advection-free) | PASS |
| exp1_constants | R2: A_r,perp (A_full-based) vanishes by eigenmode invariance | PASS |
| exp1_constants | R2: S3 D_v closed form = attained sup (round-off) >= sampled sup (within 1%) | PASS |
| exp1_constants | R2: hypothesis H0 | PASS |
| exp1_constants | R2: hypothesis H1 | PASS |
| exp1_constants | R2: hypothesis H2 | PASS |
| exp1_constants | R2: hypothesis H3 | PASS |
| exp1_constants | R2: hypothesis H4 | PASS |
| exp1_constants | R2: hypothesis H5 | PASS |
| exp1_constants | R2: hypothesis H6 | PASS |
| exp1_constants | R2: hypothesis H7 | PASS |
| exp1_constants | R2: hypothesis H_sat | PASS |
| exp1_constants | R2: hypothesis timing_ok | PASS |
| exp1_constants | R2: hypothesis codesign_window | PASS |
| exp1_constants | R2: theorem_valid (all of (H0)-(H7)) and the certificate does not refuse | PASS |
| exp1_constants | R1: Prop2 mu_2P(A_r0) = -gamma/2 | PASS |
| exp1_constants | R1: Prop2 Q = gamma P identity | PASS |
| exp1_constants | R1: Prop2 full-order mu_2,Pfe(A_full) = -gamma/2 | PASS |
| exp1_constants | R1: tail rate lambda_K+1 equals modal pencil rate | PASS |
| exp1_constants | R1: tail energy decays exactly as exp(-gamma t/2) (Q = gamma P) | PASS |
| exp1_constants | R1: Lemma quad_split h(Phi z + x_perp) = h(Phi z) - //x_perp//^2 (rel 1e-12) | PASS |
| exp1_constants | R1: (H2a) Phi^T P_eps P_perp = 0 | PASS |
| exp1_constants | R1: Q = gamma P relative error < 1e-12 with K_f = 4 > K_m = 2 (A_full advection-free) | PASS |
| exp1_constants | R1: A_r,perp (A_full-based) vanishes by eigenmode invariance | PASS |
| exp1_constants | R1: S3 D_v closed form = attained sup (round-off) >= sampled sup (within 1%) | PASS |
| exp1_constants | R1: hypothesis H0 | PASS |
| exp1_constants | R1: hypothesis H1 | PASS |
| exp1_constants | R1: hypothesis H2 | PASS |
| exp1_constants | R1: hypothesis H3 | PASS |
| exp1_constants | R1: hypothesis H4 | PASS |
| exp1_constants | R1: hypothesis H5 | PASS |
| exp1_constants | R1: hypothesis H6 | PASS |
| exp1_constants | R1: hypothesis H7 | PASS |
| exp1_constants | R1: hypothesis H_sat | PASS |
| exp1_constants | R1: hypothesis timing_ok | PASS |
| exp1_constants | R1: hypothesis codesign_window | PASS |
| exp1_constants | R1: theorem_valid (all of (H0)-(H7)) and the certificate does not refuse | PASS |
| exp1_constants | R1: prop3_valid ((H-dz): r_ss < r_min and mu_R > L_sec) | PASS |
| exp1_constants | R2: (H-dz) honestly fails (L_sec > mu_R: open-loop unstable advection); Prop. 3 not claimed | PASS |
| exp1_constants | G-D: //A_adv,r,perp//_2 > 0 in R2 (K_f > K_m) | PASS |
| exp1_constants | G-D: L_perp_P(R2) materially larger than L_perp_P(R1) (ratio >= Pe ratio / 2) | PASS |
| exp1_constants | G-D: L_perp_P = linear advective part + sin-nonlinearity part (both reported) | PASS |
| exp1_constants | G-D: tail-coupling share of delta_total reported and material in R2 (>= 10%) | PASS |
| exp1_constants | certificate refuses: state-space LUT ($m = K$ variant) | PASS |
| exp1_constants | certificate refuses: $u_{\max}$ below the (H-sat) left side | PASS |
| exp1_constants | certificate refuses: pure lattice advection (no spectral cutoff) | PASS |
| exp1_constants | certificate refuses: plain projector $\Pi_f$ ($w = 1$) | PASS |
| exp1_constants | certificate refuses: previous sampling period $T_s=5\times10^{-4}$ at $w^\star$ | PASS |
| exp1_constants | certificate refuses: $T_s$, $\tau_A$, $\tau_B$ $\times 10$ (hold term) | PASS |
| exp1_constants | R1: Prop. 3(i) contraction mu_R > L_sec holds | PASS |
| exp1_constants | R2: Prop. 3(i) fails (L_sec > mu_R), so Proposition 3 does not apply on R2 | PASS |
| exp1_constants | R2: gamma0 h_act exceeds the (H6) margin delta_cert - delta_total (Corollary 1 infeasible on R2 even at sigma = 0) | PASS |
| exp1_constants | R1: the six budget parts (tail, LUT, CROWN, lag, hold, adversary) sum to delta_total | PASS |
| exp1_constants | R2: the six budget parts (tail, LUT, CROWN, lag, hold, adversary) sum to delta_total | PASS |
| exp2_crown | b_B=3 nearest: twin-CROWN bound >= max sampled /Psi_B - T/ (100016 samples) | PASS |
| exp2_crown | b_B=3 learned: twin-CROWN bound >= max sampled /Psi_B - T/ (100016 samples) | PASS |
| exp2_crown | b_B=4 nearest: twin-CROWN bound >= max sampled /Psi_B - T/ (100016 samples) | PASS |
| exp2_crown | b_B=4 learned: twin-CROWN bound >= max sampled /Psi_B - T/ (100016 samples) | PASS |
| exp2_crown | b_B=5 nearest: twin-CROWN bound >= max sampled /Psi_B - T/ (100016 samples) | PASS |
| exp2_crown | b_B=5 learned: twin-CROWN bound >= max sampled /Psi_B - T/ (100016 samples) | PASS |
| exp2_crown | b_B=6 nearest: twin-CROWN bound >= max sampled /Psi_B - T/ (100016 samples) | PASS |
| exp2_crown | b_B=6 learned: twin-CROWN bound >= max sampled /Psi_B - T/ (100016 samples) | PASS |
| exp2_crown | b_B=7 nearest: twin-CROWN bound >= max sampled /Psi_B - T/ (100016 samples) | PASS |
| exp2_crown | b_B=7 learned: twin-CROWN bound >= max sampled /Psi_B - T/ (100016 samples) | PASS |
| exp2_crown | b_B=8 nearest: twin-CROWN bound >= max sampled /Psi_B - T/ (100016 samples) | PASS |
| exp2_crown | b_B=8 learned: twin-CROWN bound >= max sampled /Psi_B - T/ (100016 samples) | PASS |
| exp2_crown | b_B=9 nearest: twin-CROWN bound >= max sampled /Psi_B - T/ (100016 samples) | PASS |
| exp2_crown | b_B=9 learned: twin-CROWN bound >= max sampled /Psi_B - T/ (100016 samples) | PASS |
| exp2_crown | b_B=10 nearest: twin-CROWN bound >= max sampled /Psi_B - T/ (100016 samples) | PASS |
| exp2_crown | b_B=10 learned: twin-CROWN bound >= max sampled /Psi_B - T/ (100016 samples) | PASS |
| exp2_crown | b_B=11 nearest: twin-CROWN bound >= max sampled /Psi_B - T/ (100016 samples) | PASS |
| exp2_crown | b_B=11 learned: twin-CROWN bound >= max sampled /Psi_B - T/ (100016 samples) | PASS |
| exp2_crown | b_B=12 nearest: twin-CROWN bound >= max sampled /Psi_B - T/ (100016 samples) | PASS |
| exp2_crown | b_B=12 learned: twin-CROWN bound >= max sampled /Psi_B - T/ (100016 samples) | PASS |
| exp2_crown | b_B=13 nearest: twin-CROWN bound >= max sampled /Psi_B - T/ (100016 samples) | PASS |
| exp2_crown | b_B=13 learned: twin-CROWN bound >= max sampled /Psi_B - T/ (100016 samples) | PASS |
| exp2_crown | b_B=14 nearest: twin-CROWN bound >= max sampled /Psi_B - T/ (100016 samples) | PASS |
| exp2_crown | b_B=14 learned: twin-CROWN bound >= max sampled /Psi_B - T/ (100016 samples) | PASS |
| exp2_crown | nearest: eps_CROWN monotone non-increasing in b_B | PASS |
| exp2_crown | nearest: eps_CROWN scales as 2^(-b_B) (log2 slope within [-1.25,-0.75]) | PASS |
| exp2_crown | learned: eps_CROWN monotone non-increasing in b_B | PASS |
| exp2_crown | learned: eps_CROWN scales as 2^(-b_B) (log2 slope within [-1.25,-0.75]) | PASS |
| exp2_crown | port box: twin-CROWN bound >= sampled gap on the Path-A table domain | PASS |
| exp2_crown | learned rounding does not increase the certified eps_CROWN at the deployed b_B | PASS |
| exp3_delta_cert | R2: sound two-branch delta_cert <= brute-force inf over C_act | PASS |
| exp3_delta_cert | R2: (H6) delta_total <= delta_cert <= brute-force inf | PASS |
| exp3_delta_cert | R2: Lemma delta_cert_classH: c_M = lambda_min(M,P) - 2 s_nl | PASS |
| exp3_delta_cert | R1: sound two-branch delta_cert <= brute-force inf over C_act | PASS |
| exp3_delta_cert | R1: (H6) delta_total <= delta_cert <= brute-force inf | PASS |
| exp3_delta_cert | R1: Lemma delta_cert_classH: c_M = lambda_min(M,P) - 2 s_nl | PASS |
| exp3_delta_cert | class D gamma0=0.3: closed form = brute-force inf of affine bound over rho | PASS |
| exp3_delta_cert | class D gamma0=0.3: closed form <= brute-force inf of exact bracket | PASS |
| exp3_delta_cert | class D gamma0=0.3: branch 1 iff c >= gamma0 | PASS |
| exp3_delta_cert | class D gamma0=3.0: closed form = brute-force inf of affine bound over rho | PASS |
| exp3_delta_cert | class D gamma0=3.0: closed form <= brute-force inf of exact bracket | PASS |
| exp3_delta_cert | class D gamma0=3.0: branch 1 iff c >= gamma0 | PASS |
| exp4_theorem1 | R2: interior level h_r = 0.5 h_act lies strictly inside C~_r | PASS |
| exp4_theorem1 | R2 boundary: (a) tail //x_perp(t)//_Pfe <= eps_basis_P (all seeds, all nodes) | PASS |
| exp4_theorem1 | R2 boundary: (b) h~_r(z(t_k)) >= 0 at every sampling instant | PASS |
| exp4_theorem1 | R2 boundary: (b) h~_r(z(t)) >= -Delta h_dip at every node | PASS |
| exp4_theorem1 | R2 boundary: (c) h(x(t)) >= h~_r(z(t)) + Delta h_dip | PASS |
| exp4_theorem1 | R2 boundary: (c) h(x(t)) >= 0 for all t, all seeds (full-order safety) | PASS |
| exp4_theorem1 | R2 boundary: z(t_k) stays in Z_ODD | PASS |
| exp4_theorem1 | R2 boundary: min_t h = Delta h_dip at t = 0 (initial state on the certified floor) | PASS |
| exp4_theorem1 | R2 interior: (a) tail //x_perp(t)//_Pfe <= eps_basis_P (all seeds, all nodes) | PASS |
| exp4_theorem1 | R2 interior: (b) h~_r(z(t_k)) >= 0 at every sampling instant | PASS |
| exp4_theorem1 | R2 interior: (b) h~_r(z(t)) >= -Delta h_dip at every node | PASS |
| exp4_theorem1 | R2 interior: (c) h(x(t)) >= h~_r(z(t)) + Delta h_dip | PASS |
| exp4_theorem1 | R2 interior: (c) h(x(t)) >= 0 for all t, all seeds (full-order safety) | PASS |
| exp4_theorem1 | R2 interior: z(t_k) stays in Z_ODD | PASS |
| exp4_theorem1 | R2 interior: min_t h(x_t) = min_j h(x_j(0)) (attained at t = 0; h rises along every closed-loop trajectory) | PASS |
| exp4_theorem1 | R1: interior level h_r = 0.5 h_act lies strictly inside C~_r | PASS |
| exp4_theorem1 | R1 boundary: (a) tail //x_perp(t)//_Pfe <= eps_basis_P (all seeds, all nodes) | PASS |
| exp4_theorem1 | R1 boundary: (b) h~_r(z(t_k)) >= 0 at every sampling instant | PASS |
| exp4_theorem1 | R1 boundary: (b) h~_r(z(t)) >= -Delta h_dip at every node | PASS |
| exp4_theorem1 | R1 boundary: (c) h(x(t)) >= h~_r(z(t)) + Delta h_dip | PASS |
| exp4_theorem1 | R1 boundary: (c) h(x(t)) >= 0 for all t, all seeds (full-order safety) | PASS |
| exp4_theorem1 | R1 boundary: z(t_k) stays in Z_ODD | PASS |
| exp4_theorem1 | R1 boundary: min_t h = Delta h_dip at t = 0 (initial state on the certified floor) | PASS |
| exp4_theorem1 | R1 interior: (a) tail //x_perp(t)//_Pfe <= eps_basis_P (all seeds, all nodes) | PASS |
| exp4_theorem1 | R1 interior: (b) h~_r(z(t_k)) >= 0 at every sampling instant | PASS |
| exp4_theorem1 | R1 interior: (b) h~_r(z(t)) >= -Delta h_dip at every node | PASS |
| exp4_theorem1 | R1 interior: (c) h(x(t)) >= h~_r(z(t)) + Delta h_dip | PASS |
| exp4_theorem1 | R1 interior: (c) h(x(t)) >= 0 for all t, all seeds (full-order safety) | PASS |
| exp4_theorem1 | R1 interior: z(t_k) stays in Z_ODD | PASS |
| exp4_theorem1 | R1 interior: min_t h(x_t) = min_j h(x_j(0)) (attained at t = 0; h rises along every closed-loop trajectory) | PASS |
| exp4b_budget | R2 boundary: node set = every integrator step of every period and seed | PASS |
| exp4b_budget | R2 boundary: decomposition zdot = f_r + g_r u + d_r w + Delta_perp closes to roundoff (every node) | PASS |
| exp4b_budget | R2 boundary: identity r = Theta(z_k) + H + Q + A + C closes to roundoff (every node) | PASS |
| exp4b_budget | R2 boundary: (i) tail coupling realised/bound <= 1 at every node | PASS |
| exp4b_budget | R2 boundary: (ii) adversary realised/bound <= 1 at every node | PASS |
| exp4b_budget | R2 boundary: (iii) quantization/lag realised/bound <= 1 at every node | PASS |
| exp4b_budget | R2 boundary: (iv) hold shift realised/bound <= 1 at every node | PASS |
| exp4b_budget | R2 boundary: factor bounds //2Pz// <= L_h, //Delta_perp// <= L_perp //x_perp//, e_u <= ebar_u, //L_gr h~_r(z_k)// <= g_h, //w// <= eps_adv at every node | PASS |
| exp4b_budget | R2 boundary: //z - z_k// <= M_z s at every node | PASS |
| exp4b_budget | R2 boundary: (v) soundness: realised rate - RHS(s) of eq:ode_diff_ineq >= 0 at every active-shell node | PASS |
| exp4b_budget | R2 boundary: Teacher condition Theta(z_k) >= delta_cert at every active-shell sample | PASS |
| exp4b_budget | R2 boundary: interior branch h~_r(z(t)) >= h~_r(z_k) - L_h M_z s at every B_int node | PASS |
| exp4b_budget | R2 interior: node set = every integrator step of every period and seed | PASS |
| exp4b_budget | R2 interior: decomposition zdot = f_r + g_r u + d_r w + Delta_perp closes to roundoff (every node) | PASS |
| exp4b_budget | R2 interior: identity r = Theta(z_k) + H + Q + A + C closes to roundoff (every node) | PASS |
| exp4b_budget | R2 interior: (i) tail coupling realised/bound <= 1 at every node | PASS |
| exp4b_budget | R2 interior: (ii) adversary realised/bound <= 1 at every node | PASS |
| exp4b_budget | R2 interior: (iii) quantization/lag realised/bound <= 1 at every node | PASS |
| exp4b_budget | R2 interior: (iv) hold shift realised/bound <= 1 at every node | PASS |
| exp4b_budget | R2 interior: factor bounds //2Pz// <= L_h, //Delta_perp// <= L_perp //x_perp//, e_u <= ebar_u, //L_gr h~_r(z_k)// <= g_h, //w// <= eps_adv at every node | PASS |
| exp4b_budget | R2 interior: //z - z_k// <= M_z s at every node | PASS |
| exp4b_budget | R2 interior: (v) soundness: realised rate - RHS(s) of eq:ode_diff_ineq >= 0 at every active-shell node | PASS |
| exp4b_budget | R2 interior: Teacher condition Theta(z_k) >= delta_cert at every active-shell sample | PASS |
| exp4b_budget | R2 interior: interior branch h~_r(z(t)) >= h~_r(z_k) - L_h M_z s at every B_int node | PASS |
| exp4b_budget | R1 boundary: node set = every integrator step of every period and seed | PASS |
| exp4b_budget | R1 boundary: decomposition zdot = f_r + g_r u + d_r w + Delta_perp closes to roundoff (every node) | PASS |
| exp4b_budget | R1 boundary: identity r = Theta(z_k) + H + Q + A + C closes to roundoff (every node) | PASS |
| exp4b_budget | R1 boundary: (i) tail coupling realised/bound <= 1 at every node | PASS |
| exp4b_budget | R1 boundary: (ii) adversary realised/bound <= 1 at every node | PASS |
| exp4b_budget | R1 boundary: (iii) quantization/lag realised/bound <= 1 at every node | PASS |
| exp4b_budget | R1 boundary: (iv) hold shift realised/bound <= 1 at every node | PASS |
| exp4b_budget | R1 boundary: factor bounds //2Pz// <= L_h, //Delta_perp// <= L_perp //x_perp//, e_u <= ebar_u, //L_gr h~_r(z_k)// <= g_h, //w// <= eps_adv at every node | PASS |
| exp4b_budget | R1 boundary: //z - z_k// <= M_z s at every node | PASS |
| exp4b_budget | R1 boundary: (v) soundness: realised rate - RHS(s) of eq:ode_diff_ineq >= 0 at every active-shell node | PASS |
| exp4b_budget | R1 boundary: Teacher condition Theta(z_k) >= delta_cert at every active-shell sample | PASS |
| exp4b_budget | R1 boundary: interior branch h~_r(z(t)) >= h~_r(z_k) - L_h M_z s at every B_int node | PASS |
| exp4b_budget | R1 interior: node set = every integrator step of every period and seed | PASS |
| exp4b_budget | R1 interior: decomposition zdot = f_r + g_r u + d_r w + Delta_perp closes to roundoff (every node) | PASS |
| exp4b_budget | R1 interior: identity r = Theta(z_k) + H + Q + A + C closes to roundoff (every node) | PASS |
| exp4b_budget | R1 interior: (i) tail coupling realised/bound <= 1 at every node | PASS |
| exp4b_budget | R1 interior: (ii) adversary realised/bound <= 1 at every node | PASS |
| exp4b_budget | R1 interior: (iii) quantization/lag realised/bound <= 1 at every node | PASS |
| exp4b_budget | R1 interior: (iv) hold shift realised/bound <= 1 at every node | PASS |
| exp4b_budget | R1 interior: factor bounds //2Pz// <= L_h, //Delta_perp// <= L_perp //x_perp//, e_u <= ebar_u, //L_gr h~_r(z_k)// <= g_h, //w// <= eps_adv at every node | PASS |
| exp4b_budget | R1 interior: //z - z_k// <= M_z s at every node | PASS |
| exp4b_budget | R1 interior: (v) soundness: realised rate - RHS(s) of eq:ode_diff_ineq >= 0 at every active-shell node | PASS |
| exp4b_budget | R1 interior: Teacher condition Theta(z_k) >= delta_cert at every active-shell sample | PASS |
| exp4b_budget | R1 interior: interior branch h~_r(z(t)) >= h~_r(z_k) - L_h M_z s at every B_int node | PASS |
| exp5_deadzone | full method (deadzone, b_A=9): prop3_valid ((H-dz) r_ss < r_min and mu_R > L_sec) | PASS |
| exp5_deadzone | Prop3(b): //z(t_k)//_P <= max(//z(0)//_P, r_ss) for all k, all seeds | PASS |
| exp5_deadzone | Prop3: trajectories enter B_ss(r_ss) subset B_int | PASS |
| exp5_deadzone | deadzone: Path-A control exactly zero in the final window (no limit cycle) | PASS |
| exp5_deadzone | (d) b_A=3: deadzone -> u_A exactly 0 late, no switching, zero late effort | PASS |
| exp5_deadzone | (d) b_A=3: mid-rise (no deadzone) -> u_A never zero, persistent sign switching (chattering), effort > 0 | PASS |
| exp5_deadzone | (d) b_A=3: deadzone -> once u = 0, //z//_P contracts at least at rate mu_R - L_sec (Prop3(a)) | PASS |
| exp5_deadzone | (d) b_A=3: mid-rise lag error /u_A - T(z)/ persists at the LSB scale | PASS |
| exp5_deadzone | (d) b_A=5: deadzone -> u_A exactly 0 late, no switching, zero late effort | PASS |
| exp5_deadzone | (d) b_A=5: mid-rise (no deadzone) -> u_A never zero, persistent sign switching (chattering), effort > 0 | PASS |
| exp5_deadzone | (d) b_A=5: deadzone -> once u = 0, //z//_P contracts at least at rate mu_R - L_sec (Prop3(a)) | PASS |
| exp5_deadzone | (d) b_A=5: mid-rise lag error /u_A - T(z)/ persists at the LSB scale | PASS |
| exp5_deadzone | (d) b_A=7: deadzone -> u_A exactly 0 late, no switching, zero late effort | PASS |
| exp5_deadzone | (d) b_A=7: mid-rise (no deadzone) -> u_A never zero, persistent sign switching (chattering), effort > 0 | PASS |
| exp5_deadzone | (d) b_A=7: deadzone -> once u = 0, //z//_P contracts at least at rate mu_R - L_sec (Prop3(a)) | PASS |
| exp5_deadzone | (d) b_A=7: mid-rise lag error /u_A - T(z)/ persists at the LSB scale | PASS |
| exp5_deadzone | (d) b_A=9 (deployed b_A): deadzone -> u_A exactly 0 late, no switching, zero late effort | PASS |
| exp5_deadzone | (d) b_A=9 (deployed b_A): mid-rise (no deadzone) -> u_A never zero, persistent sign switching (chattering), effort > 0 | PASS |
| exp5_deadzone | (d) b_A=9: deadzone -> once u = 0, //z//_P contracts at least at rate mu_R - L_sec (Prop3(a)) | PASS |
| exp5_deadzone | (d) b_A=9: mid-rise lag error /u_A - T(z)/ persists at the LSB scale | PASS |
| exp5_deadzone | d2: b_A = 1 violates (H-dz): r_ss >= r_min, certificate refuses Prop. 3 (prop3_valid false) | PASS |
| exp6_ablations | (a) full certificate refuses K_m too small (tail criterion) | PASS |
| exp6_ablations | (a) naive certificate without tail terms would pass (H4) and (H6) | PASS |
| exp6_ablations | (a) K_m too small: full-order safety violated (some seed exits C) | PASS |
| exp6_ablations | (b) full certificate refuses coarse b_B | PASS |
| exp6_ablations | (b) naive certificate without the quantization term would accept | PASS |
| exp6_ablations | (b) coarse b_B: full-order safety violated | PASS |
| exp6_ablations | (c) single-rate at the deployed T_s: certificate refuses (lag term w1 = 1) | PASS |
| exp6_ablations | (c) multi-rate hold strictly reduces delta_total and keeps (H6) | PASS |
| exp6_ablations | (c) single-rate at T_s=0.2: full-order safety violated | PASS |
| exp6_ablations | (c) multi-rate at the same T_s=0.2: no exit (outside certificate scope, observed) | PASS |
| exp6_ablations | (e) certificate refuses eps_adv above the certified level | PASS |
| exp6_ablations | (e) adversary eps_adv=1.0: full-order safety violated | PASS |
| exp6_ablations | (e) at the largest certifiable eps_adv* the worst-case adversary does not break safety | PASS |
| exp7_stochastic | z(0) lies strictly inside C~_r and h(x_0) > 0 | PASS |
| exp7_stochastic | sigma grid: some grid sigma is accepted by the Corollary at T and meets the fluctuation criterion | PASS |
| exp7_stochastic | sigma* is the largest grid value with Corollary acceptance and the fluctuation criterion (every larger grid value fails one) | PASS |
| exp7_stochastic | the Corollary's acceptance region on the grid is a down-set (monotone in sigma) | PASS |
| exp7_stochastic | sigma*: Corollary conditions hold (margin, beta0 > h0, h_act^- > 0, T <= T_TTL, Theorem 1) | PASS |
| exp7_stochastic | sigma*: fluctuation criterion holds on the full Monte Carlo (>= fluct_path_share of the paths) | PASS |
| exp7_stochastic | sigma*: the noise itself is non-trivial: median max_t /h(x_t) - h_det(t)/ >= fluct_fraction_of_h0 * h(x_0) | PASS |
| exp7_stochastic | sigma*: MC exit frequency <= alpha (Clopper-Pearson upper limit <= alpha) | PASS |
| exp7_stochastic | sigma*: the probability bound is conservative (T_TTL >> T and zero observed exits) | PASS |
| exp7_stochastic | (f) tail-white noise: full Corollary refuses (OU tail envelope exceeds the safe set) | PASS |
| exp7_stochastic | (f) crude bound ignoring tail diffusion would accept | PASS |
| exp7_stochastic | (f) MC violates the crude probability bound (Clopper-Pearson lower limit > alpha) | PASS |
| exp7_stochastic | R2, sigma = 0: gamma0 h_act exceeds delta_cert - delta_total, so the Corollary margin fails and the Corollary refuses | PASS |
| exp7_stochastic | stochastic runs use r_W2 = 0 (no distributional drift of the adversary) | PASS |
| exp8_numerics | stiffness: at the convergence-test step explicit Euler is unstable (max/1 + h lambda/ > 1) while the exponential integrator is used | PASS |
| exp8_numerics | ETD2RK convergence: observed order in [1.8, 2.2] on step halving (error ratio ~ 4) | PASS |
| exp8_numerics | ETD1 (exponential Euler) for reference: observed order ~ 1 | PASS |
| exp8_numerics | exact latch timing: simulator == closed-form augmented exponential (kappa=0) to round-off | PASS |
| exp8_numerics | step-resolution latching would be O(step) wrong (discriminating test) | PASS |
| exp8_numerics | modal sparse propagation == dense nodal ETD2RK (closed loop) | PASS |
| exp8_numerics | bit-for-bit determinism: two runs give identical SHA-256 | PASS |
| exp9_wsweep | w grid is geometric, contains w = 1 and reaches down to w <= 0.005 with >= 8 values | PASS |
| exp9_wsweep | plain projector w = 1 is refused at the R2 Peclet number | PASS |
| exp9_wsweep | open loop exits C for every seed at every w (the plant is unsafe for every swept w) | PASS |
| exp9_wsweep | //A_adv,r,perp//_2 > 0 at every w (K_f > K_m) and proportional to w | PASS |
| exp9_wsweep | a certifiable frontier exists: some grid w accepts with (H6) margin > 0 and all ablations refused | PASS |
| exp9_wsweep | certifiable set is a down-set of the grid (every w below the frontier accepts) | PASS |
| exp9_wsweep | deployed config weight equals the sweep-selected w* (largest admissible grid value) | PASS |
| exp9_wsweep | at w*, the tail-coupling share of delta_total is material (>= wstar_min_tail_share) | PASS |
| exp9_wsweep | plain projector (w = 1): certifiable-c set and open-loop-unstable-c set are disjoint and both non-empty | PASS |
