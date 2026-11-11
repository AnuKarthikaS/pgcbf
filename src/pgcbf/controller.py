"""Isometric GroupSort-2 Teacher, b_B-bit fixed-point Student (QNN), LUT path.

Implements eq:groupsort_relu, eq:teacher_weights, eq:teacher_def, eq:LTz,
the mid-tread quantiser Q_{b_A} with deadzone (Section sec:problem, item 1) and
the learned-rounding distillation of the Student (numpy Adam).
"""
from __future__ import annotations

import numpy as np


# --------------------------------------------------------------------- helpers
def gs2(Y: np.ndarray) -> np.ndarray:
    """Implements eq:groupsort_relu  pairwise MaxMin on rows (2i, 2i+1) (0-based); Y is H x P."""
    a, b = Y[0::2], Y[1::2]
    r = np.maximum(a - b, 0)
    out = np.empty_like(Y)
    out[0::2] = b + r  # max
    out[1::2] = a - r  # min
    return out


def random_orthogonal(n: int, rng: np.random.Generator) -> np.ndarray:
    A = rng.standard_normal((n, n))
    Q, R = np.linalg.qr(A)
    return Q * np.sign(np.diag(R))


def project_U(A: np.ndarray, u_max: float) -> np.ndarray:
    """Implements eq:teacher_def: Pi_U(a) = a min(1, u_max/||a||) (eq:teacher_def), columns."""
    nrm = np.linalg.norm(A, axis=0)
    scale = np.minimum(1, u_max / np.maximum(nrm, np.finfo(float).tiny))
    return A * scale


def round_grid(X: np.ndarray, step: float) -> np.ndarray:
    return np.round(X / step) * step


# --------------------------------------------------------------------- teacher
class Teacher:
    """Implements eq:teacher_weights / Lemma lem:sector. Origin-preserving L-layer isometric GroupSort-2 network acting on v = g_r^T P z."""

    def __init__(self, m: int, H: int, L: int, theta: float, L_T: float, seed: int):
        assert H % 2 == 0 and H // 2 >= m and L >= 2, "architecture violates H even, H/2>=m, L>=2"
        rng = np.random.default_rng(seed)
        h2 = H // 2
        self.m, self.H, self.L, self.theta, self.L_T = m, H, L, theta, L_T
        Us = np.zeros((H, h2))
        Ud = np.zeros((H, h2))
        for i in range(h2):
            Us[2 * i, i] = Us[2 * i + 1, i] = 1 / np.sqrt(2)
            Ud[2 * i, i] = 1 / np.sqrt(2)
            Ud[2 * i + 1, i] = -1 / np.sqrt(2)
        self.Usum, self.Udiff = Us, Ud
        Win = random_orthogonal(h2, rng)[:, :m]
        Qmix = random_orthogonal(h2, rng)
        W = [np.cos(theta) * Us @ Win + np.sin(theta) * Ud @ Qmix @ Win]
        for _ in range(2, L):
            Qs, Qd = random_orthogonal(h2, rng), random_orthogonal(h2, rng)
            W.append(Us @ Qs @ Us.T + Ud @ Qd @ Ud.T)
        prod = np.eye(m)
        for Wl in W:
            prod = Wl @ prod
        W.append(-L_T * prod.T)
        self.W = W
        self.eta = L_T * np.cos(2 * theta)  # Lemma lem:sector

    def forward(self, V: np.ndarray) -> np.ndarray:
        S = V
        for Wl in self.W[:-1]:
            S = gs2(Wl @ S)
        return self.W[-1] @ S


# --------------------------------------------------------------------- student
class QuantizedNet:
    """Implements Proposition prop:twin_crown (quantized net). b_B-bit fixed-point realisation N_q of the Teacher (Proposition prop:twin_crown).

    Weight grid step  s_l 2^{-b_B}  (s_l = max|W_l|), so floor/ceil/nearest choices
    all satisfy |W^q - W| <= s_l 2^{-b_B}.  Input and every pre-activation are rounded
    (round-half-even) to grids of step a_l 2^{-b_B}; the rounding errors are the
    injections e_l of the Proposition with rho_l = a_l 2^{-b_B} / 2.
    """

    def __init__(self, teacher: Teacher, b_B: int, Y: float, offsets=None):
        self.t = teacher
        self.b_B = int(b_B)
        self.Y = Y
        L = teacher.L
        self.scales = [np.max(np.abs(Wl)) for Wl in teacher.W]
        self.wstep = [s * 2 ** (-float(self.b_B)) for s in self.scales]
        if offsets is None:  # nearest rounding
            self.Wq = [round_grid(Wl, st) for Wl, st in zip(teacher.W, self.wstep)]
        else:  # learned rounding: floor + offset in {0,1}
            self.Wq = [np.floor(Wl / st) * st + st * o for Wl, st, o in zip(teacher.W, self.wstep, offsets)]
        self.dW = [Wq - W for Wq, W in zip(self.Wq, teacher.W)]
        # activation ranges: ||s_l||_2 <= Y for l < L (norm preservation, ||W_l|| <= 1)
        self.act_range = [Y] + [Y] * (L - 1) + [teacher.L_T * Y]  # index 0: input, l: pre-act of layer l
        self.astep = [a * 2 ** (-float(self.b_B)) for a in self.act_range]
        self.rho = [st / 2 for st in self.astep]  # rho_0..rho_L
        self.bq = [np.zeros(Wl.shape[0]) for Wl in teacher.W]  # zero biases => Delta_b = 0
        self.widths = [teacher.m] + [teacher.H] * (L - 1) + [teacher.m]

    def forward(self, V: np.ndarray) -> np.ndarray:
        S = round_grid(V, self.astep[0])
        for l, Wl in enumerate(self.Wq[:-1], start=1):
            S = gs2(round_grid(Wl @ S, self.astep[l]))
        return round_grid(self.Wq[-1] @ S, self.astep[-1])

    def forward_injected(self, V: np.ndarray, E: list) -> np.ndarray:
        """y_l^q = W_l^q s_{l-1}^q + b_l^q + e_l with arbitrary injections (soundness tests)."""
        S = V + E[0]
        for l, Wl in enumerate(self.Wq[:-1], start=1):
            S = gs2(Wl @ S + self.bq[l - 1][:, None] + E[l])
        return self.Wq[-1] @ S + self.bq[-1][:, None] + E[-1]


def learn_rounding(teacher: Teacher, b_B: int, Y: float, V_train: np.ndarray, tcfg: dict, seed: int):
    """Distil the Teacher into the b_B-bit Student by learning floor/ceil offsets with Adam.

    Relaxation o = sigmoid(a / T) with temperature annealed from T0 to T1; loss
    mean ||N_q(v) - N(v)||^2 on the training set (activation rounding is not applied
    during training).  Returns hard offsets in {0,1} and the training history.
    """
    rng = np.random.default_rng(seed)
    base = QuantizedNet(teacher, b_B, Y)
    W, st = teacher.W, base.wstep
    frac = [W_l / s - np.floor(W_l / s) for W_l, s in zip(W, st)]
    # initialise at nearest rounding (logit of the fractional part)
    clipf = float(tcfg["init_clip"])
    a = [np.log(np.clip(f, clipf, 1 - clipf) / (1 - np.clip(f, clipf, 1 - clipf))) for f in frac]
    mA = [np.zeros_like(x) for x in a]
    vA = [np.zeros_like(x) for x in a]
    b1, b2, epsA, lr = (float(tcfg[k]) for k in ("beta1", "beta2", "adam_eps", "lr"))
    T0, T1 = float(tcfg["temp_start"]), float(tcfg["temp_end"])
    n_it, bs = int(tcfg["iterations"]), int(tcfg["batch"])
    target_all = teacher.forward(V_train)
    hist = []
    floors = [np.floor(W_l / s) * s for W_l, s in zip(W, st)]
    L = teacher.L
    for it in range(n_it):
        T = T0 * (T1 / T0) ** (it / max(n_it - 1, 1))
        idx = rng.integers(0, V_train.shape[1], bs)
        V = V_train[:, idx]
        tgt = target_all[:, idx]
        sig = [1 / (1 + np.exp(-x / T)) for x in a]
        Wq = [f + s * o for f, s, o in zip(floors, st, sig)]
        # forward with caches
        S = [V]
        Ys = []
        for l in range(L - 1):
            Yl = Wq[l] @ S[-1]
            Ys.append(Yl)
            S.append(gs2(Yl))
        out = Wq[-1] @ S[-1]
        R = out - tgt
        loss = np.mean(np.sum(R**2, axis=0))
        hist.append(loss)
        g = 2 * R / bs
        grads = [None] * L
        grads[L - 1] = g @ S[-1].T
        gS = Wq[-1].T @ g
        for l in range(L - 2, -1, -1):
            Yl = Ys[l]
            sw = (Yl[0::2] >= Yl[1::2])  # max picks the even row
            gY = np.empty_like(gS)
            gY[0::2] = np.where(sw, gS[0::2], gS[1::2])
            gY[1::2] = np.where(sw, gS[1::2], gS[0::2])
            grads[l] = gY @ S[l].T
            gS = Wq[l].T @ gY
        for l in range(L):
            ga = grads[l] * st[l] * sig[l] * (1 - sig[l]) / T
            mA[l] = b1 * mA[l] + (1 - b1) * ga
            vA[l] = b2 * vA[l] + (1 - b2) * ga**2
            mh = mA[l] / (1 - b1 ** (it + 1))
            vh = vA[l] / (1 - b2 ** (it + 1))
            a[l] = a[l] - lr * mh / (np.sqrt(vh) + epsA)
    offsets = [(x > 0).astype(float) for x in a]
    return offsets, np.array(hist)


# --------------------------------------------------------------------- policy
class Policy:
    """Teacher T(z), Student Psi_B(z) and LUT Psi_A acting on reduced states.

    lut_space = "state": Q_{b_A} on z in Z_ODD (step D_ODD/2^{b_A}, Section sec:problem item 1).
    lut_space = "port":  Q_{b_A} on v = g_r^T P z (Remark rem:lut_variants(i)), step D_v/2^{b_A}.
    Mid-tread lattice with half-open bins, outermost bins closed; the zero bin is the deadzone
    (Psi_A(0) = 0).  deadzone=False is the ablation: mid-rise lattice without a zero level.
    """

    def __init__(self, teacher: Teacher, qnet: QuantizedNet, portmap: np.ndarray, u_max: float,
                 D_ODD: float, b_A: int, deadzone: bool = True, lut_space: str = "state", D_v=None):
        self.t, self.q = teacher, qnet
        self.portmap = portmap  # g_r^T P   (m x K)
        self.u_max = u_max
        self.D_ODD = D_ODD
        self.b_A = int(b_A)
        self.lut_space = lut_space
        self.D_lut = D_ODD if lut_space == "state" else D_v
        self.dq = self.D_lut / 2**self.b_A
        self.jmax = 2 ** (self.b_A - 1)
        self.deadzone = deadzone

    def teacher(self, Z):
        return project_U(self.t.forward(self.portmap @ Z), self.u_max)

    def psi_B(self, Z):
        return project_U(self.q.forward(self.portmap @ Z), self.u_max)

    def _q(self, X):
        if self.deadzone:
            j = np.clip(np.floor(X / self.dq + 1 / 2), -self.jmax, self.jmax)
            return j * self.dq
        j = np.clip(np.floor(X / self.dq), -self.jmax, self.jmax - 1)
        return (j + 1 / 2) * self.dq

    def quantize_A(self, Z):
        """Implements eq:lut_capacity / eq:eps_A_triangle rounding: Q_{b_A} applied to the LUT index variable (z or v)."""
        return self._q(Z if self.lut_space == "state" else self.portmap @ Z)

    def psi_A(self, Z):
        Xq = self.quantize_A(Z)
        if self.lut_space == "state":
            U = self.psi_B(Xq)
        else:
            U = project_U(self.q.forward(Xq), self.u_max)
        if self.deadzone:
            zero = np.all(Xq == 0, axis=0)
            U[:, zero] = 0
        return U
