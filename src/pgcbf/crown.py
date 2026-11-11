"""Twin-CROWN bound for the Student-Teacher gap (Proposition prop:twin_crown).

Implements eq:twin_affine_relax, eq:twin_slopes, eq:twin_crown_bound and
eq:crown_explicit.  The twin network propagates (y_l, delta y_l) jointly:

    y_l      = W_l s_{l-1}                  s_l  = GS2(y_l),  s_0 = v = g_r^T P z
    delta y_l = W_l^q delta s_{l-1} + DeltaW_l s_{l-1} + e_l,  delta s_0 = e_0

and GS2 on a pair (a, b) is written with r = ReLU(a - b):
    s = [b + r ; a - r],  delta s = [delta b + nu ; delta a - nu],
    nu = ReLU(d^q) - ReLU(d),  d = a - b,  d^q = d + delta d.
The teacher ReLU r uses the standard CROWN relaxation; nu uses eq:twin_slopes.
"""
from __future__ import annotations

import itertools

import math

import numpy as np


def _relu_relax(l, u):
    """Standard single-network CROWN secant relaxation of ReLU(d), d in [l, u]."""
    aU = np.zeros_like(l); bU = np.zeros_like(l); aL = np.zeros_like(l)
    pos = l >= 0
    mix = (l < 0) & (u > 0)
    aU[pos] = 1; aL[pos] = 1
    den = np.where(mix, u - l, 1)
    aU[mix] = (u / den)[mix]
    bU[mix] = (-u * l / den)[mix]
    aL[mix] = (u[mix] >= -l[mix]).astype(float)
    return aL, np.zeros_like(l), aU, bU


def twin_slopes(dl, du, l, u, ql, qu):
    """Implements eq:twin_slopes: (alpha^L, beta^L, alpha^U, beta^U) for nu = ReLU(d+dd) - ReLU(d)."""
    aL = np.zeros_like(dl); bL = np.zeros_like(dl); aU = np.zeros_like(dl); bU = np.zeros_like(dl)
    Iplus = (l >= 0) & (ql >= 0)
    Iminus = (u <= 0) & (qu <= 0)
    star = ~(Iplus | Iminus)
    c1 = star & (dl >= 0)
    c2 = star & ~c1 & (du <= 0)
    c3 = star & ~c1 & ~c2
    aL[Iplus] = 1; aU[Iplus] = 1
    aU[c1] = 1
    aL[c2] = 1
    den = np.where(c3, du - dl, 1)
    aL[c3] = (-dl / den)[c3]
    bL[c3] = (dl * du / den)[c3]
    aU[c3] = (du / den)[c3]
    bU[c3] = (-dl * du / den)[c3]
    return aL, bL, aU, bU


class TwinCrown:
    def __init__(self, qnet, portmap):
        self.qn = qnet
        self.t = qnet.t
        self.Mz = portmap  # v = Mz z
        self.L = self.t.L

    def _backward(self, Lam, Mu, k_top, bnds, zc, zr):
        """Upper bound of Lam y_k + Mu delta y_k over the box (centre zc, radius zr)."""
        const = np.zeros(Lam.shape[0])
        W, Wq, dW, rho = self.t.W, self.qn.Wq, self.qn.dW, self.qn.rho
        for k in range(k_top, 0, -1):
            const = const + np.abs(Mu).sum(axis=1) * rho[k]
            As = Lam @ W[k - 1] + Mu @ dW[k - 1]
            Ad = Mu @ Wq[k - 1]
            if k == 1:
                C = As @ self.Mz
                return const + C @ zc + np.abs(C) @ zr + np.abs(Ad).sum(axis=1) * rho[0]
            b = bnds[k - 1]  # bounds at hidden layer k-1
            a = As[:, 0::2] - As[:, 1::2]
            c = Ad[:, 0::2] - Ad[:, 1::2]
            Lam_n = np.empty_like(As); Mu_n = np.empty_like(Ad)
            Lam_n[:, 0::2] = As[:, 1::2]; Lam_n[:, 1::2] = As[:, 0::2]
            Mu_n[:, 0::2] = Ad[:, 1::2]; Mu_n[:, 1::2] = Ad[:, 0::2]
            rl, rbl, ru, rbu = b["relu"]
            sl = np.where(a >= 0, ru, rl); sb = np.where(a >= 0, rbu, rbl)
            Lam_n[:, 0::2] += a * sl; Lam_n[:, 1::2] -= a * sl
            const = const + (a * sb).sum(axis=1)
            tl, tbl, tu, tbu = b["twin"]
            sl = np.where(c >= 0, tu, tl); sb = np.where(c >= 0, tbu, tbl)
            Mu_n[:, 0::2] += c * sl; Mu_n[:, 1::2] -= c * sl
            const = const + (c * sb).sum(axis=1)
            Lam, Mu = Lam_n, Mu_n
        raise RuntimeError("unreachable")

    def bounds_box(self, zlo, zhi):
        """Return (D_low, D_up) for the output gap over one box of the partition."""
        zc, zr = (zlo + zhi) / 2, (zhi - zlo) / 2
        H = self.t.H
        h2 = H // 2
        S = np.zeros((h2, H))
        S[np.arange(h2), 2 * np.arange(h2)] = 1
        S[np.arange(h2), 2 * np.arange(h2) + 1] = -1
        Z = np.zeros_like(S)
        bnds = {}
        for k in range(1, self.L):
            Lam = np.vstack([S, -S, Z, Z, S, -S])
            Mu = np.vstack([Z, Z, S, -S, S, -S])
            ub = self._backward(Lam, Mu, k, bnds, zc, zr)
            u, lneg, du, dlneg, qu, qlneg = np.split(ub, Lam.shape[0] // h2)  # the six stacked bound blocks
            l, dl, ql = -lneg, -dlneg, -qlneg
            qu = np.minimum(qu, u + du); ql = np.maximum(ql, l + dl)
            bnds[k] = {"d": (l, u), "dd": (dl, du), "dq": (ql, qu),
                       "relu": _relu_relax(l, u), "twin": twin_slopes(dl, du, l, u, ql, qu)}
        m = self.t.m
        I = np.eye(m)
        ub = self._backward(np.zeros((2 * m, m)), np.vstack([I, -I]), self.L, bnds, zc, zr)
        return -ub[m:], ub[:m], bnds

    def eps_crown(self, D_ODD: float, K: int, splits: int):
        """Implements eq:twin_crown_bound over a uniform partition of Z_ODD into splits^K boxes."""
        edges = np.linspace(-D_ODD / 2, D_ODD / 2, splits + 1)
        best = 0
        boxes = []
        for idx in itertools.product(range(splits), repeat=K):
            idx = np.array(idx)
            lo, hi = edges[idx], edges[idx + 1]
            Dl, Du, _ = self.bounds_box(lo, hi)
            val = np.sqrt(np.sum(np.maximum(Dl**2, Du**2)))
            boxes.append((lo, hi, Dl, Du))
            best = max(best, val)
        return best, boxes


def E_q(qnet, Y: float) -> float:
    """Implements eq:crown_explicit  E_q(b_B) = sum_l (Delta_W,l Y + Delta_b,l + sqrt(H_l) rho_l) prod_{k>l} ||W_k^q||."""
    L = qnet.t.L
    dW = [0] + [np.linalg.norm(d, 2) for d in qnet.dW]
    db = [0] + [np.linalg.norm(b) for b in qnet.bq]
    norms = [np.linalg.norm(Wq, 2) for Wq in qnet.Wq]  # index k-1 for layer k
    tot = 0
    for l in range(L + 1):
        prod = 1
        for k in range(l + 1, L + 1):
            prod *= norms[k - 1]
        tot += (dW[l] * Y + db[l] + np.sqrt(qnet.widths[l]) * qnet.rho[l]) * prod
    return tot
