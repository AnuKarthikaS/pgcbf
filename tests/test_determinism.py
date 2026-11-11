"""Bit-for-bit reproducibility: the same closed-loop (worst-case adversary) and stochastic runs executed
in two independent interpreter processes produce identical SHA-256 hashes; plus exact latch timing."""
import os
import subprocess
import sys

import numpy as np
from scipy.linalg import expm

from conftest import ROOT, load
from pgcbf.certificates import build

SCRIPT = r"""
import sys, json, hashlib, numpy as np
sys.path.insert(0, %r)
from pgcbf.certificates import build, certificate
from pgcbf.simulate import Simulator, initial_conditions
from pgcbf.stochastic import noise_matrix
cfg = json.load(open(%r))
s = build(cfg, train=True); c = certificate(cfg, s); sc = cfg["sim"]
X0 = initial_conditions(s["lat"], s["rom"], c, int(sc["n_seeds"]), np.random.default_rng(int(cfg["seeds"]["ic"])),
                        float(sc["ic_level_fraction"]), float(sc["ic_tail_fraction"]), int(sc["ic_n_tail_modes"]))
n = int(cfg["experiments"]["latch_periods"])
r = Simulator(s["lat"], s["rom"], s["policy"], cfg, adversary="worst").run(X0, n, record="substep")
Sig = noise_matrix(s["lat"], s["rom"], cfg["stochastic"]["noise_small"])
rs = Simulator(s["lat"], s["rom"], s["policy"], cfg, adversary=None, noise=Sig).run(X0, n, np.random.default_rng(int(cfg["seeds"]["mc"])), record="period")
h = hashlib.sha256()
for a in (r["h"], r["X_final"], rs["h"], rs["X_final"], *s["qnet"].Wq):
    h.update(np.ascontiguousarray(a).tobytes())
print(h.hexdigest())
"""


def test_two_processes_identical_sha256():
    code = SCRIPT % (os.path.join(ROOT, "src"), os.path.join(ROOT, "config", "default.json"))
    env = dict(os.environ)
    outs = [subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env, check=True).stdout.strip()
            for _ in range(2)]
    assert len(outs[0]) == 64 and outs[0] == outs[1]


def test_exact_latch_timing():
    from pgcbf.simulate import Simulator
    cfg = load("default")
    cfg["plant"]["kappa"] = 0
    s = build(cfg, train=True)
    lat = s["lat"]
    rng = np.random.default_rng(int(cfg["seeds"]["mc"]))
    uA, uB = rng.uniform(-1, 1, (lat.m, 1)), rng.uniform(-1, 1, (lat.m, 1))

    class Const:
        D_ODD = np.inf

        def psi_A(self, Z):
            return np.repeat(uA, Z.shape[1], axis=1)

        def psi_B(self, Z):
            return np.repeat(uB, Z.shape[1], axis=1)

    x0 = s["rom"].Pperp @ rng.standard_normal(lat.N) + s["rom"].Phi @ rng.standard_normal(s["rom"].K)
    npl = int(cfg["experiments"]["latch_periods"])
    xs = Simulator(lat, s["rom"], Const(), cfg, adversary=None).run(x0[:, None], npl, record="period")["X_final"][:, 0]
    Ts, tA, tB = (float(cfg["timing"][k]) for k in ("T_s", "tau_A", "tau_B"))

    def step(x, u, dt):
        M = np.zeros((lat.N + 1, lat.N + 1))
        M[:lat.N, :lat.N] = lat.A_lin * dt
        M[:lat.N, lat.N] = (lat.G @ u)[:, 0] * dt
        E = expm(M)
        return E[:lat.N, :lat.N] @ x + E[:lat.N, lat.N]

    x = x0.copy()
    for _ in range(npl):
        x = step(step(step(x, uB, tA), uA, tB - tA), uB, Ts - tB)
    assert np.linalg.norm(xs - x) <= float(cfg["tolerances"]["roundoff_rel"]) * npl * lat.N * np.linalg.norm(x)
