"""Experiment 0 - the benchmark is genuinely unsafe open loop (u = 0).

R2 benchmark (Pe = c ell / nu set in config/default.json).  From certified initial states
(H7) the uncontrolled full-order lattice leaves C = {h >= 0} within the horizon.  The worst-case
adversary is switched off here, so the violation is due to the plant alone."""
import numpy as np
from _common import Ledger, load_config, FIG, write_table
from pgcbf.certificates import build, certificate
from pgcbf.simulate import Simulator, initial_conditions
from pgcbf import plotting as pl

L = Ledger("exp0_openloop")
cfg = load_config("default")
sysd = build(cfg, train=True)
cert = certificate(cfg, sysd)
lat, rom, pol = sysd["lat"], sysd["rom"], sysd["policy"]
sc = cfg["sim"]
growth = float(np.linalg.eigvals(lat.A_lin).real.max())
L.check("open-loop linearisation is exponentially unstable (max Re eig > 0)", growth > 0, max_re_eig=growth, Pe=lat.peclet)
X0 = initial_conditions(lat, rom, cert, int(sc["n_seeds"]), np.random.default_rng(int(cfg["seeds"]["ic"])),
                        float(sc["ic_level_fraction"]), float(sc["ic_tail_fraction"]), int(sc["ic_n_tail_modes"]))
h0 = lat.h_barrier(X0, cert["E_safe"])
L.check("initial states lie in C (h(x0) >= 0)", np.all(h0 >= 0), min_h0=h0.min())
sim = Simulator(lat, rom, pol, cfg, adversary=None, controller=False)
n = int(round(float(sc["horizon"]) / sim.Ts))
r = sim.run(X0, n, record="substep")
mh = r["h"].min(axis=0)
first = np.array([r["t"][np.argmax(r["h"][:, j] < 0)] if mh[j] < 0 else np.nan for j in range(X0.shape[1])])
L.check("open loop exits C for every seed within the horizon", np.all(mh < 0), n_exit=int(np.sum(mh < 0)),
        n_seeds=X0.shape[1], worst_min_h=mh.min(), latest_exit_time=np.nanmax(first))
L.data.update(min_h=mh, first_exit=first, growth=growth, Pe=lat.peclet, horizon=n * sim.Ts)
write_table("exp0_openloop", ["seed", "$h(x_0)$", "$\\min_t h(x(t))$", "first exit time"],
            [[j, h0[j], mh[j], first[j]] for j in range(X0.shape[1])], wide=False,
            caption=f"Experiment 0: open-loop ($u\\equiv0$) trajectories of the R2 benchmark (Pe $= {lat.peclet:.3g}$) from certified initial states.")
fig, ax = pl.figure(1, 1, 3.5, 2.4)
for j in range(X0.shape[1]):
    ax[0, 0].plot(r["t"], r["h"][:, j])
pl.boundary(ax[0, 0])
ax[0, 0].set_xlabel(r"$t$"); ax[0, 0].set_ylabel(r"$h(x(t))$"); ax[0, 0].set_title(f"Open loop, $\\mathrm{{Pe}}={lat.peclet:.3g}$: all seeds (raw)")
ax[0, 0].legend(loc="lower left")
pl.save(fig, FIG, "exp0_openloop_h")
L.finish()
