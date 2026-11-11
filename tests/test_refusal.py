"""The certificate must refuse (raise) on every non-admissible configuration, and accept the deployed ones."""
import copy

import pytest

from conftest import load
from pgcbf.certificates import CertificateRefused, build, certificate, refuse_if_invalid


def _cert(cfg, **sections):
    c = copy.deepcopy(cfg)
    for sec, kv in sections.items():
        c[sec].update(kv)
    return certificate(c, build(c, train=True))


@pytest.mark.parametrize("name", ["default", "r1_passive"])
def test_deployed_configs_accepted(name):
    assert refuse_if_invalid(_cert(load(name)))


def test_refusals(cfg):
    ex = cfg["experiments"]
    base = certificate(cfg, build(cfg, train=True))
    cases = [
        dict(controller={"lut_space": "state"}),                                         # codesign window empty
        dict(controller={"u_max": base["H_sat_lhs"] * ex["u_max_refusal_fraction"]}),    # (H-sat) fails
        dict(timing={k: v * ex["Ts_refusal_factor"] for k, v in cfg["timing"].items()}),  # hold term, (H7)
        dict(controller={"b_B": int(ex["ablation_b_B"])}),                               # quantisation, (H7)
        dict(rom={"K_m": int(ex["ablation_K_m"])}),                                      # tail / rank
        dict(safety={"eps_adv": ex["ablation_eps_adv"]}),                                # adversary
    ]
    for kw in cases:
        with pytest.raises(CertificateRefused):
            refuse_if_invalid(_cert(cfg, **kw))


def test_timing_order_refused(cfg):
    t = cfg["timing"]
    c = _cert(cfg, timing={"tau_A": t["tau_B"], "tau_B": t["tau_A"]})
    assert not c["timing_ok"]
    with pytest.raises(CertificateRefused):
        refuse_if_invalid(c)
