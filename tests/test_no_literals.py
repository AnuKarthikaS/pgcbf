"""Hard constraint 1: no certificate quantity may be a literal.

(i)  AST scan: the numerical modules of src/pgcbf (everything except plotting.py) contain no float
     or complex literal and no float("...") call.  Integer literals are allowed only as exact
     mathematical coefficients of the implemented formulas (2 sqrt(E), /6, ...).
(ii) AST scan of experiments/: float literals may appear only as figure geometry or line-style
     keywords (cosmetics), never in a computation or an assertion.
(iii) Perturbation test: every float-valued certificate quantity changes when at least one
     problem parameter of the config file is perturbed, i.e. none of them is a hidden constant.
"""
import ast
import copy
import glob
import math
import os

import numpy as np

from conftest import ROOT
from pgcbf.certificates import build, certificate

NUMERICAL = [f for f in glob.glob(os.path.join(ROOT, "src", "pgcbf", "*.py")) if not f.endswith("plotting.py")]
COSMETIC_KW = {"alpha", "lw", "linewidth", "ms", "markersize", "fontsize", "s", "bbox_to_anchor", "linthresh"}


def _float_literals(tree, allowed=frozenset()):
    bad = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, (float, complex)) and id(node) not in allowed:
            bad.append((node.lineno, node.value))
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in ("float", "complex")
                and node.args and isinstance(node.args[0], ast.Constant) and isinstance(node.args[0].value, str)):
            bad.append((node.lineno, f'{node.func.id}("{node.args[0].value}")'))
    return bad


def test_no_float_literals_in_numerical_modules():
    assert len(NUMERICAL) >= 7
    offenders = {}
    for f in NUMERICAL:
        bad = _float_literals(ast.parse(open(f).read()))
        if bad:
            offenders[os.path.basename(f)] = bad
    assert not offenders, f"float literals in certificate code: {offenders}"


def test_literal_scanner_detects_literals():
    """The scanner itself is not vacuous."""
    assert sorted(_float_literals(ast.parse("x = 0.5 * y\nz = float('inf')")), key=str) == [(1, 0.5), (2, 'float("inf")')]


def _cosmetic_nodes(tree):
    ok = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            fn = node.func
            if isinstance(fn, ast.Attribute) and fn.attr == "figure":
                ok.update(id(a) for a in node.args)
            for kw in node.keywords:
                if kw.arg in COSMETIC_KW:
                    ok.update(id(n) for n in ast.walk(kw.value))
        if isinstance(node, ast.FunctionDef) and node.name in ("_fmt", "_cell"):  # number formatting of logs / tables
            ok.update(id(n) for n in ast.walk(node))
    return ok


def test_experiment_scripts_use_floats_only_for_cosmetics():
    offenders = {}
    for f in glob.glob(os.path.join(ROOT, "experiments", "*.py")):
        tree = ast.parse(open(f).read())
        bad = _float_literals(tree, frozenset(_cosmetic_nodes(tree)))
        if bad:
            offenders[os.path.basename(f)] = bad
    assert not offenders, offenders


PERTURB = [("plant", "nu"), ("plant", "gamma"), ("plant", "c"), ("plant", "kappa"), ("plant", "actuator_width"),
           ("plant", "disturbance_width"), ("safety", "E_safe"), ("safety", "gamma0"), ("safety", "eps_adv"),
           ("safety", "odd_margin"), ("safety", "h_act_fraction"), ("controller", "u_max"), ("controller", "L_T"),
           ("controller", "theta"), ("timing", "T_s"), ("timing", "tau_A"), ("timing", "tau_B"), ("plant", "adv_cutoff_weight")]
INT_PERTURB = [("controller", "b_A"), ("controller", "b_B"), ("plant", "N_m"), ("rom", "K_m")]
# quantities that are round-off residuals of exact identities (they must be ~0, not constant)
RESIDUALS = {"Q_identity_relerr", "H2a_residual", "H2a_invariance_residual"}


def test_every_certificate_quantity_depends_on_the_inputs(cfg):
    base = certificate(cfg, build(cfg, train=False))
    keys = [k for k, v in base.items() if isinstance(v, float) and math.isfinite(v) and k not in RESIDUALS]
    assert len(keys) >= 50
    rel = float(cfg["tolerances"]["roundoff_rel"])
    changed = set()
    for sec, key in PERTURB + INT_PERTURB:
        c = copy.deepcopy(cfg)
        v = c[sec][key]
        c[sec][key] = v + 1 if isinstance(v, int) and (sec, key) in INT_PERTURB else v * (1 + 1 / 64)
        if (sec, key) == ("timing", "tau_A"):
            c[sec][key] = v * (1 - 1 / 64)
        cc = certificate(c, build(c, train=False))
        for k in keys:
            x, y = base[k], cc.get(k)
            if isinstance(y, float) and math.isfinite(y) and abs(x - y) > rel * max(abs(x), abs(y), np.finfo(float).tiny):
                changed.add(k)
    frozen = sorted(set(keys) - changed)
    assert not frozen, f"certificate quantities insensitive to every input (hidden literals?): {frozen}"
