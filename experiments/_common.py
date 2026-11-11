"""Shared experiment harness: config loading, assertion ledger, tables, timing."""
from __future__ import annotations

import copy
import csv
import hashlib
import json
import math
import os
import re
import resource
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
RES = os.path.join(ROOT, "results")
FIG = os.path.join(RES, "figures")
TAB = os.path.join(RES, "tables")
LOG = os.path.join(RES, "logs")
for d in (FIG, TAB, LOG):
    os.makedirs(d, exist_ok=True)


def load_config(name="default"):
    with open(os.path.join(ROOT, "config", f"{name}.json")) as f:
        return json.load(f)


def override(cfg, **sections):
    """override(cfg, timing={"T_s": ...}) -> deep copy with the given keys replaced."""
    c = copy.deepcopy(cfg)
    for sec, kv in sections.items():
        c[sec].update(kv)
    return c


class Ledger:
    """Collects every assertion of an experiment; `finish` writes the log and raises if any failed."""

    def __init__(self, name):
        self.name = name
        self.t0 = time.perf_counter()
        self.items = []
        self.data = {}
        print(f"=== {name} ===", flush=True)

    def check(self, label, cond, **info):
        ok = bool(cond)
        self.items.append({"assertion": label, "passed": ok, **{k: _js(v) for k, v in info.items()}})
        print(f"[{'PASS' if ok else 'FAIL'}] {label} " + " ".join(f"{k}={_fmt(v)}" for k, v in info.items()), flush=True)
        return ok

    def finish(self):
        wall = time.perf_counter() - self.t0
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
        out = {"experiment": self.name, "wall_clock_s": wall, "peak_rss_MB": peak,
               "assertions": self.items, "data": _js(self.data)}
        with open(os.path.join(LOG, f"{self.name}.json"), "w") as f:
            json.dump(out, f, indent=1)
        n_fail = sum(not i["passed"] for i in self.items)
        print(f"--- {self.name}: {len(self.items) - n_fail}/{len(self.items)} assertions passed, "
              f"wall {wall:.1f} s, peak RSS {peak:.0f} MB", flush=True)
        assert n_fail == 0, f"{self.name}: {n_fail} assertion(s) failed"


def _js(v):
    import numpy as np
    if isinstance(v, dict):
        return {str(k): _js(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_js(x) for x in v]
    if isinstance(v, np.ndarray):
        return _js(v.tolist())
    if isinstance(v, (np.bool_,)):
        return bool(v)
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating, float)):
        v = float(v)
        return v if math.isfinite(v) else str(v)
    return v


def _fmt(v):
    try:
        return f"{float(v):.6g}"
    except (TypeError, ValueError):
        return str(v)


def write_csv(name, header, rows):
    """Machine-readable CSV only (data source of the main-paper figures; no LaTeX table)."""
    with open(os.path.join(TAB, f"{name}.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow([_cell(x) for x in r])


def write_table(name, header, rows, caption="", label="", wide=None, align=None, note=""):
    """CSV (raw values, machine-readable) + LaTeX (booktabs) table typeset through pgcbf.texnames (V1/V2):
    quantity keys -> math names, booleans -> checkmark / times, nan -> N/A, equation labels -> \\eqref / Lemma~\\ref,
    floats in math mode.  No \\resizebox; \\footnotesize; single-column float for <= 4 columns unless wide=True."""
    with open(os.path.join(TAB, f"{name}.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow([_cell(x) for x in r])
    wide = len(header) > 4 if wide is None else wide
    env = "table*" if wide else "table"
    align = align or ("l" + "r" * (len(header) - 1))
    lines = [r"\begin{%s}[t]" % env, r"\centering", r"\caption{" + caption + "}", r"\label{" + (label or "tab:" + name) + "}",
             r"\footnotesize", r"\setlength{\tabcolsep}{%dpt}" % (3 if len(header) >= 8 else 4), r"\begin{tabular}{" + align + "}", r"\toprule",
             " & ".join(head_tex(h) for h in header) + r" \\", r"\midrule"]
    for r in rows:
        if r == "midrule":
            lines.append(r"\midrule")
            continue
        lines.append(" & ".join(cell_tex(x) for x in r) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    if note:
        lines += [r"\\[2pt]", r"\parbox{\linewidth}{\scriptsize " + note + "}"]
    lines += [r"\end{%s}" % env]
    with open(os.path.join(TAB, f"{name}.tex"), "w") as f:
        f.write("\n".join(lines) + "\n")


def head_tex(h):
    from pgcbf.texnames import NAMES
    h = str(h)
    if h in NAMES:
        return NAMES[h]
    return _tex(h)


_LABEL = re.compile(r"^(eq|lem|prop|thm|rem|assump|cor|sec|def|app):[A-Za-z0-9_]+$")


def cell_tex(x):
    """One table cell -> LaTeX (V1)."""
    import numpy as np
    from pgcbf import texnames as tn
    if x is None:
        return "--"
    if isinstance(x, (bool, np.bool_)):
        return tn.boolean(x)
    if isinstance(x, (int, np.integer)):
        return "$%d$" % int(x)
    if isinstance(x, (float, np.floating)):
        return tn.num(x)
    s = str(x)
    if _LABEL.match(s):
        return tn.ref(s)
    if s in tn.NAMES and (s.startswith("H") or s in ("timing_ok", "codesign_window", "theorem_valid", "prop3_valid")):
        return tn.NAMES[s]
    return _tex(s)


def _cell(x):
    import numpy as np
    if x is None:
        return ""
    if isinstance(x, (bool, np.bool_)):
        return "True" if x else "False"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    if isinstance(x, (float, np.floating)):
        return repr(float(x))
    return str(x)


def _tex(s):
    """Text cell: math/LaTeX passes through, everything else is escaped (an underscore can never leak)."""
    s = str(s)
    if "$" in s or "\\" in s:
        return s
    if "_" in s:
        raise ValueError(f"raw identifier in table text (use a texnames key or LaTeX): {s!r}")
    esc = {"%": r"\%", "&": r"\&", "#": r"\#", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}",
           "|": r"\textbar{}", "<": r"\textless{}", ">": r"\textgreater{}", "{": r"\{", "}": r"\}", "$": r"\$"}
    return "".join(esc.get(ch, ch) for ch in s)


def array_hash(*arrs):
    import numpy as np
    h = hashlib.sha256()
    for a in arrs:
        h.update(np.ascontiguousarray(a).tobytes())
    return h.hexdigest()
