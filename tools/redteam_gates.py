"""Red-team gates for the pgcbf repository.
Usage: python redteam_gates.py <repo_root>
Gates:
  G1  Hardcoded-constant scan: numeric literals (other than 0,1,2,-1,0.5 and small ints used as
      indices/exponents) inside src/**/certificates|crown|rom|controller|stochastic modules.
  G2  Equation-label traceability: every eq:/lem:/prop:/thm: label in the theory LaTeX that defines a
      certificate quantity must appear in a docstring in src/.
  G3  Clean-room: no reference to quarantined filenames anywhere in the repo.
  G4  Assertions present: count of assert / pytest checks per experiment file >= 1.
  G5  Plot hygiene: no set_ylim/set_xlim with literal bounds that could clip, no savgol/smooth calls.
  G6  Determinism: seeds set from config, no unseeded np.random.* default_rng() calls without seed.
"""
import ast, os, re, sys, pathlib

root = pathlib.Path(sys.argv[1])
src = root / "src"
ALLOWED = {0, 1, 2, -1, 0.5, 1.0, 2.0, 0.0, 1e-12, 1e-14, 3, 4, 10, 100}
CERT_MODULES = ("certificates", "crown", "rom", "controller", "stochastic", "lattice", "simulate")

def gate1():
    hits = []
    for p in src.rglob("*.py"):
        if not any(m in p.name for m in CERT_MODULES):
            continue
        tree = ast.parse(p.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)) and not isinstance(node.value, bool):
                if node.value not in ALLOWED:
                    hits.append((p.relative_to(root), node.lineno, node.value))
    return hits

def gate2(labels):
    text = "\n".join(p.read_text() for p in src.rglob("*.py"))
    return [l for l in labels if l not in text]

def gate3():
    bad = []
    for p in root.rglob("*"):
        if p.is_file() and p.resolve() != pathlib.Path(__file__).resolve() and p.suffix in {".py", ".md", ".tex", ".sh", ".json", ".yaml", ".yml", ".txt"}:
            t = p.read_text(errors="ignore")
            if "sample_experiment" in t:
                bad.append(p.relative_to(root))
    return bad

def gate4():
    out = {}
    for p in (root / "experiments").rglob("*.py"):
        t = p.read_text()
        out[p.name] = len(re.findall(r"\bassert\b|pytest\.|raise\s+AssertionError|check\(", t))
    return out

def _code_lines(p):
    """Source lines with docstrings and comments removed (the gate inspects code, not prose)."""
    import io, tokenize
    keep = {}
    toks = list(tokenize.generate_tokens(io.StringIO(p.read_text()).readline))
    for j, t in enumerate(toks):
        if t.type == tokenize.COMMENT:
            continue
        if t.type == tokenize.STRING and (j == 0 or toks[j - 1].type in (tokenize.INDENT, tokenize.NEWLINE, tokenize.NL, tokenize.DEDENT)):
            continue  # docstring / bare string statement
        keep.setdefault(t.start[0], []).append(t.string)
    return sorted((i, " ".join(v)) for i, v in keep.items())


def _plots(p):
    t = p.read_text()
    return "matplotlib" in t or "plotting" in t or "pl." in t


def gate5():
    hits = []
    for p in list(src.rglob("*.py")) + list((root / "experiments").rglob("*.py")):
        if not _plots(p):
            continue  # plot hygiene applies to plotting code (np.clip in the LUT quantiser is not a plot clip)
        for i, line in _code_lines(p):
            if re.search(r"set_[xy]lim\(\s*-?\d|ylim\(\s*-?\d|savgol|smooth|rolling\(|clip\(", line):
                hits.append((p.relative_to(root), i, line.strip()[:100]))
    return hits

def gate6():
    hits = []
    for p in list(src.rglob("*.py")) + list((root / "experiments").rglob("*.py")):
        for i, line in enumerate(p.read_text().splitlines(), 1):
            if re.search(r"default_rng\(\s*\)|np\.random\.(rand|randn|normal|uniform|choice|seed)\(", line):
                hits.append((p.relative_to(root), i, line.strip()[:100]))
    return hits

if __name__ == "__main__":
    labels = ["eq:delta_total_exact", "eq:delta_dip", "eq:exp_weights", "eq:eps_A_triangle", "eq:twin_crown_bound",
              "eq:eps_basis_def", "eq:codesign_window", "eq:mu_p0", "eq:rss_ball", "eq:delta_stoch", "eq:teacher_cbf",
              "eq:control_authority", "eq:tail_dissipation", "eq:reduced_dynamics", "eq:multirate_hold", "lem:delta_cert_classH"]
    print("G1 hardcoded literals:", *gate1(), sep="\n  ")
    print("G2 missing labels in src docstrings:", gate2(labels))
    print("G3 clean-room violations:", gate3())
    print("G4 assertions per experiment:", gate4())
    print("G5 plot hygiene hits:", *gate5(), sep="\n  ")
    print("G6 unseeded RNG hits:", *gate6(), sep="\n  ")
    g4 = gate4()
    n_hits = len(gate1()) + len(gate2(labels)) + len(gate3()) + sum(v < 1 for k, v in g4.items() if k.startswith("exp")) \
        + len(gate5()) + len(gate6())
    print("red-team gate hits:", n_hits)
    sys.exit(1 if n_hits else 0)
