"""Build a LaTeX document from a {{key}} template (story bible Section 9; generalises tools/build_section.py).

    python tools/build_main.py --template paper/main_template.tex --out paper/main.tex

Every {{key}} placeholder is resolved from results/section_numbers.json (the "tex" rendering written by
tools/build_section.py from the logs of the last run); nothing is recomputed here.  The build fails (exit 1) on
  * any unresolved placeholder,
  * any unused new key: a key whose name starts with one of --new-prefixes (default: the Phase 4.1 keys e4b_, util_
    and e6_c_delta) that occurs neither in --template nor in any --also template (default: the Section V template
    results/experimental_section_template.tex, which feeds the supplement),
  * any texlint finding (tools/texlint.py lint_files) or any prose numeral of the output that does not occur among
    the values of section_numbers.json (texlint.unsourced_numerals).
"""
import argparse
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tools"))
from texlint import lint_files, unsourced_numerals  # noqa: E402

PLACEHOLDER = re.compile(r"\{\{([A-Za-z0-9_]+)\}\}")
NEW_PREFIXES = "e4b_,util_,e6_c_delta,budget_,cfg_,tauA_R1,tauB_R1"
SECTION_TEMPLATE = os.path.join(ROOT, "results", "experimental_section_template.tex")


def load_numbers(path):
    with open(path) as fh:
        d = json.load(fh)
    return {k: v["tex"] for k, v in d.items()}


def build(template_src, numbers, also_srcs=(), new_prefixes=()):
    """Return (output, issues): output with every resolvable placeholder replaced, issues the list of failures."""
    keys = set(PLACEHOLDER.findall(template_src))
    issues = [f"unresolved placeholder {{{{{k}}}}}" for k in sorted(keys - set(numbers))]
    used = set(keys)
    for s in also_srcs:
        used |= set(PLACEHOLDER.findall(s))
    for k in sorted(numbers):
        if any(k.startswith(p) for p in new_prefixes) and k not in used:
            issues.append(f"unused new key {k!r}")
    out = PLACEHOLDER.sub(lambda m: numbers.get(m.group(1), m.group(0)), template_src)
    issues += [f"numeral {t!r} in the prose does not occur in section_numbers.json"
               for t in unsourced_numerals(out, list(numbers.values()))]
    return out, issues


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--template", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--numbers", default=os.path.join(ROOT, "results", "section_numbers.json"))
    ap.add_argument("--also", action="append", default=None,
                    help="further templates whose placeholders count as uses of new keys (default: Section V template)")
    ap.add_argument("--new-prefixes", default=NEW_PREFIXES, help="comma-separated; empty string disables the check")
    a = ap.parse_args(argv)
    numbers = load_numbers(a.numbers)
    with open(a.template) as fh:
        tpl = fh.read()
    also = a.also if a.also is not None else ([SECTION_TEMPLATE] if os.path.exists(SECTION_TEMPLATE) else [])
    also_srcs = []
    for p in also:
        with open(p) as fh:
            also_srcs.append(fh.read())
    prefixes = tuple(p for p in a.new_prefixes.split(",") if p)
    out, issues = build(tpl, numbers, also_srcs, prefixes)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write(out)
    issues += lint_files([a.out])
    for x in issues:
        print(x)
    n = len(set(PLACEHOLDER.findall(tpl)))
    print(f"build_main: {a.template} -> {a.out}: {n} placeholders, {len(issues)} issue(s)")
    return 1 if issues else 0


if __name__ == "__main__":
    sys.exit(main())
