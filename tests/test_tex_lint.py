"""V1: generated LaTeX must not leak raw identifiers, booleans, nan or plain-text labels into the PDF."""
import glob
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tools"))
from texlint import lint_text, lint_files  # noqa: E402

TABLES = sorted(glob.glob(os.path.join(ROOT, "results", "tables", "*.tex")))
SECTION = os.path.join(ROOT, "results", "experimental_section.tex")


def test_linter_is_not_vacuous():
    bad = ["eps_basis_P = 0.1", "valid: True", "value nan here", "see eq:teacher_cbf", "None", "hypothesis (H8)",
           r"\textit{K_m}"]
    for s in bad:
        assert lint_text(s), s
    good = [r"$\varepsilon_{\mathrm{basis},P}=0.1$ \checkmark", r"see \eqref{eq:teacher_cbf} and Lemma~\ref{lem:delta_cert_classH}",
            r"\label{tab:exp1_constants}", r"N/A", r"$K_m$ and \(K_f\)", r"\includegraphics[width=\linewidth]{figures/exp4_theorem1.pdf}"]
    for s in good:
        assert not lint_text(s), (s, lint_text(s))


def test_generated_tables_exist():
    assert len(TABLES) >= 15, TABLES


@pytest.mark.parametrize("path", TABLES, ids=[os.path.basename(p) for p in TABLES])
def test_generated_table_is_clean(path):
    issues = lint_files([path])
    assert not issues, "\n".join(issues)


def test_tables_have_no_resizebox_and_footnotesize_minimum():
    for p in TABLES:
        s = open(p).read()
        assert "\\resizebox" not in s, p
        assert "\\scriptsize" not in s.replace("\\parbox{\\linewidth}{\\scriptsize", ""), p
        assert "\\tiny" not in s, p


def test_tables_have_at_most_seven_columns_unless_wide():
    import re
    for p in TABLES:
        s = open(p).read()
        m = re.search(r"\\begin\{tabular\}\{([^}]*)\}", s)
        ncol = len(re.sub(r"[^lcrp]", "", m.group(1)))
        assert ncol <= 7 or "table*" in s, (p, ncol)
        assert ncol <= 10, (p, ncol)


@pytest.mark.skipif(not os.path.exists(SECTION), reason="Section V not built yet (tools/build_section.py)")
def test_experimental_section_is_clean():
    issues = lint_files([SECTION])
    assert not issues, "\n".join(issues)


def test_numeral_provenance_checker_is_not_vacuous():
    from texlint import unsourced_numerals
    vals = ["0.3725", "2", r"9.994\times10^{-4}", "1"]
    good = r"(H0)--(H7) hold in R2; ETD2RK; \(c_M = 0.3725\), margin \(9.994\times10^{-4}\), \(P_0^{-1}\), \(m = 2\)"
    assert unsourced_numerals(good, vals) == [], unsourced_numerals(good, vals)
    for bad in (r"\(c_M = 0.373\)", r"a margin of 3.1\%", r"\(9.99\times10^{-4}\)", r"about 46.5 percent"):
        assert unsourced_numerals(bad, vals), bad


@pytest.mark.skipif(not os.path.exists(SECTION), reason="Section V not built yet (tools/build_section.py)")
def test_every_prose_numeral_occurs_in_section_numbers():
    """Review item 8: every numeral of the generated Section V prose occurs in results/section_numbers.json."""
    import json
    from texlint import unsourced_numerals
    nums = json.load(open(os.path.join(ROOT, "results", "section_numbers.json")))
    bad = unsourced_numerals(open(SECTION).read(), [v["tex"] for v in nums.values()])
    assert not bad, bad
