"""V1 linter for generated LaTeX (results/tables/*.tex, results/experimental_section.tex).

Fails on any raw identifier that would be visible in the PDF:
  * an underscore outside math mode (escaped \\_ is also refused: identifiers must be typeset as math),
  * the words True, False, nan, None (booleans must be \\checkmark / $\\times$, missing values N/A or --),
  * a plain-text cross-reference label (eq:..., lem:..., thm:..., ...) outside \\ref/\\eqref/\\label,
  * the retired hypothesis name (H8) and the retired keys theorem_abc_valid / theorem_full_valid.
Usage: python tools/texlint.py file1.tex [file2.tex ...]  (exit 1 on any finding)."""
import re
import sys

_ARG_CMDS = r"label|ref|eqref|cref|Cref|cite|includegraphics|pgcfig|input|include|begin|end|usepackage|documentclass|newlabel|url|texttt"
_LABEL_WORD = re.compile(r"\b(eq|lem|prop|thm|rem|assump|cor|sec|def|app|tab|fig):[A-Za-z0-9_]+")
_BAD_WORD = re.compile(r"\b(True|False|nan|NaN|None)\b")


def strip_math(s: str) -> str:
    """Remove comments, math segments and arguments of reference-like commands; keep visible text."""
    s = re.sub(r"(?<!\\)%.*", "", s)
    for env in ("equation", "equation*", "align", "align*", "gather", "gather*", "multline", "eqnarray"):
        s = re.sub(r"\\begin\{%s\}.*?\\end\{%s\}" % (re.escape(env), re.escape(env)), " ", s, flags=re.S)
    s = re.sub(r"\\\[.*?\\\]", " ", s, flags=re.S)
    s = re.sub(r"\\\(.*?\\\)", " ", s, flags=re.S)
    s = re.sub(r"(?<!\\)\$\$.*?(?<!\\)\$\$", " ", s, flags=re.S)
    s = re.sub(r"(?<!\\)\$.*?(?<!\\)\$", " ", s, flags=re.S)
    s = re.sub(r"\\(%s)\*?(\[[^\]]*\])?\{[^}]*\}" % _ARG_CMDS, " ", s)
    s = re.sub(r"\\begin\{tabular\}\{[^}]*\}", " ", s)
    return s


def lint_text(src: str, name: str = "<string>"):
    out = []
    vis = strip_math(src)
    for i, line in enumerate(vis.splitlines(), 1):
        if "_" in line:
            out.append(f"{name}:{i}: underscore outside math: {line.strip()[:100]}")
        m = _BAD_WORD.search(line)
        if m:
            out.append(f"{name}:{i}: raw value {m.group(0)!r}: {line.strip()[:100]}")
        m = _LABEL_WORD.search(line)
        if m:
            out.append(f"{name}:{i}: plain-text label {m.group(0)!r}: {line.strip()[:100]}")
    for bad in ("(H8)", "theorem_abc_valid", "theorem_full_valid", r"\textrm{(H8)}", r"\mathrm{(H8)}"):
        if bad in src:
            out.append(f"{name}: retired name {bad!r}")
    return out


# ---------------------------------------------------------------- numeral provenance (review item 8)
_SCI = re.compile(r"\d+(?:\.\d+)?\\times10\^\{-?\d+\}")
_NUM = re.compile(r"(?<![A-Za-z0-9.])\d+(?:\.\d+)?(?![A-Za-z0-9])")
# visible tokens that contain digits but are names, not numbers
_NAME_TOKENS = re.compile(
    r"\((?:H\d[a-z]?|H-dz|H-sat)\)(?:--\((?:H\d[a-z]?)\))?"          # hypothesis labels (H0)--(H7), (H2a)
    r"|\bR[12]\b"                                                    # regime names R1 / R2
    r"|\bExperiments?~?\s*\d+(?:\s*(?:and|,|--)\s*\d+)*"              # experiment numbers
    r"|\b[A-Za-z]+\d+[A-Za-z0-9]*\b|\b[A-Z]+-\d+\b"                     # ETD2RK, float64, ETD1, SHA-256
    r"|[LH]\^\{?\d\}?"                                                # L^2 / H^1 function spaces
    r"|\((?:c[12]|[a-f])\)")                                          # ablation tags


def _strip_for_numerals(src: str) -> str:
    s = re.sub(r"(?<!\\)%.*", "", src)
    s = re.sub(r"\\(?:providecommand|newcommand|renewcommand)\*?.*", " ", s)  # macro definitions (layout, not prose)
    s = re.sub(r"\\(?:setlength|ifdim|typeout|@ifundefined|newlength|newsavebox|sbox|setcounter)\b.*", " ", s)
    s = re.sub(r"\\(%s)\*?(\[[^\]]*\])?\{[^}]*\}" % _ARG_CMDS, " ", s)
    s = re.sub(r"\\(?:sub)*section\*?\{", " ", s)
    return s


def numerals(src: str, prose: bool = True):
    """Numeric tokens of LaTeX source: 'a\\times10^{b}' as one token, then plain decimals.  With prose=True the
    source is stripped of comments, macro definitions, reference/label/input/figure arguments, name tokens with
    digits, and sub-/superscripts (indices and exponents of symbols such as P_0^{-1}, \\omega_{\\min})."""
    s = _strip_for_numerals(src) if prose else src
    out = _SCI.findall(s)
    s = _SCI.sub(" ", s)
    if prose:
        s = _NAME_TOKENS.sub(" ", s)
        s = re.sub(r"[_^]\{[^{}]*\}", " ", s)
        s = re.sub(r"[_^]\\?[A-Za-z0-9]", " ", s)
    out += _NUM.findall(s)
    return out


def unsourced_numerals(section_src: str, number_values) -> list:
    """Numerals of the generated prose that do not occur among the numerals of section_numbers.json values."""
    allowed = set()
    for v in number_values:
        allowed.update(numerals(str(v), prose=False))
        allowed.add(str(v))
    return sorted({t for t in numerals(section_src) if t not in allowed})


def lint_files(paths):
    out = []
    for p in paths:
        with open(p) as f:
            out += lint_text(f.read(), p)
    return out


if __name__ == "__main__":
    issues = lint_files(sys.argv[1:])
    for x in issues:
        print(x)
    print(f"texlint: {len(issues)} issue(s) in {len(sys.argv) - 1} file(s)")
    sys.exit(1 if issues else 0)
