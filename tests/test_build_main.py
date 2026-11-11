"""tools/build_main.py (story bible Section 9): placeholder resolution from section_numbers.json, failure on an
unresolved placeholder, on an unused new key and on an unsourced prose numeral."""
import json
import os
import subprocess
import sys

from conftest import ROOT

TOOL = os.path.join(ROOT, "tools", "build_main.py")


def _run(tmp_path, template, numbers, extra=()):
    nj = tmp_path / "numbers.json"
    nj.write_text(json.dumps({k: {"tex": v, "raw": None} for k, v in numbers.items()}))
    tp = tmp_path / "tpl.tex"
    tp.write_text(template)
    out = tmp_path / "out.tex"
    r = subprocess.run([sys.executable, TOOL, "--template", str(tp), "--out", str(out), "--numbers", str(nj),
                        "--also", str(tp), *extra], capture_output=True, text=True)
    return r.returncode, r.stdout, out.read_text() if out.exists() else ""


def test_resolves_from_json(tmp_path):
    code, log, out = _run(tmp_path, r"The ratio is \({{util_R2}}\) and \({{e4b_x}}\).", {"util_R2": "0.931", "e4b_x": "0.25"})
    assert code == 0, log
    assert out == r"The ratio is \(0.931\) and \(0.25\)."


def test_unresolved_placeholder_fails(tmp_path):
    code, log, _ = _run(tmp_path, r"Value \({{missing_key}}\).", {"util_R2": "0.931"}, ("--new-prefixes", ""))
    assert code == 1 and "unresolved placeholder" in log


def test_unused_new_key_fails(tmp_path):
    code, log, _ = _run(tmp_path, r"Value \({{a}}\).", {"a": "0.5", "e4b_unused": "0.25"})
    assert code == 1 and "unused new key 'e4b_unused'" in log


def test_unsourced_numeral_fails(tmp_path):
    code, log, _ = _run(tmp_path, r"Value \({{a}}\) and 0.37 typed by hand.", {"a": "0.5"})
    assert code == 1 and "0.37" in log


def test_build_section_unchanged_interface():
    """build_section.py keeps its own entry point (no arguments) and is not replaced by build_main."""
    src = open(os.path.join(ROOT, "tools", "build_section.py")).read()
    assert "def main():" in src and "experimental_section_template.tex" in src
