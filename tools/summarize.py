"""Collect the final run: every assertion (pass/fail), runtimes, peak RAM -> results/summary.md + tables.
Exits nonzero if any assertion failed, any script failed, or the suite exceeded the 20-minute budget."""
import csv
import re
import glob
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOG, TAB = os.path.join(ROOT, "results", "logs"), os.path.join(ROOT, "results", "tables")
BUDGET_S = 20 * 60

rows, n_fail = [], 0
for f in sorted(glob.glob(os.path.join(LOG, "exp*.json"))):
    d = json.load(open(f))
    for a in d["assertions"]:
        rows.append([d["experiment"], a["assertion"], "PASS" if a["passed"] else "FAIL"])
        n_fail += not a["passed"]
with open(os.path.join(TAB, "assertions.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["experiment", "assertion", "result"])
    w.writerows(rows)
rt = list(csv.reader(open(os.path.join(LOG, "runtimes.csv"))))[1:]
_m = re.search(r"(\d+) passed", open(os.path.join(LOG, "pytest.log")).read()) if os.path.exists(os.path.join(LOG, "pytest.log")) else None
N_TESTS = int(_m.group(1)) if _m else 0
total = [r for r in rt if r[0] == "TOTAL"]
scripts = [r for r in rt if r[0] != "TOTAL"]
script_fail = [r[0] for r in scripts if r[1] != "0"]
total_s = float(total[-1][2]) if total else sum(float(r[2]) for r in scripts)
peak = max(float(r[3]) for r in scripts)
with open(os.path.join(TAB, "runtimes.tex"), "w") as fh:
    fh.write("\\begin{table}[t]\n\\centering\n\\caption{Measured wall-clock and peak resident memory per script (2 cores); experiment assertions passed: "
             f"${len(rows) - n_fail}$ of ${len(rows)}$; pytest: ${N_TESTS}$ passed.}}\n"
             "\\label{tab:runtimes}\n\\footnotesize\n\\begin{tabular}{lrr}\n\\toprule\nscript & wall [s] & peak RSS [MB] \\\\\n\\midrule\n")
    for r in scripts:
        fh.write(f"{r[0].replace('_', ' ')} & ${float(r[2]):.1f}$ & ${float(r[3]):.0f}$ \\\\\n")
    fh.write(f"\\midrule\ntotal & ${total_s:.1f}$ & ${peak:.0f}$ \\\\\n\\bottomrule\n\\end{{tabular}}\n\\end{{table}}\n")
with open(os.path.join(TAB, "runtimes.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["script", "wall_s", "peak_rss_MB"])
    w.writerows([[r[0], r[2], r[3]] for r in scripts] + [["total", total_s, peak]])
with open(os.path.join(ROOT, "results", "summary.md"), "w") as fh:
    fh.write(f"# Final run summary\n\nAssertions: {len(rows) - n_fail}/{len(rows)} passed in experiments; "
             f"script failures: {script_fail or 'none'}.\n\nTotal wall-clock: {total_s:.1f} s (budget {BUDGET_S} s); "
             f"peak RSS of any script: {peak:.0f} MB.\n\n| script | exit | wall [s] | peak RSS [MB] |\n|---|---|---|---|\n")
    for r in scripts:
        fh.write(f"| {r[0]} | {r[1]} | {r[2]} | {r[3]} |\n")
    fh.write("\n| experiment | assertion | result |\n|---|---|---|\n")
    for r in rows:
        fh.write(f"| {r[0]} | {r[1].replace('|', '/')} | {r[2]} |\n")
print(f"summary: {len(rows) - n_fail}/{len(rows)} experiment assertions passed; script failures {script_fail}; "
      f"total {total_s:.1f} s; peak RSS {peak:.0f} MB")
assert n_fail == 0 and not script_fail and total_s < BUDGET_S
