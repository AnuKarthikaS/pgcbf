#!/usr/bin/env bash
# Full validation suite: every experiment, then the pytest suite.  Exits nonzero if ANY assertion fails.
# Hardware budget: 2 CPU cores (BLAS threads pinned to 2; taskset to 2 cores where available), no GPU.
set -u
cd "$(dirname "$0")"
ROOT="$(pwd)"
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 PYTHONHASHSEED=0 MPLBACKEND=Agg
PY="${PYTHON:-python3}"
PIN=""
if command -v taskset >/dev/null 2>&1; then PIN="taskset -c 0,1"; fi
LOGS="$ROOT/results/logs"
mkdir -p "$LOGS"
echo "name,exit_code,wall_clock_s,peak_rss_MB" > "$LOGS/runtimes.csv"
{
  echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //')"
  echo "cores on machine: $(grep -c ^processor /proc/cpuinfo); used: 2 (BLAS threads=2${PIN:+, $PIN})"
  echo "memory: $(free -m | awk '/Mem/ {print $2" MB"}')"
  $PY -c "import sys,numpy,scipy,matplotlib;print('python',sys.version.split()[0],'numpy',numpy.__version__,'scipy',scipy.__version__,'matplotlib',matplotlib.__version__)"
} > "$LOGS/hardware.txt"
cat "$LOGS/hardware.txt"
# fresh outputs only: no stale table or figure can survive into Section V
rm -f results/tables/*.tex results/tables/*.csv results/figures/*.pdf results/figures/*.png results/figures/layout_report.json results/main.tex \
      results/experimental_section.tex results/section_numbers.json
rm -f "$LOGS"/exp*.json "$LOGS"/exp*.log "$LOGS"/pytest.log
T0=$(date +%s.%N)
FAIL=0
for e in exp0_openloop exp1_constants exp2_crown exp3_delta_cert exp4_theorem1 exp4b_budget exp5_deadzone exp6_ablations exp7_stochastic exp8_numerics exp9_wsweep; do
  $PIN $PY tools/timed_run.py "$e" "$LOGS" "$ROOT/experiments" -- $PY "$e.py" || { echo "!!! $e FAILED"; FAIL=1; }
done
$PIN $PY tools/timed_run.py pytest "$LOGS" "$ROOT" -- $PY -m pytest -q tests -p no:cacheprovider || { echo "!!! pytest FAILED"; FAIL=1; }
T1=$(date +%s.%N)
TOTAL=$(echo "$T1 - $T0" | bc)
echo "TOTAL,$FAIL,$TOTAL," >> "$LOGS/runtimes.csv"
$PY tools/summarize.py || FAIL=1
# V4: Section V prose generated from the logs/CSVs of this run (fails on unresolved placeholders or lint findings)
$PY tools/build_section.py || { echo "!!! build_section FAILED"; FAIL=1; }
# main-paper figures from the CSVs / logs of this run
$PY tools/main_figures.py || { echo "!!! main_figures FAILED"; FAIL=1; }
if [ -f paper/main_template.tex ]; then
  $PY tools/build_main.py --template paper/main_template.tex --out results/main.tex || { echo "!!! build_main FAILED"; FAIL=1; }
fi
# V1 lint of the final tables and Section V
$PY -m pytest -q tests/test_tex_lint.py tests/test_figures.py -p no:cacheprovider > "$LOGS/lint_final.log" 2>&1 \
  && tail -1 "$LOGS/lint_final.log" || { cat "$LOGS/lint_final.log"; echo "!!! final lint FAILED"; FAIL=1; }
# V6: orchestrator red-team gates and the standalone compile + page sweep (G-A, G-I)
$PY tools/redteam_gates.py "$ROOT" > "$LOGS/redteam_gates.log" 2>&1 && tail -1 "$LOGS/redteam_gates.log" \
  || { cat "$LOGS/redteam_gates.log"; echo "!!! redteam_gates FAILED"; FAIL=1; }
$PY tools/standalone_pdf.py || { echo "!!! standalone compile / page sweep FAILED"; FAIL=1; }
T2=$(date +%s.%N)
echo "total wall-clock: $TOTAL s (suite) + $(echo "$T2 - $T1" | bc) s (section build, lint, gates, PDF); status: $([ $FAIL -eq 0 ] && echo PASS || echo FAIL)"
exit $FAIL
