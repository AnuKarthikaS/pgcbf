"""Run one command, tee its output to a log, and record wall-clock and peak RSS of the child (os.wait4).

Usage: python tools/timed_run.py <name> <logdir> <cwd> -- cmd args...
Appends "name,exit_code,wall_clock_s,peak_rss_MB" to <logdir>/runtimes.csv and exits with the child's code."""
import os
import subprocess
import sys
import time

name, logdir, cwd = sys.argv[1:4]
cmd = sys.argv[sys.argv.index("--") + 1:]
os.makedirs(logdir, exist_ok=True)
t0 = time.perf_counter()
with open(os.path.join(logdir, f"{name}.log"), "w") as log:
    p = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    for line in p.stdout:
        sys.stdout.write(line)
        log.write(line)
    _, status, ru = os.wait4(p.pid, 0)
wall = time.perf_counter() - t0
code = os.waitstatus_to_exitcode(status)
peak = ru.ru_maxrss / 1024
with open(os.path.join(logdir, "runtimes.csv"), "a") as f:
    f.write(f"{name},{code},{wall:.2f},{peak:.1f}\n")
print(f"### {name}: exit {code}, wall {wall:.2f} s, peak RSS {peak:.1f} MB", flush=True)
sys.exit(code)
