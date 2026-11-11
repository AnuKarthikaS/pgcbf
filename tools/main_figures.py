"""Render the two main-paper figures from the CSV tables / JSON logs of the current run.
Usage: python tools/main_figures.py   (writes results/figures/main_fig1_ablation.{pdf,png}, main_fig2_budget.{pdf,png}, main_fig_results.{pdf,png}; the paper uses main_fig_results)"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from pgcbf import plotting as pl  # noqa: E402

if __name__ == "__main__":
    out = os.path.join(ROOT, "results", "figures")
    for f in (pl.main_fig1_ablation, pl.main_fig2_budget, pl.main_fig_results):
        print("wrote", *f(ROOT, out))
