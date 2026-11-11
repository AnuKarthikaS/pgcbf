"""V3: every generated figure uses constrained layout at a fixed width (3.5 in single panel / 4.4 in multi-panel, both
<= 0.62 x 7.16 in, the one-column figure cap), keeps every artist (titles, labels, tick labels, legends) inside the
saved canvas, has no overlapping panels or in-axes legends wider than their panel, and has titles of at most 60 characters.
The layout report is written by pgcbf.plotting.save at render time; the bbox checker itself is tested here."""
import json
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
import matplotlib  # noqa: E402
matplotlib.use("Agg")
from pgcbf import plotting as pl  # noqa: E402

FIGDIR = os.path.join(ROOT, "results", "figures")
REPORT = os.path.join(FIGDIR, "layout_report.json")
TEXTWIDTH_IN, TEXTHEIGHT_IN = 7.16, 9.0  # one-column IEEEtran text block (height: conservative lower bound)
FIG_CAP_W, FIG_CAP_H = 0.62, 0.72  # review item 1: width <= 0.62 textwidth, height <= 0.72 textheight
EXPECTED = ["exp0_openloop_h", "exp2_crown", "exp3_delta_cert", "exp4_theorem1", "exp5_deadzone_sweep", "exp5_deadzone_z",
            "exp6_ablations", "exp7_stochastic", "exp8_convergence", "exp9_wsweep"]


def test_bbox_checker_detects_clipped_artists(tmp_path):
    fig, ax = pl.figure(1, 1, 3.5, 2.4)
    ax[0, 0].plot([0, 1], [0, 1])
    assert pl.out_of_canvas(fig) == []
    fig.text(0.95, 0.5, "a long label that runs off the right edge of the canvas")
    assert pl.out_of_canvas(fig), "text past the right edge must be reported"
    with pytest.raises(AssertionError):
        pl.save(fig, str(tmp_path), "clipped")


def test_save_refuses_long_titles_and_non_ieee_width(tmp_path):
    fig, ax = pl.figure(1, 1, 3.5, 2.4)
    ax[0, 0].set_title("x" * (pl.MAX_TITLE + 1))
    with pytest.raises(AssertionError):
        pl.save(fig, str(tmp_path), "long")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(5, 3), layout="constrained")
    fig.add_subplot(111).plot([0, 1])
    with pytest.raises(AssertionError):
        pl.save(fig, str(tmp_path), "width")


@pytest.mark.skipif(not os.path.exists(REPORT), reason="figures not rendered yet")
def test_every_figure_rendered_clean():
    rep = json.load(open(REPORT))
    for name in EXPECTED:
        assert name in rep, f"missing figure {name}"
        r = rep[name]
        assert r["out_of_canvas"] == [], (name, r["out_of_canvas"])
        assert r["long_titles"] == [], (name, r["long_titles"])
        assert r.get("overlapping_axes", None) == [], (name, r.get("overlapping_axes"))
        assert r.get("legend_outside_axes", None) == [], (name, r.get("legend_outside_axes"))
        assert min(abs(r["width_in"] - pl.SINGLE), abs(r["width_in"] - pl.WIDE)) < 1e-9, (name, r["width_in"])
        assert r["width_in"] <= FIG_CAP_W * TEXTWIDTH_IN and r["height_in"] <= FIG_CAP_H * TEXTHEIGHT_IN, (name, r)
        for ext in ("pdf", "png"):
            assert os.path.exists(os.path.join(FIGDIR, f"{name}.{ext}")), (name, ext)


def test_overlap_checker_detects_colliding_panels(tmp_path):
    fig, ax = pl.figure(2, 1, 3.58, 2.3)
    for a in ax[0]:
        a.plot([0, 1], [0, 1])
    assert pl.overlapping_axes(fig) == []
    ax[0, 0].set_title("a deliberately long title that runs into the next panel", loc="left")
    ax[0, 0].title.set_wrap(False)
    fig.canvas.draw()
    assert pl.overlapping_axes(fig), "a title running into the neighbouring panel must be reported"


def test_legend_checker_detects_legend_wider_than_panel(tmp_path):
    fig, ax = pl.figure(2, 1, 3.58, 2.3)
    a = ax[0, 1]
    a.plot([0, 1], [0, 1], label="x")
    a.legend(loc="best")
    assert pl.legend_outside_axes(fig) == []
    a.plot([0, 1], [1, 0], label="a deliberately long legend entry wider than the panel")
    a.legend(loc="best")
    assert pl.legend_outside_axes(fig), "a legend running over the tick labels must be reported"
    a.legend(loc="upper center", bbox_to_anchor=(1 / 2, -1 / 4))
    assert pl.legend_outside_axes(fig) == [], "legends placed below the axes on purpose are not in-axes legends"


def test_wrap_titles_breaks_titles_wider_than_their_axes(tmp_path):
    fig, ax = pl.figure(2, 1, 3.58, 2.3)
    ax[0, 0].set_title(r"a fairly long panel title with $b_A$ math inside it")
    assert pl.wrap_titles(fig)
    assert "\n" in ax[0, 0].get_title() and ax[0, 0].get_title().count("$") == 2


MAIN_FIGS = ["main_fig1_ablation", "main_fig2_budget", "main_fig_results"]


@pytest.mark.skipif(not (os.path.exists(REPORT) and all(n in json.load(open(REPORT)) for n in MAIN_FIGS)),
                    reason="main-paper figures not rendered yet (tools/main_figures.py)")
def test_main_figures_rendered_clean():
    """The main-paper figures are 6.0 in wide (jmlr text width), clean, untitled, and saved as PDF + PNG."""
    rep = json.load(open(REPORT))
    for name in MAIN_FIGS:
        r = rep[name]
        assert abs(r["width_in"] - pl.MAIN) < 1e-9, (name, r["width_in"])
        for k in ("out_of_canvas", "long_titles", "overlapping_axes", "legend_outside_axes"):
            assert r[k] == [], (name, k, r[k])
        for ext in ("pdf", "png"):
            assert os.path.exists(os.path.join(FIGDIR, f"{name}.{ext}")), (name, ext)


def test_main_figure_code_has_no_titles_and_reads_data_files():
    import inspect
    for f in (pl.main_fig1_ablation, pl.main_fig2_budget, pl.main_fig_results):
        src = inspect.getsource(f)
        assert "set_title" not in src and "suptitle" not in src, f.__name__
        # the figure reads the run's CSV tables, directly or through its _ablation_rows / _sweep_data loaders
        assert "_read_csv(" in src or "_ablation_rows(" in src or "_sweep_data(" in src, f.__name__
