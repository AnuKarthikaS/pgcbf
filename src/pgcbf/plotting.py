"""Plotting helpers: raw trajectories only (no smoothing, no clipping, autoscaled axes that
always include the plotted data and the h = 0 safety boundary).  Saves PDF (vector) + PNG.

V3: constrained_layout, no bbox_inches='tight' (the saved canvas is exactly the declared width: 3.5 in for a single
column, or WIDE = 4.4 in for multi-panel figures; review round: WIDE <= 0.62 x 7.16 in, the figure cap of the
one-column IEEEtran build, so every figure is displayed at its natural size (scale 1, 8 pt fonts) in both the
one-column and the two-column build), titles <= 60 characters, LaTeX math axis labels; the tight bounding boxes of
different axes (with titles, labels and legends) must not overlap.  `save` renders the
figure and refuses (raises) if any text / legend / axes artist leaves the canvas; the per-figure report is written to
figures/layout_report.json and checked by tests/test_figures.py."""
from __future__ import annotations

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams.update({"font.size": 8, "axes.grid": True, "grid.alpha": 0.3, "figure.dpi": 150,
                     "pdf.fonttype": 42, "lines.linewidth": 0.8, "axes.titlesize": 8, "legend.fontsize": 8,
                     "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8, "mathtext.fontset": "cm"})

SINGLE, WIDE = 3.5, 4.4  # IEEE column width; multi-panel width <= 0.62 x 7.16 in (one-column figure cap)
MAIN = 6.0  # main-paper figures, full text width of the jmlr (PMLR / L4DC) class
DOUBLE = WIDE  # backwards-compatible alias
MAX_TITLE = 60


def figure(ncols=1, nrows=1, width=3.5, height=2.4):
    """ncols x nrows grid; total width snapped to 3.5 in (if ncols * width fits a column) or WIDE = 4.4 in."""
    W = SINGLE if ncols * width <= SINGLE + 0.1 else WIDE
    fig, ax = plt.subplots(nrows, ncols, figsize=(W, height * nrows), squeeze=False, layout="constrained")
    return fig, ax


def boundary(ax, y=0, label="safe-set boundary"):
    ax.axhline(y, color="k", lw=0.8, ls="--", label=label)


def out_of_canvas(fig, pad=0.5):
    """Render and list every visible artist (figure texts, legends, and each axes with its drawn decorations:
    title, axis labels, tick labels, offset text, legend) whose tight bbox leaves the canvas by more than `pad` px."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, H = fig.canvas.get_width_height()
    bad = []
    arts = [t for t in fig.texts if t.get_text()] + list(fig.legends) + [a for a in fig.axes if a.get_visible() and a.axison]
    for art in arts:
        if not art.get_visible():
            continue
        bb = art.get_tightbbox(r) if hasattr(art, "get_tightbbox") else art.get_window_extent(r)
        if bb is None or (bb.width == 0 and bb.height == 0):
            continue
        if bb.x0 < -pad or bb.y0 < -pad or bb.x1 > W + pad or bb.y1 > H + pad:
            bad.append({"artist": type(art).__name__, "title": getattr(art, "get_title", lambda: "")(),
                        "bbox": [float(bb.x0), float(bb.y0), float(bb.x1), float(bb.y1)], "canvas": [W, H]})
    return bad


def overlapping_axes(fig, pad=0.5):
    """Render and list pairs of visible axes whose tight bboxes (axes + title + labels + tick labels + legend) overlap
    by more than `pad` px in both directions (e.g. a long title running into the neighbouring panel)."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    axs = [a for a in fig.axes if a.get_visible() and a.axison]
    bbs = [a.get_tightbbox(r) for a in axs]
    bad = []
    for i in range(len(axs)):
        for j in range(i + 1, len(axs)):
            b1, b2 = bbs[i], bbs[j]
            ox = min(b1.x1, b2.x1) - max(b1.x0, b2.x0)
            oy = min(b1.y1, b2.y1) - max(b1.y0, b2.y0)
            if ox > pad and oy > pad:
                bad.append({"axes": [i, j], "titles": [axs[i].get_title(), axs[j].get_title()], "overlap_px": [float(ox), float(oy)]})
    return bad


def legend_outside_axes(fig, pad=0.5):
    """Render and list every in-axes legend (no explicit bbox_to_anchor) that leaves its axes frame by more than
    `pad` px, i.e. a legend wider than a narrow panel running over the tick labels (review round, item 1)."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    bad = []
    for a in fig.axes:
        lg = a.get_legend()
        if lg is None or not lg.get_visible() or lg.get_bbox_to_anchor() is not a.bbox:
            continue
        lb, ab = lg.get_window_extent(r), a.get_window_extent(r)
        if lb.x0 < ab.x0 - pad or lb.x1 > ab.x1 + pad or lb.y0 < ab.y0 - pad or lb.y1 > ab.y1 + pad:
            bad.append({"title": a.get_title(), "legend": [float(lb.x0), float(lb.y0), float(lb.x1), float(lb.y1)],
                        "axes": [float(ab.x0), float(ab.y0), float(ab.x1), float(ab.y1)]})
    return bad


def _split_title(t):
    """Split a title at the blank (outside $...$ math) closest to its middle; None if there is no such blank."""
    depth, cands = 0, []
    for i, ch in enumerate(t):
        if ch == "$":
            depth ^= 1
        elif ch == " " and depth == 0:
            cands.append(i)
    if not cands:
        return None
    i = min(cands, key=lambda k: abs(k - len(t) / 2))
    return t[:i] + "\n" + t[i + 1:]


def wrap_titles(fig):
    """Review round (narrower figures, item 1): break every axes title that is wider than its axes onto two lines
    (at a blank outside math), then let constrained_layout re-solve. Returns the list of wrapped titles."""
    wrapped = []
    for _ in range(2):
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        changed = False
        for a in fig.axes:
            t = a.get_title()
            if not t or not a.get_visible():
                continue
            tw = a.title.get_window_extent(r).width
            if tw > a.get_window_extent(r).width:
                lines = t.split("\n")
                k = max(range(len(lines)), key=lambda j: len(lines[j]))
                new = _split_title(lines[k])
                if new is None:
                    continue
                lines[k] = new
                a.set_title("\n".join(lines))
                wrapped.append(t)
                changed = True
        if not changed:
            break
    return wrapped


def save(fig, outdir, name, widths=(SINGLE, WIDE)):
    os.makedirs(outdir, exist_ok=True)
    wrap_titles(fig)
    titles = [a.get_title().replace("\n", " ") for a in fig.axes if a.get_title()]
    long = [t for t in titles if len(t) > MAX_TITLE]
    bad = out_of_canvas(fig)
    ovl = overlapping_axes(fig)
    lgo = legend_outside_axes(fig)
    wi, hi = fig.get_size_inches()
    rep_path = os.path.join(outdir, "layout_report.json")
    rep = json.load(open(rep_path)) if os.path.exists(rep_path) else {}
    rep[name] = {"width_in": float(wi), "height_in": float(hi), "out_of_canvas": bad, "long_titles": long,
                 "overlapping_axes": ovl, "legend_outside_axes": lgo}
    with open(rep_path, "w") as f:
        json.dump(rep, f, indent=1)
    if bad or long or ovl or lgo or min(abs(wi - w) for w in widths) > 1e-9:
        raise AssertionError(f"figure {name}: layout defect {rep[name]}")
    paths = []
    for ext in ("pdf", "png"):
        p = os.path.join(outdir, f"{name}.{ext}")
        fig.savefig(p)
        paths.append(p)
    plt.close(fig)
    return paths


# ----------------------------------------------------------------------------------------------------------------------
# The two main-paper figures (6.0 in wide, sans-serif, Okabe-Ito colourblind-safe palette, no titles).
# Every plotted number is read from the CSV tables / JSON logs of the run; nothing below is a result value.
# ----------------------------------------------------------------------------------------------------------------------
OKABE_ITO = {"blue": "#0072B2", "orange": "#E69F00", "vermillion": "#D55E00", "green": "#009E73", "sky": "#56B4E9",
             "yellow": "#F0E442", "purple": "#CC79A7", "black": "#000000", "grey": "#999999"}
MAIN_RC = {"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"], "mathtext.fontset": "dejavusans",
           "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "legend.fontsize": 7.5, "xtick.labelsize": 7.5,
           "ytick.labelsize": 7.5, "axes.grid": True, "grid.alpha": 0.3, "pdf.fonttype": 42}


def _read_csv(path):
    import csv
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def _col(row, needle):
    """Value of the first column whose header contains `needle` (headers carry LaTeX decorations)."""
    return next(v for k, v in row.items() if needle in k)


def _boolcell(v):
    return {"True": True, "False": False}.get(v, None)


def _sig(x, d=3):
    return f"{x:.{d}g}"


def _ablation_rows(root):
    """Rows of the ablation panel (deployed R2 + ablations (a), (b), (c), (e)) and the seed count, read from the run."""
    tab, logs = os.path.join(root, "results", "tables"), os.path.join(root, "results", "logs")
    abl = _read_csv(os.path.join(tab, "exp6_ablations.csv"))
    e4 = next(r for r in _read_csv(os.path.join(tab, "exp4_theorem1.csv")) if r["config"] == "R2" and r["start"] == "boundary")
    hyp = json.load(open(os.path.join(logs, "exp1_constants.json")))["data"]["hypotheses"]["R2"]
    cfg = json.load(open(os.path.join(root, "config", "default.json")))
    TsV = float(cfg["experiments"]["single_rate_Ts_violation"])
    n_seed = int(e4["seeds"])
    ex_head = next(k for k in abl[0] if k.startswith("exits"))
    rows = [{"label": "deployed R2", "naive": None, "full": bool(hyp["theorem_valid"]),
             "exits": 0 if float(_col(e4, "min_t h")) >= 0 else None, "minh": float(_col(e4, "min_t h"))}]
    for r in abl:
        nt, nc = r["naive $\\delta_{\\rm total}$"], r["naive $\\delta_{\\rm cert}$"]
        lab = r["ablation"].split(",")[0]
        if lab.startswith("(c)"):
            lab = "(c) single rate\n" + f"(simulated at $T_s={TsV:g}$)"
        rows.append({"label": lab, "naive": (float(nt) <= float(nc)) if nt and nc else None,
                     "full": not _boolcell(r["refused"]), "exits": int(r[ex_head]), "minh": float(r["$\\min_t h$"])})
    assert rows[0]["exits"] == 0, "deployed configuration: min_t h < 0 on the boundary start"
    return rows, n_seed


def _draw_ablation(at, ae, ah, rows, n_seed, cx=None, all_value_labels=True):
    """Draw the verdict table (at), the exit counts (ae) and min_t h (ah) for the rows of _ablation_rows.
    cx: x positions of the label / naive / full columns in the text axis; all_value_labels=False prints only the
    (otherwise invisible) non-negative minima next to the min_t h bars."""
    n = len(rows)
    y = list(range(n))[::-1]
    if True:  # (kept as one block so that the body matches main_fig1_ablation's original layout code)
        at.set_xlim(0, 1)
        at.set_ylim(-0.6, n - 0.4)
        at.axis("off")
        cx = cx or {"label": 0.0, "naive": 0.7, "full": 0.9}
        at.text(cx["naive"], n - 0.45, "naive\ncert.", ha="center", va="bottom", fontweight="bold")
        at.text(cx["full"], n - 0.45, "full\ncert.", ha="center", va="bottom", fontweight="bold")
        word = {True: ("accepts", OKABE_ITO["blue"]), False: ("refuses", OKABE_ITO["vermillion"]),
                None: ("n/a", OKABE_ITO["grey"])}
        for yi, r in zip(y, rows):
            at.text(cx["label"], yi, r["label"], ha="left", va="center")
            for c in ("naive", "full"):
                w, col = word[r[c]]
                at.text(cx[c], yi, w, ha="center", va="center", color=col, fontweight="bold" if r[c] is not None else "normal")
        cols = [OKABE_ITO["vermillion"] if r["exits"] > 0 else OKABE_ITO["blue"] for r in rows]
        ae.barh(y, [r["exits"] for r in rows], color=cols, height=0.6)
        for yi, r in zip(y, rows):
            ae.text(r["exits"], yi, f" {r['exits']}/{n_seed}", va="center", ha="left")
        ae.set_xlim(0, n_seed * 1.5)
        ae.set_yticks(y, [""] * n)
        ae.set_xlabel(f"seeds leaving the safe set\n(of {n_seed}, full-order plant)")
        ah.barh(y, [r["minh"] for r in rows], color=[OKABE_ITO["vermillion"] if r["minh"] < 0 else OKABE_ITO["blue"] for r in rows],
                height=0.6)
        ah.axvline(0, color="k", lw=0.8, ls="--")
        for yi, r in zip(y, rows):  # print the value, so that a positive (invisible) minimum is readable
            if all_value_labels or r["minh"] >= 0:
                ah.text(min(r["minh"], 0), yi, f"{r['minh']:+.2g} ", va="center", ha="right", fontsize=6.5 if all_value_labels else 6)
        lo = min(r["minh"] for r in rows)
        ah.set_xlim(lo * 1.45, -lo * 0.12)
        ah.set_yticks(y, [""] * n)
        ah.set_xlabel("min$_t\\, h(x(t))$ over seeds\n(full-order plant)")
        for a in (ae, ah):
            a.set_ylim(-0.6, n - 0.4)
            a.grid(axis="y", visible=False)




def main_fig1_ablation(root, outdir):
    """Main-paper figure, ablation panel: deployed certificate and ablations (a), (b), (c), (e).  Per row: the naive
    verdict (the certificate with the ablated term dropped; defined only for (a) and (b)), the full verdict, the
    full-order exits out of the seeds, and min_t h over the seeds.  Sources: exp6_ablations.csv, exp4_theorem1.csv
    (deployed R2, boundary start), exp1_constants.json (deployed verdict), config/default.json (T_s of the (c) run)."""
    rows, n_seed = _ablation_rows(root)
    n = len(rows)
    with plt.rc_context(MAIN_RC):
        fig = plt.figure(figsize=(MAIN, 0.29 * n + 0.6), layout="constrained")
        gs = fig.add_gridspec(1, 3, width_ratios=[1.8, 1, 1])
        at, ae, ah = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
        _draw_ablation(at, ae, ah, rows, n_seed)
        fig.get_layout_engine().set(w_pad=0.12)
        return save(fig, outdir, "main_fig1_ablation", widths=(MAIN,))


def _draw_sweep(ar, w, dtot, dcer, fail, wstar, legend_fontsize=7):
    """Exp 9 cutoff-weight sweep: delta_total and delta_cert against w with the (H4) and (H6)-only failure bands."""
    # right: w-sweep with failure bands at the geometric midpoints between grid values
    import math
    order = sorted(range(len(w)), key=lambda i: w[i])
    wl = [w[i] for i in order]
    edges = [wl[0] / math.sqrt(wl[1] / wl[0])] + [math.sqrt(wl[i] * wl[i + 1]) for i in range(len(wl) - 1)] + \
            [wl[-1] * math.sqrt(wl[-1] / wl[-2])]
    seen = set()
    for k, i in enumerate(order):
        f = fail[i]
        if "H4" in f:
            key, col, lab = "H4", OKABE_ITO["vermillion"], "(H4) fails"
        elif f:
            key, col, lab = "H6", OKABE_ITO["orange"], "only (H6) fails"
        else:
            continue
        ar.axvspan(edges[k], edges[k + 1], color=col, alpha=0.18, lw=0, label=None if key in seen else lab)
        seen.add(key)
    ar.plot(w, dtot, "o-", color=OKABE_ITO["blue"], ms=3, label=r"$\delta_{\mathrm{total}}$")
    pos = [(wi, di) for wi, di in zip(w, dcer) if di > 0]
    ar.plot([p[0] for p in pos], [p[1] for p in pos], "s-", color=OKABE_ITO["black"], ms=3,
            label=r"$\delta_{\mathrm{cert}}$ (where $>0$)")
    ar.axvline(wstar, color=OKABE_ITO["green"], lw=1.2, ls="--", label=r"$w^\star$ (deployed)")
    ar.set_xscale("log")
    ar.set_yscale("log")
    ar.set_xlim(edges[0], edges[-1])
    ar.set_xlabel(r"cutoff-mode weight $w$")
    ar.set_ylabel(r"$\delta$ (units of $h$)")
    ar.legend(loc="upper left", fontsize=legend_fontsize, frameon=True, framealpha=0.9, handlelength=1.4)


def main_fig2_budget(root, outdir):
    """Main-paper Fig. 2.  Left: the six parts of delta_total (eq:delta_total_exact; quantization split into the CROWN
    part g_h eps_CROWN and the LUT-rounding part g_h w_2 (eps_A - eps_CROWN)) stacked for R1 and R2, against delta_cert.
    Right: Exp 9 cutoff-weight sweep, delta_total and delta_cert against w with the (H4) and (H6)-only failure bands and
    w* marked.  Sources: exp1_budget.csv, exp9_wsweep_deltas.csv."""
    tab = os.path.join(root, "results", "tables")
    bud = {r["config"]: r for r in _read_csv(os.path.join(tab, "exp1_budget.csv"))}
    parts = ["tail", "LUT", "CROWN", "lag", "hold", "adversary"]
    pcol = dict(zip(parts, [OKABE_ITO[k] for k in ("blue", "yellow", "sky", "purple", "orange", "green")]))
    w, dtot, dcer, fail, wstar = _sweep_data(root)
    with plt.rc_context(MAIN_RC):
        fig, (al, ar) = plt.subplots(1, 2, figsize=(MAIN, 2.25), layout="constrained", width_ratios=[1, 1.35])
        names = ["R1", "R2"]
        for xi, nme in enumerate(names):
            b = bud[nme]
            bottom = 0 * float(b["delta_total"])
            for p in parts:
                v = float(b[p])
                al.bar(xi, v, bottom=bottom, width=0.5, color=pcol[p], edgecolor="white", lw=0.4, label=p if xi == 0 else None)
                bottom += v
            dc = float(b["delta_cert"])
            al.hlines(dc, xi - 0.35, xi + 0.35, color="k", lw=1.6, label=r"$\delta_{\mathrm{cert}}$" if xi == 0 else None)
            al.text(xi, max(dc, bottom), "\n" + r"$\delta_{\mathrm{total}}/\delta_{\mathrm{cert}}$" + f" = {_sig(bottom / dc)}",
                    ha="center", va="bottom")
        al.set_xticks(range(len(names)), names)
        al.set_xlim(-0.85, len(names) - 0.4)
        top = max(float(bud[n]["delta_cert"]) for n in names)
        al.set_ylim(0, top * 1.3)
        al.set_ylabel(r"budget term (units of $h$)")
        al.grid(axis="x", visible=False)
        h, lab = al.get_legend_handles_labels()
        al.legend(h[1:] + h[:1], lab[1:] + lab[:1], loc="upper center", bbox_to_anchor=(1 / 2, -0.16), fontsize=7, ncol=4, frameon=False,
                  handlelength=1.0, columnspacing=1.0)
        _draw_sweep(ar, w, dtot, dcer, fail, wstar)
        return save(fig, outdir, "main_fig2_budget", widths=(MAIN,))


def _sweep_data(root):
    """Exp 9 sweep columns: w, delta_total, delta_cert, failing hypotheses per grid value, and the selected w*."""
    ws = _read_csv(os.path.join(root, "results", "tables", "exp9_wsweep_deltas.csv"))
    w = [float(r["w"]) for r in ws]
    dtot = [float(r["delta_total"]) for r in ws]
    dcer = [float(r["delta_cert"]) for r in ws]
    fail = [[x for x in r["failing"].split() if x != "none"] for r in ws]
    wstar = next(float(r["w"]) for r in ws if r["w_star"] == "True")
    return w, dtot, dcer, fail, wstar


def main_fig_results(root, outdir):
    """Main-paper Figure 1 (single float): (a) the ablation panel of main_fig1_ablation on the left 62% of the width and
    (b) the cutoff-weight sweep of main_fig2_budget (right panel) on the remaining width, in one row.  Same data files
    as the two standalone figures: exp6_ablations.csv, exp4_theorem1.csv, exp1_constants.json, exp9_wsweep_deltas.csv."""
    rows, n_seed = _ablation_rows(root)
    for r in rows:  # the narrower text column: shorten the (c) label, keeping the period it was simulated at
        if r["label"].startswith("(c)"):
            r["label"] = r["label"].replace("(simulated at ", "(at ")
    w, dtot, dcer, fail, wstar = _sweep_data(root)
    n = len(rows)
    with plt.rc_context(MAIN_RC):
        fig = plt.figure(figsize=(MAIN, 0.29 * n + 0.62), layout="constrained")
        gs = fig.add_gridspec(1, 4, width_ratios=[2.6, 1.0, 0.75, 1.9])
        at, ae, ah, ar = (fig.add_subplot(gs[i]) for i in range(4))
        _draw_ablation(at, ae, ah, rows, n_seed, cx={"label": 0.0, "naive": 0.56, "full": 0.9}, all_value_labels=False)
        ae.set_xlim(0, n_seed * 2.25)  # room for the "k/16" labels in the narrower panel
        lo = min(r["minh"] for r in rows)
        ah.set_xlim(lo * 2.0, -lo * 0.12)  # room for the "+k" label of the non-negative minimum
        ae.set_xticks([0, n_seed])
        ae.set_xlabel(f"seeds leaving $\\mathcal{{C}}$\n(of {n_seed}, full plant)")  # shorter labels for the narrow panels
        ah.set_xlabel("min$_t\\, h(x(t))$\n(over seeds)")
        _draw_sweep(ar, w, dtot, dcer, fail, wstar, legend_fontsize=6.5)
        fig.get_layout_engine().set(w_pad=0.14)
        fig.text(0.0, 1.0, "(a)", ha="left", va="top", fontweight="bold")
        fig.text(0.64, 1.0, "(b)", ha="left", va="top", fontweight="bold")
        return save(fig, outdir, "main_fig_results", widths=(MAIN,))
