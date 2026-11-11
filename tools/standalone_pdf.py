"""Gates G-A / G-I: quick standalone compile of Section V alone, then the orchestrator's page sweep.

* builds paper_build/section_v.tex (IEEEtran journal, ONE-COLUMN: the layout of the combined manuscript, the primary
  gate since the review round) and paper_build/section_v_twocolumn.tex (IEEEtran journal, two-column), both with
  amsmath/amssymb/booktabs/graphicx and the shorthand macros of
  paper/theory/sections_II_IV.tex and \\newlabel stubs for every label that Section V cites from Sections II-IV /
  Appendices A-D (numbered in the order they appear in the theory files; labels introduced by the theory revision
  spec, e.g. lem:delta_cert_classH, get the next free number),
* compiles with tectonic (or pdflatex if tectonic is absent), twice,
* fails on any LaTeX error, undefined reference or '??' in the PDF text,
* runs tools/page_sweep.py on the PDF and fails if it reports any issue,
* checks every figure size written to the log by the \\pgcfig macro (width <= 0.62 textwidth, height <= 0.72
  textheight; review item 1) and that no text span or drawing reaches the bottom page margin,
* renders every page to results/pdf_pages/<variant>/page-NN.png for visual inspection.
Exit status 1 on any failure."""
import glob
import json
import os
import re
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BUILD = os.path.join(ROOT, "paper_build")
RES = os.path.join(ROOT, "results")
THEORY = sorted(glob.glob(os.path.join(ROOT, "paper", "theory", "*.tex")))
KIND = {"eq": None, "lem": "Lemma", "prop": "Proposition", "thm": "Theorem", "rem": "Remark", "assump": "Assumption",
        "cor": "Corollary", "sec": "Section", "app": "Appendix", "def": "Definition"}


def theory_labels():
    labs = []
    macros = []
    for p in THEORY:
        s = open(p).read()
        labs += re.findall(r"\\label\{([^}]+)\}", s)
        macros += [ln for ln in s.splitlines() if ln.startswith("\\providecommand")]
    return labs, macros


def stubs(cited, internal, labs):
    counters, out = {}, []
    order = [l for l in labs if l in cited] + sorted(l for l in cited if l not in labs)
    roman = ["I", "II", "III", "IV", "V", "VI"]
    for lab in order:
        if lab in internal:
            continue
        pre = lab.split(":")[0]
        counters[pre] = counters.get(pre, 0) + 1
        n = counters[pre]
        if pre == "sec":
            val = roman[min(n, len(roman) - 1)]
        elif pre == "app":
            val = chr(ord("A") + n - 1)
        else:
            val = str(n)
        out.append("\\newlabel{%s}{{%s}{1}}" % (lab, val))
    return out


VARIANTS = (("onecolumn", "section_v", r"\documentclass[journal,onecolumn]{IEEEtran}"),
            ("twocolumn", "section_v_twocolumn", r"\documentclass[journal]{IEEEtran}"))


def _pt(x):
    return float(x.replace("pt", ""))


def compile_variant(variant, stem, docclass, cited, internal, labs, macros):
    import pymupdf
    main_tex = "\n".join([
        docclass,
        r"\usepackage{amsmath,amssymb,booktabs,graphicx}",
        r"\usepackage[T1]{fontenc}",
        r"\makeatletter\input{external_labels.aux}\makeatother",
        *macros,
        r"\begin{document}",
        r"\setcounter{section}{4}",
        r"\input{experimental_section}",
        r"\end{document}", ""])
    with open(os.path.join(BUILD, f"{stem}.tex"), "w") as f:
        f.write(main_tex)
    eng = shutil.which("tectonic")
    if eng:
        cmd = [eng, "--keep-logs", "--chatter", "minimal", f"{stem}.tex"]
        runs = 1
    elif shutil.which("pdflatex"):
        cmd = ["pdflatex", "-interaction=nonstopmode", f"{stem}.tex"]
        runs = 2
    else:
        print("standalone_pdf: no TeX engine (tectonic/pdflatex) found")
        return os.environ.get("PDF_GATE") == "skip"
    for _ in range(runs):
        p = subprocess.run(cmd, cwd=BUILD, capture_output=True, text=True)
    logtxt = p.stdout + p.stderr
    lp = os.path.join(BUILD, f"{stem}.log")
    if os.path.exists(lp):
        logtxt += open(lp, errors="replace").read()
    errors = [l for l in logtxt.splitlines() if l.startswith("!") or "error:" in l.lower()]
    undef = [l for l in logtxt.splitlines() if re.search(r"(Reference|Citation) .* undefined|undefined references", l)]
    overfull = [l for l in logtxt.splitlines() if l.startswith("Overfull \\hbox")]
    pdf = os.path.join(BUILD, f"{stem}.pdf")
    ok = p.returncode == 0 and os.path.exists(pdf) and not errors and not undef
    with open(os.path.join(RES, "logs", f"standalone_compile_{variant}.log"), "w") as f:
        f.write(logtxt)
    print(f"standalone_pdf [{variant}]: engine {os.path.basename(cmd[0])}, exit {p.returncode}, {len(errors)} errors, "
          f"{len(undef)} undefined references, {len(overfull)} overfull hboxes")
    for l in errors[:10] + undef[:10]:
        print("  ", l)
    if not os.path.exists(pdf):
        return False
    # figure-size gate (review item 1): sizes typed out by \pgcfig
    # TeX wraps log lines at 79 characters: un-wrap before parsing the \typeout lines
    flat = logtxt.replace("\n", "")
    figs = re.findall(r"PGCFIG (\S+?\.pdf) wd=([\d.]+pt) ht=([\d.]+pt) tw=([\d.]+pt) th=([\d.]+pt) lw=([\d.]+pt)", flat)
    n_fig_src = len(re.findall(r"\\pgcfig\{", open(os.path.join(BUILD, "experimental_section.tex")).read()))
    fig_bad = [(f, wd, ht) for f, wd, ht, tw, th, lw in figs
               if _pt(wd) > 0.62 * _pt(tw) + 0.5 or _pt(ht) > 0.72 * _pt(th) + 0.5 or _pt(wd) > _pt(lw) + 0.5]
    fig_rows = [{"figure": f, "width_pt": _pt(wd), "height_pt": _pt(ht), "width_cap_pt": 0.62 * _pt(tw),
                 "height_cap_pt": 0.72 * _pt(th)} for f, wd, ht, tw, th, lw in figs]
    with open(os.path.join(RES, "logs", f"figure_sizes_{variant}.json"), "w") as f:
        json.dump(fig_rows, f, indent=1)
    print(f"standalone_pdf [{variant}]: {len(figs)} figures sized by \\pgcfig (source has {n_fig_src}), "
          f"{len(fig_bad)} above the caps (w <= 0.62 textwidth, h <= 0.72 textheight)")
    for b in fig_bad:
        print("   figure above cap:", b)
    sweep = subprocess.run([sys.executable, os.path.join(ROOT, "tools", "page_sweep.py"), pdf], capture_output=True, text=True)
    print(sweep.stdout.strip())
    with open(os.path.join(RES, "logs", f"page_sweep_{variant}.log"), "w") as f:
        f.write(sweep.stdout + sweep.stderr)
    m = re.search(r"issues (\d+)", sweep.stdout)
    n_issues = int(m.group(1)) if m else -1
    outdir = os.path.join(RES, "pdf_pages", variant)
    shutil.rmtree(outdir, ignore_errors=True)
    os.makedirs(outdir)
    doc = pymupdf.open(pdf)
    txt = "".join(pg.get_text() for pg in doc)
    qq = txt.count("??")
    for i, pg in enumerate(doc, 1):
        pg.get_pixmap(dpi=110).save(os.path.join(outdir, f"page-{i:02d}.png"))
    shutil.copy(pdf, os.path.join(RES, "section_v_standalone.pdf" if variant == "onecolumn" else f"section_v_standalone_{variant}.pdf"))
    print(f"standalone_pdf [{variant}]: {doc.page_count} pages, page-sweep issues {n_issues}, '??' occurrences {qq}")
    # extra layout gates (beyond page_sweep): no overfull box wider than 1 pt, no two text spans overlapping, and
    # nothing (text or drawing) in the bottom margin strip (a float taller than the page is pushed past the margin)
    big = [l for l in overfull if float(re.search(r"\(([\d.]+)pt too wide", l).group(1)) > 1]
    overlaps, bottom = [], []
    for i, pg in enumerate(doc, 1):
        H = pg.rect.height
        spans = [(bi, pymupdf.Rect(sp["bbox"])) for bi, b in enumerate(pg.get_text("dict")["blocks"])
                 for ln in b.get("lines", []) for sp in ln["spans"] if len(sp["text"].strip()) > 1]
        bottom += [(i, "text", tuple(round(v) for v in r)) for _, r in spans if r.y1 > H - 18]
        bottom += [(i, "drawing", tuple(round(v) for v in d["rect"])) for d in pg.get_drawings() if d["rect"].y1 > H - 18][:1]
        spans.sort(key=lambda r: r[1].y0)
        for a_i, (ba, ra) in enumerate(spans):
            for bb, rb in spans[a_i + 1:]:
                if rb.y0 >= ra.y1:
                    break
                if ba == bb:
                    continue
                inter = ra & rb
                if not inter.is_empty and inter.get_area() > 0.3 * min(ra.get_area(), rb.get_area()):
                    overlaps.append((i, tuple(round(v) for v in ra), tuple(round(v) for v in rb)))
    for l in big:
        print("   overfull:", l)
    for o in overlaps[:20]:
        print("   overlapping text:", o)
    for o in bottom[:20]:
        print("   bottom margin:", o)
    print(f"standalone_pdf [{variant}]: {len(big)} overfull boxes > 1pt, {len(overlaps)} overlapping text spans, "
          f"{len(bottom)} objects in the bottom margin")
    return (ok and n_issues == 0 and qq == 0 and not big and not overlaps and not bottom and not fig_bad
            and len(figs) == n_fig_src)


def main():
    sec = os.path.join(RES, "experimental_section.tex")
    src = open(sec).read()
    tabs = {os.path.basename(p)[:-4]: open(p).read() for p in glob.glob(os.path.join(RES, "tables", "*.tex"))}
    used = re.findall(r"\\input\{tables/([^}]+)\}", src)
    allsrc = src + "".join(tabs.get(u, "") for u in used)
    cited = set(re.findall(r"\\(?:eq)?ref\{([^}]+)\}", allsrc))
    internal = set(re.findall(r"\\label\{([^}]+)\}", allsrc))
    labs, macros = theory_labels()
    os.makedirs(BUILD, exist_ok=True)
    for d in ("tables", "figures"):
        dst = os.path.join(BUILD, d)
        shutil.rmtree(dst, ignore_errors=True)
        shutil.copytree(os.path.join(RES, d), dst)
    shutil.copy(sec, os.path.join(BUILD, "experimental_section.tex"))
    with open(os.path.join(BUILD, "external_labels.aux"), "w") as f:
        f.write("\\relax\n" + "\n".join(stubs(cited, internal, labs)) + "\n")
    for old in glob.glob(os.path.join(RES, "pdf_pages", "page-*.png")):
        os.remove(old)
    res = {v: compile_variant(v, stem, dc, cited, internal, labs, macros) for v, stem, dc in VARIANTS}
    print("standalone_pdf: " + ", ".join(f"{v} {'PASS' if r else 'FAIL'}" for v, r in res.items()))
    sys.exit(0 if all(res.values()) else 1)


if __name__ == "__main__":
    main()
