"""Gate G-A: page-by-page sweep of a PDF for layout/content defects."""
import pymupdf, re, sys
d = pymupdf.open(sys.argv[1]); issues = []
RAW = re.compile(r"\b(eps|theorem|H\d_|delta_|mu_|L_|M_|r_|T_|w\d|naive_|cp_|u_max|b_[AB]|K_m)\w*_\w+\b|\b(True|False|nan|None)\b|\?\?")
for pno, page in enumerate(d, 1):
    W, H = page.rect.width, page.rect.height
    margin = 36
    for img in page.get_image_info():
        r = pymupdf.Rect(img["bbox"])
        if r.x0 < margin-4 or r.x1 > W-margin+4 or r.y1 > H-margin+4: issues.append((pno, "image-outside-type-area", tuple(round(v) for v in r)))
    for dr in page.get_drawings():
        r = dr["rect"]
        if r.x1 > W - 20 or r.x0 < 20: issues.append((pno, "drawing-at-edge", tuple(round(v) for v in r))); break
    td = page.get_text("dict")
    small = 0; raw = []
    for b in td["blocks"]:
        for l in b.get("lines", []):
            for s in l["spans"]:
                if s["size"] < 5.5 and s["text"].strip(): small += 1
                if s["bbox"][2] > W - 20: issues.append((pno, "text-past-right-edge", s["text"][:40]))
                m = RAW.search(s["text"])
                if m and s["font"].lower().find("mono") < 0: raw.append(s["text"].strip()[:50])
    if small: issues.append((pno, "spans-under-5.5pt", small))
    for r in sorted(set(raw))[:12]: issues.append((pno, "raw-identifier-or-flag", r))
for i in issues: print(i)
print("pages", d.page_count, "issues", len(issues))
