#!/usr/bin/env python3
"""
Summarise nf-core/differentialabundance HTML reports as landscape A4 pages in ONE pdf
(one page per paramset / contrast).

Panels:  A  PCA (2D)            D  hierarchical-clustering heatmap (variance stabilised)
         B  scree plot          E  outlier detection (MAD score)
         C  PC / metadata       F  volcano plot (adjusted p) + number of up / down genes
            association heatmap

Everything is read from the self-contained report .html (the plotly / d3heatmap
JSON and the embedded PNG), so no pipeline output folder is needed.
Dependencies: matplotlib, numpy, pillow (pyyaml optional, used to read params.yml).

Usage
-----
 Many runs, from a paramsheet (paramset_name of every entry is one run; reports are
 looked up at <workDir>/diffabund/report/<paramset_name>/study_differentialabundance_report.html):

    python summarise_da_report.py --params params.yml --workdir /Volumes/meister.data/.../my_dataset
    python summarise_da_report.py --params params.yml --workdir DIR -o all_runs.pdf
    python summarise_da_report.py --params params.yml --workdir DIR --only 01_noShrink_minAbund10_minSamples3

 Several paramsheets (pages in the order given):

    python summarise_da_report.py --params params.yml params_noOut.yml --workdir DIR

 A single report:

    python summarise_da_report.py study_differentialabundance_report.html -o summary.pdf

 Other options:  --assay "Variance stabilised"  --pca-colour lane  --outlier-variable lane
                 --report-subdir diffabund/report  --report-name study_differentialabundance_report.html

Runs whose report is missing (or cannot be parsed) are skipped with a warning and listed
at the end.

Up / down counts
----------------
Shown as  "N up (M)"  /  "N down (M)":
  N = the report's own numbers (padj <= q and abs(log2FC) >= the pipeline's fold-change
      threshold, e.g. min fold change 2 => abs(log2FC) >= 1); identical to the report's
      "Differential gene counts" table.
  M = in brackets, all genes with padj <= q and abs(log2FC) > 0.
Both are counted from the volcano plot data.  "Up" = higher in the second group of the
contrast (right of the volcano), "down" = higher in the first group (left).
"""
import argparse
import base64
import html as htmllib
import io
import json
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from PIL import Image

A4_LANDSCAPE = (11.69, 8.27)  # inches
plt.rcParams.update({
    "font.size": 6.5, "axes.titlesize": 7.5, "axes.labelsize": 6.5,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 5.5,
    "axes.linewidth": 0.6, "pdf.fonttype": 42,
})


# --------------------------------------------------------------------------- parsing
def strip_tags(t):
    return htmllib.unescape(re.sub(r"<[^>]+>", "", t)).strip()


def get_section(doc, title_re, level, start=0):
    """Return (html_slice, heading_text) for the first <h{level}> whose text matches."""
    for m in re.finditer(r"<h%d[^>]*>(.*?)</h%d>" % (level, level), doc[start:], flags=re.S):
        text = strip_tags(m.group(1))
        if re.search(title_re, text, flags=re.I):
            a = start + m.start()
            nxt = re.search(r"<h[1-%d][ >]" % level, doc[start + m.end():])
            b = start + m.end() + nxt.start() if nxt else len(doc)
            return doc[a:b], text
    return None, None


def tab_panes(sec):
    """{tab label: pane html} for the (first level of) tabset(s) in a section."""
    labels = {}
    for m in re.finditer(r'<a class="nav-link[^"]*" id="(tabset-\d+-\d+)-tab"[^>]*>(.*?)</a>', sec, flags=re.S):
        labels[m.group(1)] = strip_tags(m.group(2))
    ids = list(labels)
    panes = {}
    for k, tid in enumerate(ids):
        a = sec.find('<div id="%s" class="tab-pane' % tid)
        if a < 0:
            continue
        nxt = [sec.find('<div id="%s" class="tab-pane' % o, a + 10) for o in ids[k + 1:]]
        nxt = [n for n in nxt if n > 0]
        panes[labels[tid]] = sec[a:min(nxt)] if nxt else sec[a:]
    return panes


def pick_pane(sec, wanted):
    """Pane for tab `wanted` (case-insensitive substring); whole section if untabbed."""
    panes = tab_panes(sec)
    if not panes:
        return sec
    for lab, p in panes.items():
        if wanted.lower() in lab.lower():
            return p
    raise SystemExit("Tab %r not found; available: %s" % (wanted, list(panes)))


def widgets(pane):
    out = []
    for m in re.finditer(r'<script type="application/json" data-for="(htmlwidget-[^"]+)">(.*?)</script>', pane, flags=re.S):
        out.append(json.loads(m.group(2)))
    return out


def plotly_2d(pane):
    for w in widgets(pane):
        d = w["x"].get("data")
        if d and d[0].get("type") == "scatter":
            return w["x"]
    raise SystemExit("no 2D plotly widget found")


def first_png(pane):
    m = re.search(r"data:image/png;base64,([A-Za-z0-9+/=]+)", pane)
    return Image.open(io.BytesIO(base64.b64decode(m.group(1)))).convert("RGB")


def parameters(doc):
    for w in re.finditer(r'<script type="application/json" data-for="htmlwidget-[^"]+">(\{"x":\{"tag".*?)</script>', doc, flags=re.S):
        try:
            d = json.loads(w.group(1))["x"]["tag"]["attribs"]["data"]
        except Exception:
            continue
        if isinstance(d, dict) and "Parameter" in d:
            vcol = [k for k in d if k != "Parameter"][0]
            return dict(zip(d["Parameter"], d[vcol]))
    return {}


def color(c, default="#888888"):
    if c is None:
        return default
    if isinstance(c, list):
        c = c[0] if c else default
    m = re.match(r"rgba?\(([^)]+)\)", str(c))
    if m:
        v = [float(x) for x in m.group(1).split(",")]
        return (v[0] / 255, v[1] / 255, v[2] / 255, v[3] if len(v) > 3 else 1)
    return c


def title_of(ax):
    t = ax if isinstance(ax, str) else (ax or {}).get("text", "")
    return htmllib.unescape(t or "")


def panel_label(ax, letter, text):
    ax.set_title("%s  %s" % (letter, text), loc="left", fontweight="bold", fontsize=7.5)


# --------------------------------------------------------------------------- panels
def draw_pca(ax, pane, colour_by):
    x = plotly_2d(pane)
    lay = x["layout"]
    for t in x["data"]:
        if t.get("legendgroup", "").lower() != colour_by.lower():
            continue
        ax.scatter(t["x"], t["y"], s=22, color=color(t["marker"].get("color")),
                   edgecolor=color(t["marker"].get("line", {}).get("color"), "k"),
                   linewidth=0.5, label=t["name"], zorder=3)
    ax.axhline(0, color="#dddddd", lw=0.5, zorder=0)
    ax.axvline(0, color="#dddddd", lw=0.5, zorder=0)
    ax.set_xlabel(title_of(lay["xaxis"]["title"]))
    ax.set_ylabel(title_of(lay["yaxis"]["title"]))
    ax.legend(frameon=False, title=colour_by.capitalize(), title_fontsize=6, loc="best", handletextpad=0.2)
    panel_label(ax, "A", "PCA")


def draw_scree(ax, pane):
    x = plotly_2d(pane)
    t = x["data"][0]
    ax.plot(t["x"], t["y"], "o--", color="#1f77b4", ms=3, lw=0.8)
    ax.set_xlabel("PC")
    ax.set_ylabel("Percent variance explained")
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.25, lw=0.4)
    panel_label(ax, "B", "Scree plot")


def draw_pc_heatmap(ax, pane):
    w = [w for w in widgets(pane) if "matrix" in w["x"]]
    if not w:
        ax.text(0.5, 0.5, "PC/metadata heatmap not found", ha="center", transform=ax.transAxes)
        return
    h = w[0]["x"]
    mx = h["matrix"]
    nr, nc = mx["dim"]
    vals = np.array([float(v) for v in mx["data"]]).reshape((nc, nr)).T  # column-major in R
    img = np.array(Image.open(io.BytesIO(base64.b64decode(h["image"].split(",")[1]))).convert("RGB"))
    if img.shape[:2] == (nc, nr):  # d3heatmap stores the image transposed
        img = np.transpose(img, (1, 0, 2))
    ax.imshow(img, aspect="auto", interpolation="nearest")
    for i in range(nr):
        for j in range(nc):
            ax.text(j, i, "%.2g" % vals[i, j], ha="center", va="center", fontsize=5,
                    color="white" if vals[i, j] < 0.05 else "black")
    ax.set_xticks(range(nc))
    ax.set_xticklabels([htmllib.unescape(c) for c in mx["cols"]], rotation=55, ha="right", fontsize=5)
    ax.set_yticks(range(nr))
    ax.set_yticklabels(mx["rows"])
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    panel_label(ax, "C", "PC / metadata association (p value; red = low)")


def draw_cluster(ax, pane):
    img = first_png(pane)
    ax.imshow(img, interpolation="lanczos", aspect="equal")
    ax.set_anchor("N")
    ax.axis("off")
    ax.set_title("D  Hierarchical clustering (variance stabilised)", loc="left", fontweight="bold", fontsize=7.5)


def draw_outliers(ax, pane, mad_thr, variable):
    x = plotly_2d(pane)
    lay = x["layout"]
    pts = [t for t in x["data"] if t.get("mode") == "markers" and t.get("name") in ("TRUE", "FALSE")]
    cats = []
    for t in pts:
        for c in t["x"]:
            if c not in cats:
                cats.append(c)
    rng = np.random.default_rng(1)
    handles = []
    for t in pts:
        xs = np.array([cats.index(c) for c in t["x"]], float) + rng.uniform(-0.15, 0.15, len(t["x"]))
        col = color(t["marker"]["color"])
        ax.scatter(xs, t["y"], s=14, color=col, zorder=3, linewidth=0)
        handles.append(Line2D([], [], marker="o", ls="", color=col, ms=3.5, label=t["name"]))
        if t["name"] == "TRUE":
            for k, (xx, yy, lab) in enumerate(zip(xs, t["y"], t.get("text", []))):
                ax.annotate(lab, (xx, yy), xytext=(5, 6 if k % 2 == 0 else -6), textcoords="offset points",
                            fontsize=4.5, color=col, va="center")
    if mad_thr:
        for s in (-mad_thr, mad_thr):
            ax.axhline(s, color="k", lw=0.6, ls=":")
    ax.axhline(0, color="#cccccc", lw=0.5, zorder=0)
    ax.set_xticks(range(len(cats)))
    ax.set_xticklabels(cats)
    ax.set_xlim(-0.6, len(cats) - 0.4)
    ax.set_xlabel(title_of(lay["xaxis"]["title"]))
    ax.set_ylabel(title_of(lay["yaxis"]["title"]))
    ax.legend(handles=handles, title="Outlier status", title_fontsize=6, frameon=False, loc="lower left")
    panel_label(ax, "E", "Outlier detection (by %s)" % variable)


def volcano_counts(vol):
    """Counts from the volcano traces.
    Returns dict: up/down (padj <= q, abs(lfc) > 0), up_fc/down_fc (padj <= q & abs(lfc) >= fc,
    i.e. the report's own numbers), q and fc (threshold text)."""
    r = dict(up=0, down=0, up_fc=0, down_fc=0, q=None, fc=None)
    for t in vol["data"]:
        nm = t.get("name", "") or ""
        m = re.search(r"Adjusted\s*<=\s*([0-9.eE-]+)", nm)
        if not (m and isinstance(t.get("x"), list)):
            continue
        r["q"] = m.group(1)
        u = sum(1 for v in t["x"] if v > 0)
        d = sum(1 for v in t["x"] if v < 0)
        r["up"] += u
        r["down"] += d
        f = re.search(r"abs\(logFC\)\s*>=\s*([0-9.eE-]+)", nm)
        if f:
            r["fc"] = f.group(1)
            r["up_fc"] += u
            r["down_fc"] += d
    return r


def multicolour_text(ax, loc, xy, pieces, fontsize=8):
    """Text made of differently coloured pieces, anchored inside the axes."""
    from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea
    boxes = [TextArea(t, textprops=dict(color=c, fontsize=fontsize, fontweight="bold")) for t, c in pieces]
    box = AnchoredOffsetbox(loc=loc, child=HPacker(children=boxes, align="baseline", pad=0, sep=0),
                            pad=0, borderpad=0, frameon=False, bbox_to_anchor=xy, bbox_transform=ax.transAxes)
    ax.add_artist(box)


def draw_volcano(ax, vol, cnt):
    lay = vol["layout"]
    for t in vol["data"]:
        nm = t.get("name", "")
        if t.get("mode") == "markers" and isinstance(t.get("x"), list) and t["x"] and nm != "unselected rows":
            ax.scatter(t["x"], t["y"], s=3 if nm == "Not significant" else 5, color=color(t["marker"]["color"]),
                       linewidth=0, label=nm, rasterized=True, zorder=2 if nm == "Not significant" else 3)
        elif t.get("mode") == "lines" and isinstance(t.get("x"), list):
            ln = t.get("line", {})
            ax.plot(t["x"], t["y"], color="k", lw=0.6,
                    ls={"dash": "--", "dot": ":", "solid": "-"}.get(ln.get("dash", "solid"), "-"),
                    label=nm if nm.startswith("padj") else None, zorder=1)
    xr = lay["xaxis"].get("range")
    if xr:
        ax.set_xlim(*xr)
    ax.set_ylim(bottom=0)
    xt = title_of(lay["xaxis"]["title"])
    m = re.match(r"higher in (.+?)\s+<<.*>>\s+higher in (.+)$", xt.strip())
    ax.set_xlabel("log2 fold change")
    if m:
        ax.annotate("← higher in %s" % m.group(1), (0, -0.2), xycoords="axes fraction", ha="left", va="top", fontsize=6)
        ax.annotate("higher in %s →" % m.group(2), (1, -0.2), xycoords="axes fraction", ha="right", va="top", fontsize=6)
    ax.set_ylabel("-log10 adjusted p")
    red, blue = "#b2182b", "#2166ac"
    multicolour_text(ax, "upper left", (0.02, 0.985),
                     [("\u25BC %d down " % cnt["down_fc"], red), ("(%d)" % cnt["down"], blue)])
    multicolour_text(ax, "upper right", (0.98, 0.985),
                     [("%d up " % cnt["up_fc"], red), ("(%d)" % cnt["up"], blue), (" \u25B2", red)])
    ax.legend(frameon=True, framealpha=0.9, edgecolor="none", loc="upper center", bbox_to_anchor=(0.5, 0.9),
              markerscale=1.5, handletextpad=0.2)
    panel_label(ax, "F", "Volcano plot (adjusted p)")


# --------------------------------------------------------------------------- main
def pages_for_report(doc, pdf, assay, pca_colour, outlier_var, label=None):
    """Add one page per contrast of one report to an open PdfPages. Returns number of pages."""
    params = parameters(doc)
    mad = params.get("exploratory_mad_threshold")
    try:
        mad = float(mad)
    except (TypeError, ValueError):
        mad = None

    pca_sec, _ = get_section(doc, r"^Principal components plots", 4)
    scree_sec, _ = get_section(doc, r"^Scree plot", 4)
    pcmeta_sec, _ = get_section(doc, r"metadata associations", 4)
    clus_sec, _ = get_section(doc, r"^Clustering heatmaps", 4)
    out_sec, _ = get_section(doc, r"^Outlier detection", 3)
    det_sec, _ = get_section(doc, r"^Differential gene details", 3)
    for nm, s_ in [("PCA", pca_sec), ("scree", scree_sec), ("PC heatmap", pcmeta_sec), ("clustering", clus_sec),
                   ("outlier", out_sec), ("differential details", det_sec)]:
        if s_ is None:
            raise SystemExit("Section not found in report: " + nm)

    title = strip_tags(re.search(r"<h1[^>]*>(.*?)</h1>", doc, flags=re.S).group(1)) if re.search(r"<h1[^>]*>", doc) else ""
    study = params.get("study_name", "")

    # one page per contrast (= h4 inside "Differential gene details")
    contrasts = []
    for m in re.finditer(r"<h4[^>]*>(.*?)</h4>", det_sec, flags=re.S):
        a = m.start()
        nxt = re.search(r"<h4[ >]", det_sec[m.end():])
        b = m.end() + nxt.start() if nxt else len(det_sec)
        contrasts.append((strip_tags(m.group(1)), det_sec[a:b]))
    if not contrasts:
        raise SystemExit("no contrasts found")

    n = 0
    for cname, csec in contrasts:
        vol = plotly_2d(pick_pane(csec, "Adjusted"))
        cnt = volcano_counts(vol)
        print("  %s\n     padj <= %s & abs(log2FC) >= %s : %d up, %d down\n"
              "     padj <= %s & abs(log2FC) > 0   : %d up, %d down"
              % (cname, cnt["q"], cnt["fc"], cnt["up_fc"], cnt["down_fc"], cnt["q"], cnt["up"], cnt["down"]))

        fig = plt.figure(figsize=A4_LANDSCAPE)
        gs = fig.add_gridspec(3, 3, width_ratios=[3.3, 3.6, 3.9], height_ratios=[3.0, 2.0, 1.55],
                              left=0.045, right=0.985, top=0.835, bottom=0.07, wspace=0.26, hspace=0.62)
        gs_r = gs[:, 2].subgridspec(2, 1, height_ratios=[2.2, 3.6], hspace=0.38)
        axD = fig.add_subplot(gs[:, 0])
        axA = fig.add_subplot(gs[0, 1])
        axB = fig.add_subplot(gs[1, 1])
        axC = fig.add_subplot(gs[2, 1])
        axE = fig.add_subplot(gs_r[0])
        axF = fig.add_subplot(gs_r[1])

        draw_pca(axA, pick_pane(pca_sec, assay), pca_colour)
        draw_scree(axB, pick_pane(scree_sec, assay))
        draw_pc_heatmap(axC, pcmeta_sec)
        draw_cluster(axD, pick_pane(clus_sec, assay))
        draw_outliers(axE, pick_pane(out_sec, outlier_var), mad, outlier_var)
        draw_volcano(axF, vol, cnt)

        # header: paramset name (if given), contrast, details, counts
        y = 0.975
        if label:
            fig.text(0.045, y, label, fontsize=12, fontweight="bold", va="top")
            y -= 0.03
        fig.text(0.045, y, cname, fontsize=10 if label else 12, fontweight="bold", va="top")
        y -= 0.03 if label else 0.028
        sub = "  |  ".join(x for x in [title, study, "PCA / scree / clustering: %s" % assay] if x)
        fig.text(0.045, y, sub, fontsize=7, color="#555555", va="top")
        y -= 0.028
        fig.text(0.045, y,
                 "Differentially expressed genes: %d up, %d down  (adjusted p \u2264 %s, abs(log2FC) > 0)"
                 % (cnt["up"], cnt["down"], cnt["q"]), fontsize=8.5, va="top")
        pdf.savefig(fig)
        plt.close(fig)
        n += 1
    return n


def read_paramset_names(params_file):
    """paramset_name of every entry of a paramsheet (yaml list of mappings), in file order."""
    text = open(params_file, encoding="utf-8").read()
    try:
        import yaml
        data = yaml.safe_load(text)
        if isinstance(data, dict):
            data = [data]
        names = [d["paramset_name"] for d in data if isinstance(d, dict) and "paramset_name" in d]
    except ImportError:
        names = re.findall(r"^[\s-]*paramset_name:\s*['\"]?([^'\"\s#]+)", text, flags=re.M)
    if not names:
        raise SystemExit("no paramset_name found in " + params_file)
    return [str(n) for n in names]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("report", nargs="?", help="a single study_differentialabundance_report.html")
    ap.add_argument("--params", nargs="+", metavar="PARAMS.YML",
                    help="one or more params.yml (paramsheets); every paramset_name is one run, "
                         "pages follow file order then entry order (duplicate names are used once)")
    ap.add_argument("--workdir", help="work directory containing <report-subdir>/<paramset_name>/<report-name>")
    ap.add_argument("--report-subdir", default="diffabund/report", help="default: %(default)s")
    ap.add_argument("--report-name", default="study_differentialabundance_report.html", help="default: %(default)s")
    ap.add_argument("--only", nargs="+", help="only these paramset_name(s)")
    ap.add_argument("-o", "--output", help="output PDF")
    ap.add_argument("--assay", default="Variance stabilised",
                    help="tab used for PCA / scree / clustering heatmap (default: %(default)s)")
    ap.add_argument("--pca-colour", default="group", help="PCA colouring: group (default) or lane/other variable")
    ap.add_argument("--outlier-variable", default="group", help="outlier tab to use (default: group)")
    a = ap.parse_args()

    import os
    jobs = []  # (label, path)
    if a.params:
        if not a.workdir:
            ap.error("--params needs --workdir")
        names = []
        for pf in a.params:
            for n in read_paramset_names(pf):
                if n not in names:
                    names.append(n)
        if a.only:
            names = [n for n in names if n in a.only]
        for n in names:
            jobs.append((n, os.path.join(a.workdir, a.report_subdir, n, a.report_name)))
        out = a.output or os.path.join(a.workdir, "diffabund_summary.pdf")
    elif a.report:
        jobs.append((None, a.report))
        out = a.output or re.sub(r"\.html?$", "", a.report) + "_summary.pdf"
    else:
        ap.error("give a report .html, or --params and --workdir")

    skipped, pages = [], 0
    with PdfPages(out) as pdf:
        for label, path in jobs:
            print("[%s]" % (label or path))
            if not os.path.isfile(path):
                print("  ! report not found, skipped:", path)
                skipped.append((label, "missing report"))
                continue
            try:
                doc = open(path, encoding="utf-8", errors="replace").read()
                pages += pages_for_report(doc, pdf, a.assay, a.pca_colour, a.outlier_variable, label)
            except (SystemExit, Exception) as e:  # keep going with the other runs
                print("  ! could not summarise, skipped:", e)
                skipped.append((label, str(e)))
    if pages == 0:
        try:
            os.remove(out)
        except OSError:
            pass
        raise SystemExit("No pages written.")
    print("\nWrote %s (%d pages from %d of %d runs)" % (out, pages, len(jobs) - len(skipped), len(jobs)))
    for label, why in skipped:
        print("  skipped:", label, "-", why)


if __name__ == "__main__":
    main()
