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

 Other options:  --assay "Variance stabilised"  --pca-colour group  --pca-shape replicate  --outlier-variable lane
                 --volcano-xlim-shrink 1  --volcano-xlim-noshrink 2  --volcano-ymax auto|<number>|none
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
    """All Parameter/Value tables of the report merged (the methods sections have small ones,
    the appendix has the full list)."""
    out = {}
    for w in re.finditer(r'<script type="application/json" data-for="htmlwidget-[^"]+">(\{"x":\{"tag".*?)</script>', doc, flags=re.S):
        try:
            d = json.loads(w.group(1))["x"]["tag"]["attribs"]["data"]
        except Exception:
            continue
        if isinstance(d, dict) and "Parameter" in d and len(d) > 1:
            vcol = [k for k in d if k != "Parameter"][0]
            out.update(zip(d["Parameter"], d[vcol]))
    return out


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
def sample_table(doc):
    """The report's sample table (section "Samples") as {column: values}, or {}."""
    sec, _ = get_section(doc, r"^Samples$", 2)
    if sec is None:
        return {}
    for w in widgets(sec):
        try:
            d = w["x"]["tag"]["attribs"]["data"]
        except (KeyError, TypeError):
            continue
        if isinstance(d, dict):
            return d
    return {}


def draw_pca(ax, pane, colour_by, shape_by=None, meta=None):
    """PCA with colour by `colour_by` ("sample" or a variable, e.g. group) and marker shape by
    `shape_by` (a variable, e.g. replicate, or None). Variables are taken from the report's PCA
    or, failing that, from its sample table (column names are matched case-insensitively)."""
    x = plotly_2d(pane)
    lay = x["layout"]
    meta = meta or {}
    traces = [t for t in x["data"] if t.get("legendgroup")]
    pcaVars = list(dict.fromkeys(t["legendgroup"] for t in traces))

    # sample -> coordinates, in report order
    pts = {}
    for t in traces:
        if t["legendgroup"] == pcaVars[0]:
            for xx, yy, smp in zip(t["x"], t["y"], t.get("text", [])):
                pts[smp] = (xx, yy)

    def values_of(var):
        """sample -> level of var, and report colours per level if the PCA has that variable."""
        if var is None:
            return None, {}
        if var.lower() == "sample":
            return {smp: smp for smp in pts}, {}
        for v in pcaVars:
            if v.lower() == var.lower():
                vals, cols = {}, {}
                for t in traces:
                    if t["legendgroup"] == v:
                        cols[t["name"]] = color(t["marker"].get("color"))
                        for smp in t.get("text", []):
                            vals[smp] = t["name"]
                return vals, cols
        smpCol = next((k for k in meta if k.lower() == "sample"), None)
        col = next((k for k in meta if k.lower() == var.lower()), None)
        if smpCol and col:
            return dict(zip(meta[smpCol], [str(v) for v in meta[col]])), {}
        print("  ! %r not found in the report (PCA: %s; sample table: %s)"
              % (var, ", ".join(pcaVars), ", ".join(meta)))
        return None, {}

    cvals, ccols = values_of(colour_by)
    svals, _ = values_of(shape_by)
    if cvals is None:
        cvals, ccols, colour_by = values_of(pcaVars[0]) + (pcaVars[0],)
    clev = list(dict.fromkeys(cvals[s_] for s_ in pts if s_ in cvals))
    slev = list(dict.fromkeys(svals[s_] for s_ in pts if s_ in svals)) if svals else []
    cmap = plt.get_cmap("tab20" if len(clev) > 10 else "tab10")
    colOf = {l: ccols.get(l, cmap(k % cmap.N)) for k, l in enumerate(clev)}
    markers = ["o", "^", "s", "D", "v", "P", "X", "*"]
    shpOf = {l: markers[k % len(markers)] for k, l in enumerate(slev)}

    for smp, (xx, yy) in pts.items():
        ax.scatter([xx], [yy], s=22, color=colOf.get(cvals.get(smp), "grey"),
                   marker=shpOf.get(svals.get(smp), "o") if svals else "o",
                   edgecolor="k", linewidth=0.4, zorder=3)

    # one legend: colour levels, then shape levels
    blank = Line2D([], [], ls="", marker="")
    handles, labels = [blank], [colour_by.capitalize()]
    for l in clev:
        handles.append(Line2D([], [], ls="", marker="o", color=colOf[l], markeredgecolor="k", markeredgewidth=0.4, ms=4))
        labels.append(l)
    if slev:
        handles.append(blank)
        labels.append(shape_by.capitalize())
        for l in slev:
            handles.append(Line2D([], [], ls="", marker=shpOf[l], color="grey", markeredgecolor="k", markeredgewidth=0.4, ms=4))
            labels.append(l)
    many = len(labels) > 8
    leg = ax.legend(handles, labels, frameon=False, loc="best", fontsize=4.3 if many else 5.5,
                    ncol=2 if len(labels) > 12 else 1, handletextpad=0.2, borderaxespad=0.2)
    for txt, lab in zip(leg.get_texts(), labels):
        if lab in (colour_by.capitalize(), (shape_by or "").capitalize()):
            txt.set_fontweight("bold")

    ax.axhline(0, color="#dddddd", lw=0.5, zorder=0)
    ax.axvline(0, color="#dddddd", lw=0.5, zorder=0)
    ax.set_xlabel(title_of(lay["xaxis"]["title"]))
    ax.set_ylabel(title_of(lay["yaxis"]["title"]))
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
    # d3heatmap stores both the values and the image pixels row by row (variable x PC),
    # whatever width/height the png header claims
    vals = np.array([float(v) for v in mx["data"]]).reshape((nr, nc))
    img = np.asarray(Image.open(io.BytesIO(base64.b64decode(h["image"].split(",")[1]))).convert("RGB"))
    img = img.reshape(-1, 3)[: nr * nc].reshape((nr, nc, 3))
    ax.imshow(img, aspect="auto", interpolation="nearest")
    for i in range(nr):
        for j in range(nc):
            lum = (0.299 * img[i, j, 0] + 0.587 * img[i, j, 1] + 0.114 * img[i, j, 2]) / 255
            ax.text(j, i, "%.2g" % vals[i, j], ha="center", va="center", fontsize=5,
                    color="white" if lum < 0.5 else "black")
    ax.set_xticks(range(nc))
    ax.set_xticklabels([htmllib.unescape(c) for c in mx["cols"]], rotation=55, ha="right", fontsize=5)
    ax.set_yticks(range(nr))
    ax.set_yticklabels(mx["rows"])
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    panel_label(ax, "C", "PC / metadata association (p value; red = most significant)")


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


def auto_ymax(ys, max_out=5):
    """Cap for the y axis that leaves out only isolated extreme points: the top point is
    dropped while it is more than twice the next highest (at most max_out points).
    Returns None if nothing needs capping."""
    v = sorted((y for y in ys if y is not None and np.isfinite(y)), reverse=True)
    k = 0
    while k < min(max_out, len(v) - 1) and v[k] > 2 * v[k + 1]:
        k += 1
    return 1.1 * v[k] if k > 0 else None


def draw_volcano(ax, vol, cnt, xlim=None, ymax="auto"):
    """xlim: symmetric x limit (None = the report's range); ymax: y cap, "auto" (only
    isolated extreme points), a number, or None. Points beyond the axes are drawn at the
    edge as open triangles; the ones above the y cap are labelled with their gene ID."""
    lay = vol["layout"]
    pts = [t for t in vol["data"] if t.get("mode") == "markers" and isinstance(t.get("x"), list)
           and t["x"] and t.get("name", "") != "unselected rows"]
    if xlim is None:
        xr = lay["xaxis"].get("range")
    else:
        xr = [-xlim, xlim]
    if ymax == "auto":
        ymax = auto_ymax([y for t in pts for y in t["y"]])

    clipped = False
    for t in pts:
        nm = t.get("name", "")
        col = color(t["marker"]["color"])
        x = np.array([np.nan if v is None else v for v in t["x"]], float)
        y = np.array([np.nan if v is None else v for v in t["y"]], float)
        outx = (np.abs(x) > xr[1]) if xr else np.zeros(len(x), bool)
        outy = (y > ymax) if ymax else np.zeros(len(y), bool)
        inside = ~(outx | outy)
        z = 2 if nm == "Not significant" else 3
        ax.scatter(x[inside], y[inside], s=3 if nm == "Not significant" else 5, color=col,
                   linewidth=0, label=nm, rasterized=True, zorder=z)
        if (outx | outy).any():
            clipped = True
            for k in np.where(outx | outy)[0]:
                xc = np.clip(x[k], xr[0] * 0.985, xr[1] * 0.985) if xr else x[k]
                yc = min(y[k], ymax * 0.985) if ymax else y[k]
                mk = "^" if outy[k] else (">" if x[k] > 0 else "<")
                ax.scatter([xc], [yc], s=14, marker=mk, facecolors="none", edgecolors=col, linewidths=0.7, zorder=4)
                if outy[k]:
                    lab = t.get("text", [None] * len(x))[k]
                    ax.annotate("%s (%.0f)" % (lab, y[k]) if lab else "%.0f" % y[k], (xc, yc), xytext=(4, -12),
                                textcoords="offset points", fontsize=4.5, color="#333333", va="top",
                                ha="left" if xc < 0 else "right")
    for t in vol["data"]:
        if t.get("mode") == "lines" and isinstance(t.get("x"), list):
            ln = t.get("line", {})
            nm = t.get("name", "")
            ls = {"dash": "--", "dot": ":", "solid": "-"}.get(ln.get("dash", "solid"), "-")
            if len(set(t["y"])) == 1:   # horizontal threshold line: span the whole axis
                ax.axhline(t["y"][0], color="k", lw=0.6, ls=ls, label=nm if nm.startswith("padj") else None, zorder=1)
            else:
                ax.plot(t["x"], t["y"], color="k", lw=0.6, ls=ls,
                        label=nm if nm.startswith("padj") else None, zorder=1)
    if xr:
        ax.set_xlim(*xr)
    if ymax:
        ax.set_ylim(0, ymax)
    else:
        ax.set_ylim(bottom=0)
    xt = title_of(lay["xaxis"]["title"])
    m = re.match(r"higher in (.+?)\s+<<.*>>\s+higher in (.+)$", xt.strip())
    ax.set_xlabel("log2 fold change")
    if m:
        ax.annotate("\u2190 higher in %s" % m.group(1), (0, -0.2), xycoords="axes fraction", ha="left", va="top", fontsize=6)
        ax.annotate("higher in %s \u2192" % m.group(2), (1, -0.2), xycoords="axes fraction", ha="right", va="top", fontsize=6)
    ax.set_ylabel("-log10 adjusted p")
    red, blue = "#b2182b", "#2166ac"
    multicolour_text(ax, "upper left", (0.02, 0.985),
                     [("\u25BC %d down " % cnt["down_fc"], red), ("(%d)" % cnt["down"], blue)])
    multicolour_text(ax, "upper right", (0.98, 0.985),
                     [("%d up " % cnt["up_fc"], red), ("(%d)" % cnt["up"], blue), (" \u25B2", red)])
    ax.legend(frameon=True, framealpha=0.9, edgecolor="none", loc="upper center", bbox_to_anchor=(0.5, 0.9),
              markerscale=1.5, handletextpad=0.2)
    panel_label(ax, "F", "Volcano plot (adjusted p)" + ("; open triangles = beyond axis" if clipped else ""))


# --------------------------------------------------------------------------- main
def pages_for_report(doc, pdf, assay, pca_colour, outlier_var, label=None, pca_shape=None, volcano=None):
    """Add one page per contrast of one report to an open PdfPages. Returns number of pages."""
    params = parameters(doc)
    meta = sample_table(doc)
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

        draw_pca(axA, pick_pane(pca_sec, assay), pca_colour, pca_shape, meta)
        draw_scree(axB, pick_pane(scree_sec, assay))
        draw_pc_heatmap(axC, pcmeta_sec)
        draw_cluster(axD, pick_pane(clus_sec, assay))
        draw_outliers(axE, pick_pane(out_sec, outlier_var), mad, outlier_var)
        volcano = volcano or {}
        shrunk = str(params.get("deseq2_shrink_lfc", "")).upper() in ("TRUE", "T", "1")
        vx = volcano.get("xlim_shrink" if shrunk else "xlim_noshrink")
        draw_volcano(axF, vol, cnt, xlim=vx if vx else None, ymax=volcano.get("ymax", "auto"))

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
    ap.add_argument("--pca-colour", default="group",
                    help="PCA colour: a variable such as group (default) or lane, or sample")
    ap.add_argument("--pca-shape", default="replicate",
                    help="PCA marker shape: a variable such as replicate (default) or lane, or none")
    ap.add_argument("--outlier-variable", default="group", help="outlier tab to use (default: group)")
    ap.add_argument("--volcano-xlim-shrink", type=float, default=1.0,
                    help="volcano x axis +/- limit for runs with LFC shrinkage (0 = report's range; default %(default)s)")
    ap.add_argument("--volcano-xlim-noshrink", type=float, default=2.0,
                    help="volcano x axis +/- limit for runs without shrinkage (0 = report's range; default %(default)s)")
    ap.add_argument("--volcano-ymax", default="auto",
                    help="volcano y axis cap: auto (only isolated extreme points are capped; default), "
                         "a number, or none")
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

    ym = a.volcano_ymax.lower()
    volcano = dict(xlim_shrink=a.volcano_xlim_shrink, xlim_noshrink=a.volcano_xlim_noshrink,
                   ymax=None if ym == "none" else ("auto" if ym == "auto" else float(ym)))
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
                pages += pages_for_report(doc, pdf, a.assay, a.pca_colour, a.outlier_variable, label,
                                          None if a.pca_shape.lower() == "none" else a.pca_shape, volcano)
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
