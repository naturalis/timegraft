"""Figure 2: overview of timegraft.

Writes fig2.svg next to this script. Convert with
    rsvg-convert -f pdf -o fig2.pdf fig2.svg
    rsvg-convert -d 300 -p 300 -o fig2.png fig2.svg
"""
from pathlib import Path

W, H = 1000, 640
BLUE = "#2f6b9a"    # published syntheses (OToL, TimeTree)
ORANGE = "#d0782a"  # barcodes
GREY = "#9a9a9a"
INK = "#222222"
SOFT = "#f3f4f6"
FONT = "Helvetica, Arial, 'Liberation Sans', sans-serif"

out: list[str] = []


def text(x, y, s, size=11, anchor="start", weight="normal", fill=INK, style="normal"):
    out.append(
        f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" text-anchor="{anchor}" '
        f'font-weight="{weight}" font-style="{style}" fill="{fill}">{s}</text>'
    )


def line(x1, y1, x2, y2, stroke=INK, width=1.6, dash=None, cap="round"):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    out.append(
        f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{stroke}" '
        f'stroke-width="{width}" stroke-linecap="{cap}"{d}/>'
    )


def rect(x, y, w, h, fill=SOFT, stroke="none", rx=6, width=1, dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    out.append(
        f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{rx}" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="{width}"{d}/>'
    )


def dot(x, y, r=3.6, fill=INK, stroke=None):
    s = f' stroke="{stroke}" stroke-width="1.4"' if stroke else ""
    out.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{fill}"{s}/>')


def arrow(x1, y1, x2, y2, stroke=INK, width=1.6):
    out.append(
        f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{stroke}" '
        f'stroke-width="{width}" marker-end="url(#arr)"/>'
    )


def elbow(px, py, cx, cy, stroke, width=1.8, dash=None):
    """Rectangular cladogram branch from parent (px, py) to child (cx, cy)."""
    line(px, py, px, cy, stroke, width, dash, cap="square")
    line(px, cy, cx, cy, stroke, width, dash, cap="butt")


def axis(x_of, ticks, y, label="Ma"):
    line(x_of(max(ticks)), y, x_of(min(ticks)), y, "#666", 1)
    for t in ticks:
        line(x_of(t), y, x_of(t), y + 4, "#666", 1)
        text(x_of(t), y + 15, f"{t:g}", 9.5, "middle", fill="#555")
    text(x_of(min(ticks)) + 8, y + 15, label, 9.5, fill="#555")


def panel(x, y, w, h, title):
    rect(x, y, w, h, fill="#ffffff", stroke="#d5d8dc", rx=8)
    text(x + 12, y + 22, title, 12.5, weight="bold")


# ---------------------------------------------------------------- sources + legend
text(15, 22, "Inputs", 11, weight="bold", fill="#555")
rect(15, 32, 225, 50, fill="#fbeee2", stroke=ORANGE, rx=6)
text(127, 52, "BOLD (COI) or UNITE (ITS)", 11.5, "middle", "bold", ORANGE)
text(127, 69, "reference barcodes, BINs or SHs", 10.5, "middle")
rect(260, 32, 112, 50, fill="#e4eef6", stroke=BLUE, rx=6)
text(316, 52, "Open Tree", 11.5, "middle", "bold", BLUE)
text(316, 69, "topology", 10.5, "middle")
rect(380, 32, 112, 50, fill="#e4eef6", stroke=BLUE, rx=6)
text(436, 52, "TimeTree", 11.5, "middle", "bold", BLUE)
text(436, 69, "divergence times", 10.5, "middle")
arrow(127, 82, 127, 106, ORANGE)
arrow(316, 82, 345, 106, BLUE)
arrow(436, 82, 405, 106, BLUE)

lx, ly = 560, 40
line(lx, ly, lx + 28, ly, BLUE, 2.4)
text(lx + 36, ly + 4, "from published syntheses", 10.5)
line(lx, ly + 18, lx + 28, ly + 18, ORANGE, 2.4)
text(lx + 36, ly + 22, "from barcodes", 10.5)
line(lx, ly + 36, lx + 28, ly + 36, GREY, 1.8, "3 3")
text(lx + 36, ly + 40, "related genera, for rooting only", 10.5)
dot(lx + 230, ly, 4, BLUE)
text(lx + 240, ly + 4, "age from TimeTree", 10.5)
dot(lx + 230, ly + 18, 4, "#fff", BLUE)
text(lx + 240, ly + 22, "age interpolated", 10.5)

PY, PH = 108, 318
PX = [15, 260, 505, 750]
PW = 235

# ---------------------------------------------------------------- (i) reference units
x0 = PX[0]
panel(x0, PY, PW, PH, "(i) Reference units")
text(x0 + 12, PY + 40, "one representative sequence per BIN", 10.5, fill="#444")
seqs = [
    [(0, 150, False), (6, 146, True), (2, 140, False), (10, 150, False)],
    [(0, 148, True), (4, 150, False), (8, 138, False)],
    [(3, 150, False), (0, 144, True), (6, 149, False), (1, 147, False)],
]
yy = PY + 62
for i, group in enumerate(seqs):
    gh = len(group) * 11 + 10
    rect(x0 + 52, yy - 8, 172, gh, fill=SOFT, rx=5)
    text(x0 + 14, yy + gh / 2 - 4, f"BIN {i + 1}", 10.5, fill="#444")
    for j, (off, ln, rep) in enumerate(group):
        sy = yy + 2 + j * 11
        if rep:
            line(x0 + 62 + off, sy, x0 + 62 + off + ln, sy, ORANGE, 3.4)
        else:
            line(x0 + 62 + off, sy, x0 + 62 + off + ln, sy, "#b9bec4", 2.2)
    yy += gh + 8
# a filtered sequence
line(x0 + 62, yy + 2, x0 + 120, yy + 2, "#c9cdd1", 2.2, "4 3")
text(x0 + 128, yy + 6, "× removed", 10.5, fill="#a33")
yy += 24
for s in ("filters: length, ambiguous bases,", "reading frame (COI)",
          "genus = majority identification"):
    text(x0 + 14, yy, s, 10.5, fill="#444")
    yy += 15

# ---------------------------------------------------------------- (ii) dated backbone
x0 = PX[1]
panel(x0, PY, PW, PH, "(ii) Dated backbone")
text(x0 + 12, PY + 40, "one tip per genus, ages in Ma", 10.5, fill="#444")
xt = x0 + 180
xb = lambda a: xt - a * 1.5  # noqa: E731  100 Ma -> 150 px
tips = {f"G{i}": PY + 70 + i * 36 for i in range(5)}
nodes = {  # name: (age, children, source)
    "n45": (30, ["G3", "G4"], "tt"),
    "n12": (40, ["G0", "G1"], "tt"),
    "n345": (65, ["G2", "n45"], "interp"),
    "root": (100, ["n12", "n345"], "tt"),
}
pos = {k: (xt, y) for k, y in tips.items()}


def place(n):
    if n in pos:
        return pos[n]
    age, kids, _ = nodes[n]
    ys = [place(k)[1] for k in kids]
    pos[n] = (xb(age), sum(ys) / len(ys))
    return pos[n]


place("root")
for n, (age, kids, src) in nodes.items():
    px, py = pos[n]
    for k in kids:
        elbow(px, py, *pos[k], BLUE, 2.2)
line(pos["root"][0] - 10, pos["root"][1], pos["root"][0], pos["root"][1], BLUE, 2.2)
for n, (age, kids, src) in nodes.items():
    px, py = pos[n]
    dot(px, py, 4, BLUE) if src == "tt" else dot(px, py, 4, "#fff", BLUE)
for k, y in tips.items():
    focus = k == "G3"
    text(xt + 7, y + 4, "genus" if focus else f"G{int(k[1]) + 1}", 10.5,
         weight="bold" if focus else "normal", fill=ORANGE if focus else INK)
sx, sy = pos["n45"]
text(sx - 5, sy - 7, "S", 11, "end", "bold", BLUE)
yax = PY + 262
axis(xb, [100, 50, 0], yax)
line(sx, sy + 6, sx, yax, BLUE, 1, "2 3")
text(x0 + 12, PY + PH - 14, "S: stem age of the genus", 10.5, fill="#444")

# ---------------------------------------------------------------- (iii) genus subtree
x0 = PX[2]
panel(x0, PY, PW, PH, "(iii) Dated subtree per genus")
text(x0 + 12, PY + 40, "barcode tree, rooted on related genera", 10.5, fill="#444")
xt = x0 + 182
S, MIN, MAX = 30, 10, 27
C = (MIN * MAX) ** 0.5  # no TimeTree age for this genus: the midpoint of the range (#42)
xg = lambda a: xt - a * 5.0  # noqa: E731  30 Ma -> 150 px
gy = [PY + 108 + i * 22 for i in range(4)]
ry = [PY + 70, PY + 88]
crown = (xg(C), (gy[0] + gy[3]) / 2)
inner1 = (xg(0.6 * C), (gy[0] + gy[1]) / 2)  # nodes below the crown at fixed fractions of it
inner2 = (xg(0.4 * C), (gy[2] + gy[3]) / 2)
elbow(*crown, *inner1, ORANGE, 2.2)
elbow(*crown, *inner2, ORANGE, 2.2)
for i in (0, 1):
    elbow(*inner1, xt, gy[i], ORANGE, 2.2)
for i in (2, 3):
    elbow(*inner2, xt, gy[i], ORANGE, 2.2)
rel = (xg(18), (ry[0] + ry[1]) / 2)
stem = (xg(S), (rel[1] + crown[1]) / 2)
elbow(*stem, *rel, GREY, 1.8, "3 3")
for y in ry:
    elbow(*rel, xt, y, GREY, 1.8, "3 3")
elbow(*stem, *crown, ORANGE, 2.2)
dot(*stem, 4.5, "#fff", INK)
text(stem[0] - 6, stem[1] - 8, "root", 10, "end", fill="#444")
dot(*crown, 4, ORANGE)
text(xt + 6, ry[1] + 4, "relatives", 10, fill="#777")
text(xt + 6, gy[3] + 4, "BINs", 10, fill=ORANGE)
# crown age range
by = gy[3] + 32
line(xg(MAX), by, xg(MIN), by, "#444", 5, cap="butt")
line(xg(MAX), by - 6, xg(MAX), by + 6, "#444", 1.4)
line(xg(MIN), by - 6, xg(MIN), by + 6, "#444", 1.4)
line(xg(C), crown[1] + 6, xg(C), by - 8, ORANGE, 1, "2 3")
dot(xg(C), by, 3.5, ORANGE)
text((xg(MAX) + xg(MIN)) / 2, by - 9, "crown age range", 10, "middle", fill="#444")
text((xg(MAX) + xg(MIN)) / 2, by + 18, "C: TimeTree crown if any,", 10, "middle", fill="#444")
text((xg(MAX) + xg(MIN)) / 2, by + 31, "else midpoint of range", 10, "middle", fill="#444")
yax = PY + 262
axis(xg, [30, 20, 10, 0], yax)
line(xg(S), stem[1] + 6, xg(S), yax, BLUE, 1, "2 3")
text(xg(S) + 3, yax - 5, "S", 10.5, weight="bold", fill=BLUE)
text(x0 + 12, PY + PH - 28, "min: BIN divergence ÷ rate", 10.5, fill="#444")
text(x0 + 12, PY + PH - 13, "max: 0.9 S and TimeTree crown", 10.5, fill="#444")

# ---------------------------------------------------------------- (iv) grafting
x0 = PX[3]
panel(x0, PY, PW, PH, "(iv) Grafting")
text(x0 + 12, PY + 40, "genus tip replaced by its subtree", 10.5, fill="#444")
xt = x0 + 200
xf = lambda a: xt - a * 5.0  # noqa: E731
gy = [PY + 120 + i * 22 for i in range(4)]  # the subtree of (iii), crown at C
sis_y = PY + 78
crown = (xf(C), (gy[0] + gy[3]) / 2)
node = (xf(S), (sis_y + crown[1]) / 2)
line(node[0] - 14, node[1], node[0], node[1], BLUE, 2.2)
elbow(*node, xt, sis_y, BLUE, 2.2)
elbow(*node, *crown, BLUE, 2.2)
inner1 = (xf(0.6 * C), (gy[0] + gy[1]) / 2)
inner2 = (xf(0.4 * C), (gy[2] + gy[3]) / 2)
elbow(*crown, *inner1, ORANGE, 2.2)
elbow(*crown, *inner2, ORANGE, 2.2)
for i in (0, 1):
    elbow(*inner1, xt, gy[i], ORANGE, 2.2)
for i in (2, 3):
    elbow(*inner2, xt, gy[i], ORANGE, 2.2)
dot(*node, 4, BLUE)
dot(*crown, 4, ORANGE)
text(xt + 6, sis_y + 4, "sister", 10.5)
text(xt + 6, gy[3] + 4, "BINs", 10.5, fill=ORANGE)
text(node[0] - 6, node[1] - 8, "S", 11, "end", "bold", BLUE)
text(crown[0] + 4, crown[1] - 8, "C", 11, weight="bold", fill=ORANGE)
text((node[0] + crown[0]) / 2, crown[1] + 16, "S − C", 10.5, "middle", fill=BLUE)
yax = PY + 262
axis(xf, [30, 20, 10, 0], yax)
text(x0 + 12, PY + PH - 28, "no rescaling: both parts in Ma;", 10.5, fill="#444")
text(x0 + 12, PY + PH - 13, "genera run in parallel", 10.5, fill="#444")

# arrows between panels
for a, b in zip(PX[:-1], PX[1:]):
    arrow(a + PW + 1, PY + PH / 2, b - 3, PY + PH / 2, "#555", 1.6)

# ---------------------------------------------------------------- (v) samples
y0 = PY + PH + 22
panel(15, y0, 970, 170, "Samples: placement and diversity")
cy = y0 + 92
rect(32, cy - 26, 150, 52, fill=SOFT)
text(107, cy - 6, "metabarcoding", 11, "middle")
text(107, cy + 10, "sequence (ASV)", 11, "middle")
arrow(182, cy, 222, cy, "#555")
rect(226, cy - 26, 150, 52, fill=SOFT)
text(301, cy - 6, "best match among", 11, "middle")
text(301, cy + 10, "reference sequences", 11, "middle")
outcomes = [
    ("≥ 97%", "assigned to that BIN", ORANGE, cy - 50),
    ("80–97%", "placed on the tree by nearest BIN", ORANGE, cy),
    ("&lt; 80%", "left out and reported", "#a33", cy + 50),
]
for pct, what, col, oy in outcomes:
    out.append(
        f'<path d="M376,{cy} C400,{cy} 400,{oy} 424,{oy}" fill="none" stroke="#555" '
        f'stroke-width="1.6" marker-end="url(#arr)"/>'
    )
    rect(428, oy - 17, 250, 34, fill="#fff", stroke=col, rx=5)
    text(440, oy + 4, pct, 11, weight="bold", fill=col)
    text(500, oy + 4, what, 10.5)
for oy in (cy - 50, cy):
    out.append(
        f'<path d="M678,{oy} C705,{oy} 705,{cy} 732,{cy}" fill="none" stroke="#555" '
        f'stroke-width="1.6" marker-end="url(#arr)"/>'
    )
rect(736, cy - 46, 234, 92, fill="#fbeee2", rx=6)
text(853, cy - 24, "Diversity in Ma of history", 11.5, "middle", "bold")
text(853, cy - 4, "Faith's PD, mean pairwise", 10.5, "middle")
text(853, cy + 12, "distance; between samples:", 10.5, "middle")
text(853, cy + 28, "UniFrac, PhyloSor", 10.5, "middle")

svg = (
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" '
    f'font-family="{FONT}">\n'
    '<defs><marker id="arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" '
    'markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="#555"/>'
    "</marker></defs>\n"
    f'<rect width="{W}" height="{H}" fill="#fff"/>\n' + "\n".join(o for o in out if o) + "\n</svg>\n"
)
Path(__file__).with_name("fig2.svg").write_text(svg)
