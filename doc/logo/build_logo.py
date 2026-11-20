#!/usr/bin/env python3
"""Build the KARST logo files in this folder.

Writes the SVG sources (text converted to outlines, so they do not depend on
installed fonts), then renders transparent PNGs with headless Google Chrome
and downsizes them with ImageMagick.

    python3 doc/logo/build_logo.py

Fonts: Raleway Bold / Regular (SIL Open Font License 1.1). Set KARST_FONT_BOLD
and KARST_FONT_REGULAR to their .ttf files if the defaults below are missing,
and CHROME to the Chrome binary.
"""

import math
import os
import random
import shutil
import subprocess
import tempfile

from fontTools.pens.boundsPen import BoundsPen
from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.pens.transformPen import TransformPen
from fontTools.ttLib import TTFont

HERE = os.path.dirname(os.path.abspath(__file__))
_RMD_FONTS = ("/Library/Frameworks/R.framework/Versions/Current/Resources/"
              "library/rmarkdown/rmd/h/bootstrap/css/fonts")
FONT_BOLD = os.environ.get("KARST_FONT_BOLD", f"{_RMD_FONTS}/RalewayBold.ttf")
FONT_REGULAR = os.environ.get("KARST_FONT_REGULAR", f"{_RMD_FONTS}/Raleway.ttf")
CHROME = os.environ.get(
    "CHROME", "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")

# Palette -------------------------------------------------------------------
INK = "#0B2F3A"        # border, wordmark text
SKY_TOP = "#0E3946"
SKY_LOW = "#5FA597"
FAR = "#A9D3C7"        # farthest peaks (mist)
MID = "#5E9E91"
FRONT = "#14444C"
WATER_TOP = "#1A525B"
WATER_LOW = "#0B2F3A"
MIST = "#D4EAE3"       # waterline, reflection bars
LIME = "#EEF6F2"       # limestone white, text on dark
SUN = "#F2B04A"
SUN_GLOW = "#F7D58C"
TEAL_TEXT = "#3D6C70"  # subtitle on light


# Text to outlines ----------------------------------------------------------
def _num(v):
    t = f"{v:.2f}".rstrip("0").rstrip(".")
    return "0" if t == "-0" else t


class Font:
    def __init__(self, path):
        self.tt = TTFont(path)
        self.glyphs = self.tt.getGlyphSet()
        self.cmap = self.tt.getBestCmap()
        self.upm = self.tt["head"].unitsPerEm
        self.hmtx = self.tt["hmtx"]

    def draw(self, text, size, x, y, tracking, pens):
        """Draw text with its baseline at y; return the advance width."""
        scale = size / self.upm
        cx = x
        for i, ch in enumerate(text):
            name = self.cmap[ord(ch)]
            for pen in pens:
                self.glyphs[name].draw(
                    TransformPen(pen, (scale, 0, 0, -scale, cx, y)))
            cx += self.hmtx[name][0] * scale
            if i < len(text) - 1:
                cx += tracking * size
        return cx - x


BOLD = Font(FONT_BOLD)
REGULAR = Font(FONT_REGULAR)


def text_runs(runs, size, y, tracking=0.0, x=0.0, anchor="start",
              space_between_runs=True):
    """runs: list of (text, Font, fill). Returns (svg, ink_bounds).

    anchor "middle" centres the ink (not the advance box) on x.
    """
    drawn = []
    bounds = BoundsPen(None)
    cx = 0.0
    for k, (text, font, fill) in enumerate(runs):
        pen = SVGPathPen(font.glyphs, ntos=_num)
        cx += font.draw(text, size, cx, y, tracking, [pen, bounds])
        if space_between_runs and k < len(runs) - 1:
            cx += tracking * size
        drawn.append((pen.getCommands(), fill))
    x0, y0, x1, y1 = bounds.bounds
    dx = x - x0 if anchor == "start" else x - (x0 + x1) / 2
    svg = "".join(f'<path d="{d}" fill="{fill}"/>' for d, fill in drawn if d)
    svg = f'<g transform="translate({dx:.2f},0)">{svg}</g>'
    return svg, (x0 + dx, y0, x1 + dx, y1)


# Landscape -----------------------------------------------------------------
def tower(x, c, h, w, p=2.8, lean=0.0):
    """One tower-karst peak: steep sides, rounded top; lean skews it."""
    u = (x - c) / w
    u *= (1 + lean) if u > 0 else (1 - lean)
    return h * math.exp(-abs(u) ** p)


def ridge(x, towers, base, q=6.0):
    """Soft maximum of a low base and the towers, so saddles stay rounded."""
    s = base ** q + sum(tower(x, *t) ** q for t in towers)
    return s ** (1.0 / q)


# (centre, height, half-width[, shape, lean]) in px; heights above the water.
FAR_TOWERS = [(28, 150, 20), (92, 196, 22, 3.0, 0.15), (148, 166, 18),
              (206, 232, 24, 2.8, -0.1), (262, 150, 17), (304, 188, 20),
              (350, 214, 22, 3.0, 0.2), (468, 204, 24, 2.8, -0.15),
              (518, 158, 20)]
MID_TOWERS = [(0, 122, 22), (58, 150, 22, 3.0, -0.2), (118, 112, 17),
              (176, 168, 23, 2.9, 0.15), (252, 138, 20), (318, 160, 22),
              (402, 98, 22), (452, 146, 21, 3.0, 0.2), (515, 124, 19)]
FRONT_TOWERS = [(38, 96, 27), (128, 128, 29, 2.9, 0.12), (204, 70, 21),
                (272, 112, 25, 2.9, -0.1), (352, 84, 23), (438, 122, 29),
                (506, 76, 23)]
LAYER_BASE = {"far": 32.0, "mid": 22.0, "front": 12.0}


def layer_path(towers, base, water_y, hscale, x_from=-20, x_to=540, step=1.5):
    pts = [(x_from, water_y)]
    x = x_from
    while x <= x_to:
        pts.append((x, water_y - hscale * ridge(x, towers, base)))
        x += step
    pts.append((x_to, water_y))
    return "M" + " L".join(f"{px:.1f},{py:.1f}" for px, py in pts) + " Z"


def reflection_bars(water_y, hscale, seed=11):
    """Sparse vertical bars hanging from the waterline: the observed
    abundances, a broken reflection of the front ridge."""
    rng = random.Random(seed)
    bars = []
    pitch, width = 9.0, 5.2
    x = 22.0
    while x <= 498:
        keep = rng.random() < 0.58
        h = 0.62 * hscale * ridge(x, FRONT_TOWERS, LAYER_BASE["front"])
        # keep each bar clear of the hexagon's lower edges
        edge_y = 446.5 + math.tan(math.radians(30)) * min(x - 6.25, 513.75 - x)
        h = min(h, edge_y - 22 - (water_y + 5))
        if keep and h >= 9:
            bars.append(f'<rect x="{x - width / 2:.1f}" y="{water_y + 5:.1f}" '
                        f'width="{width}" height="{h:.1f}" rx="1.4"/>')
        x += pitch
    return "".join(bars)


def hexagon(cx=260.0, cy=300.0, r=293.0):
    pts = [(cx + r * math.cos(math.radians(a)),
            cy + r * math.sin(math.radians(a)))
           for a in (-90, -30, 30, 90, 150, 210)]
    return " ".join(f"{x:.2f},{y:.2f}" for x, y in pts)


def scene(uid, water_y, hscale, sun_r=25):
    """Sky, sun, three peak layers, water and reflection, in 520 x 600."""
    sx, sy, sr = 404, water_y - 128 * hscale, sun_r * hscale
    return f"""
  <defs>
    <linearGradient id="{uid}-sky" x1="0" y1="0" x2="0" y2="{water_y}" gradientUnits="userSpaceOnUse">
      <stop offset="0" stop-color="{SKY_TOP}"/>
      <stop offset="0.55" stop-color="#1F5E66"/>
      <stop offset="1" stop-color="{SKY_LOW}"/>
    </linearGradient>
    <radialGradient id="{uid}-glow" cx="{sx}" cy="{sy}" r="{sr * 3.2}" gradientUnits="userSpaceOnUse">
      <stop offset="0" stop-color="{SUN_GLOW}" stop-opacity="0.55"/>
      <stop offset="1" stop-color="{SUN_GLOW}" stop-opacity="0"/>
    </radialGradient>
    <linearGradient id="{uid}-water" x1="0" y1="{water_y}" x2="0" y2="600" gradientUnits="userSpaceOnUse">
      <stop offset="0" stop-color="{WATER_TOP}"/>
      <stop offset="1" stop-color="{WATER_LOW}"/>
    </linearGradient>
    <linearGradient id="{uid}-bars" x1="0" y1="{water_y}" x2="0" y2="{water_y + 0.62 * hscale * 130}" gradientUnits="userSpaceOnUse">
      <stop offset="0" stop-color="{MIST}" stop-opacity="0.95"/>
      <stop offset="1" stop-color="{MIST}" stop-opacity="0.08"/>
    </linearGradient>
  </defs>
  <rect x="0" y="0" width="520" height="{water_y}" fill="url(#{uid}-sky)"/>
  <circle cx="{sx}" cy="{sy}" r="{sr * 3.2}" fill="url(#{uid}-glow)"/>
  <circle cx="{sx}" cy="{sy}" r="{sr}" fill="{SUN}"/>
  <path d="{layer_path(FAR_TOWERS, LAYER_BASE['far'], water_y, hscale)}" fill="{FAR}"/>
  <path d="{layer_path(MID_TOWERS, LAYER_BASE['mid'], water_y, hscale)}" fill="{MID}"/>
  <path d="{layer_path(FRONT_TOWERS, LAYER_BASE['front'], water_y, hscale)}" fill="{FRONT}"/>
  <rect x="0" y="{water_y}" width="520" height="{600 - water_y}" fill="url(#{uid}-water)"/>
  <g fill="url(#{uid}-bars)">{reflection_bars(water_y, hscale)}</g>
  <rect x="0" y="{water_y - 1.2}" width="520" height="2.4" fill="{MIST}"/>"""


def hex_badge(uid, with_text=True):
    hx = hexagon()
    if with_text:
        body = scene(uid, water_y=436, hscale=1.0)
        title, _ = text_runs([("KARST", BOLD, LIME)], size=72, y=152,
                             tracking=0.12, x=260, anchor="middle")
    else:
        body = scene(uid, water_y=404, hscale=1.2)
        title = ""
    return f"""<g>
  <clipPath id="{uid}-clip"><polygon points="{hx}"/></clipPath>
  <g clip-path="url(#{uid}-clip)">{body}
  </g>
  {title}
  <polygon points="{hx}" fill="none" stroke="{INK}" stroke-width="14" stroke-linejoin="round"/>
</g>"""


def svg_doc(w, h, inner, title, desc):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
            f'viewBox="0 0 {w} {h}" role="img" aria-labelledby="t d">\n'
            f'<title id="t">{title}</title>\n<desc id="d">{desc}</desc>\n'
            f'{inner}\n</svg>\n')


DESC = ("Tower-karst peaks over water at dusk. The peaks are the fitness "
        "landscape; the sparse bars below the waterline, its broken "
        "reflection, are the observed karyotype abundances it is "
        "reconstructed from.")


def wordmark(dark):
    """Hex icon, KARST, and the name spelled out on two lines."""
    ink = LIME if dark else INK
    sub = "#A9CFC6" if dark else TEAL_TEXT
    accent = SUN if dark else "#C47A12"
    icon_scale = 300 / 600
    icon = (f'<g transform="translate(20,20) scale({icon_scale})">'
            f'{hex_badge("wm" + ("d" if dark else "l"), with_text=False)}</g>')
    x0 = 20 + 520 * icon_scale + 40
    title, tb = text_runs([("KARST", BOLD, ink)], size=150, y=178,
                          tracking=0.10, x=x0)
    lines = [
        [("K", BOLD, accent), ("aryotype ", REGULAR, sub),
         ("A", BOLD, accent), ("bundance-based ", REGULAR, sub),
         ("R", BOLD, accent), ("econstruction", REGULAR, sub)],
        [("of fitness landscape from ", REGULAR, sub),
         ("S", BOLD, accent), ("parse ", REGULAR, sub),
         ("T", BOLD, accent), ("ime points", REGULAR, sub)],
    ]
    subs = []
    right = tb[2]
    for k, runs in enumerate(lines):
        s, b = text_runs(runs, size=27, y=236 + k * 38, x=x0 + 4,
                         space_between_runs=False)
        subs.append(s)
        right = max(right, b[2])
    w = int(math.ceil(right + 30))
    inner = icon + title + "".join(subs)
    return svg_doc(w, 340, inner, "KARST", DESC)


def main():
    files = {
        "karst_hex.svg": svg_doc(520, 600, hex_badge("hx"), "KARST", DESC),
        "karst_icon.svg": svg_doc(520, 600, hex_badge("ic", with_text=False),
                                  "KARST icon", DESC),
        "karst_wordmark.svg": wordmark(dark=False),
        "karst_wordmark_dark.svg": wordmark(dark=True),
    }
    for name, svg in files.items():
        with open(os.path.join(HERE, name), "w") as fh:
            fh.write(svg)
        print("wrote", name)

    with tempfile.TemporaryDirectory() as tmp:
        hex2x = os.path.join(tmp, "hex2x.png")
        render("karst_hex.svg", hex2x, tmp)
        shutil.copy(hex2x, os.path.join(HERE, "karst_hex.png"))
        downsize(hex2x, "karst_hex_small.png", "240x")
        # the package logo (README, pkgdown)
        figures = os.path.join(HERE, "..", "..", "man", "figures")
        os.makedirs(figures, exist_ok=True)
        shutil.copy(os.path.join(HERE, "karst_hex_small.png"),
                    os.path.join(figures, "logo.png"))
        icon2x = os.path.join(tmp, "icon2x.png")
        render("karst_icon.svg", icon2x, tmp)
        shutil.copy(icon2x, os.path.join(HERE, "karst_icon.png"))
        for size in (512, 180, 64, 32):
            downsize(icon2x, f"karst_icon_{size}.png", f"x{size}",
                     square=size)
        for name in ("karst_wordmark", "karst_wordmark_dark"):
            render(f"{name}.svg", os.path.join(HERE, f"{name}.png"), tmp)
    for name in sorted(os.listdir(HERE)):
        if name.endswith(".png"):
            print("wrote", name)


# PNG rendering -------------------------------------------------------------
def svg_size(name):
    with open(os.path.join(HERE, name)) as fh:
        head = fh.read(400)
    w = float(head.split('width="', 1)[1].split('"', 1)[0])
    h = float(head.split('height="', 1)[1].split('"', 1)[0])
    return int(round(w)), int(round(h))


def render(name, png, tmp, scale=2):
    """Screenshot the SVG in headless Chrome on a transparent page."""
    w, h = svg_size(name)
    page = os.path.join(tmp, "page.html")
    src = "file://" + os.path.join(HERE, name)
    with open(page, "w") as fh:
        fh.write("<!doctype html><html><head><style>html,body{margin:0;"
                 "padding:0;background:transparent}img{display:block;"
                 f"width:{w}px;height:{h}px}}</style></head><body>"
                 f'<img src="{src}"></body></html>')
    subprocess.run(
        [CHROME, "--headless=new", "--disable-gpu", "--hide-scrollbars",
         "--allow-file-access-from-files",
         f"--user-data-dir={os.path.join(tmp, 'chrome')}",
         "--default-background-color=00000000",
         f"--force-device-scale-factor={scale}", f"--window-size={w},{h}",
         f"--screenshot={png}", "file://" + page],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        timeout=120)


def downsize(src, out, geometry, square=None):
    cmd = ["magick", src, "-filter", "Lanczos", "-resize", geometry]
    if square:
        cmd += ["-background", "none", "-gravity", "center",
                "-extent", f"{square}x{square}"]
    subprocess.run(cmd + ["-strip", os.path.join(HERE, out)], check=True)


if __name__ == "__main__":
    main()
