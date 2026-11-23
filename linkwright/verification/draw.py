"""Draw the nominal motion of the inspected kit from its dimensions and receipt."""

import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = json.loads((HERE.parent / "kit" / "design.json").read_text())
receipt = json.loads((HERE / "receipt.json").read_text())["motion"]
g, a, b, c = [spec["spacings"][s] for s in ("ground", "crank", "coupler", "rocker")]
black, gray, orange, teal = "#18343c", "#50656d", "#a04819", "#086a67"
svg = [
    '<svg xmlns="http://www.w3.org/2000/svg" width="700" height="380" viewBox="0 0 700 380" role="img" aria-labelledby="title desc">',
    '<title id="title">One crank revolution, two rocker strokes</title>',
    '<desc id="desc">The rocker reverses at crank angles 56.612 and 288.078 degrees. The long stroke spans 231.466 degrees. The return takes 128.534 degrees.</desc>',
    '<rect width="700" height="380" rx="8" fill="#fafbf7"/><g font-family="DejaVu Sans,Arial,sans-serif">',
]


def label(x, y, text, size=14, color=black, weight=400, anchor="start"):
    svg.append(
        f'<text x="{x:.2f}" y="{y:.2f}" text-anchor="{anchor}" fill="{color}" font-size="{size}" font-weight="{weight}">{text}</text>'
    )


def line(x, y, x2, y2, color, width, dash=""):
    extra = f' stroke-dasharray="{dash}"' if dash else ""
    svg.append(
        f'<path d="M {x:.3f} {y:.3f} L {x2:.3f} {y2:.3f}" fill="none" stroke="{color}" stroke-width="{width}"{extra}/>'
    )


def x(t):
    return 66 + 585 * t / 360


def y(phi):
    return 258 - 1.9 * (phi - 100)


def rocker(t):
    t = math.radians(t)
    ax, ay = a * math.cos(t), a * math.sin(t)
    ray = math.atan2(ay, ax - g)
    d = math.hypot(ax - g, ay)
    phi = ray - math.acos((c * c + d * d - b * b) / (2 * c * d))
    return math.degrees(phi) % 360


label(29, 35, "One crank revolution, two rocker strokes", 24, black, 700)
label(
    29, 59, "Rounded LW-159-R1  •  counter-clockwise input  •  upper assembly", 15, gray
)
label(66, 99, "Rocker angle from the ground line", 15, gray, 600)
for t in (100, 120, 140, 160):
    line(66, y(t), 651, y(t), "#d1dad4", 1)
    label(55, y(t) + 5, f"{t}°", 13, anchor="end")
for t in (0, 90, 180, 270, 360):
    line(x(t), 268, x(t), 273, gray, 1)
    label(x(t), 292, f"{t}°", 14, anchor="middle")
line(66, 268, 651, 268, gray, 1.2)
lo, hi = receipt["nominal_crank_reversal_deg_ccw"]
for left, right, color in ((0, lo, teal), (lo, hi, orange), (hi, 360, teal)):
    path = []
    for i in range(161):
        coord = left + (right - left) * i / 160
        path.append(f'{"M" if i==0 else "L"}{x(coord):.2f} {y(rocker(coord)):.2f}')
    svg.append(
        '<path d="'
        + " ".join(path)
        + f'" fill="none" stroke="{color}" stroke-width="3.5" stroke-linecap="round"/>'
    )
for t in (lo, hi):
    px, py = x(t), y(rocker(t))
    line(px, py, px, 268, gray, 1, "3 4")
    svg.append(f'<circle cx="{px:.2f}" cy="{py:.2f}" r="4.5" fill="{black}"/>')
line(x(lo), y(rocker(lo)) - 26, x(lo), y(rocker(lo)) - 7, gray, 1)
label(x(lo), y(rocker(lo)) - 32, f"{lo:.2f}° crank", 14, black, 600, "middle")
label(x(hi), y(rocker(hi)) - 14, f"{hi:.2f}° crank", 14, black, 600, "middle")
label(650, 318, "Crank angle", 14, gray, 600, "end")
line(45, 349, 66, 349, orange, 4)
label(73, 354, "Long stroke   231.47°", 16, orange, 600)
line(353, 349, 374, 349, teal, 4)
label(382, 354, "Short stroke   128.53°", 16, teal, 600)
svg.append("</g></svg>")
(HERE / "motion.svg").write_text("\n".join(svg) + "\n", encoding="utf-8")
