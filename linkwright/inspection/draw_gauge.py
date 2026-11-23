"""Rebuild the 1:1 commissioning dial and check sightlines at both nominal reversals."""

import json
import math
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.units import mm
from reportlab.pdfgen import canvas

HERE = Path(__file__).resolve().parent
SPEC = json.loads((HERE.parent / "kit/design.json").read_text())
BASE = SPEC["base"]
G, A, B, C = (SPEC["spacings"][k] for k in ("ground", "crank", "coupler", "rocker"))
RADIUS = SPEC["links"]["width"] / 2
EXT = SPEC["links"]["pointer_extension"]
CRANK_RING = 55.5
ROCKER_RING = 138.0
INK, BLUE, OCHRE, FAINT = "#183347", "#204F7B", "#A44319", "#ADBBC1"


def xy(origin, radius, angle):
    theta = math.radians(angle)
    return origin[0] + radius * math.cos(theta), origin[1] + radius * math.sin(theta)


def distance_to_segment(p, a, b):
    dx, dy = b[0] - a[0], b[1] - a[1]
    fraction = min(1.0, max(0.0, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / (dx * dx + dy * dy)))
    return math.hypot(p[0] - a[0] - fraction * dx, p[1] - a[1] - fraction * dy)


def check_visible():
    """The extended reversal covers the direct crank graduation with the coupler."""
    o, d = (0, 0), (G, 0)
    for extended, radius_ob in [(True, B + A), (False, B - A)]:
        alpha = math.degrees(math.acos((G * G + radius_ob**2 - C * C) / (2 * G * radius_ob)))
        theta = alpha if extended else alpha + 180
        a, b = xy(o, A, theta), xy(o, radius_ob, alpha)
        phi = math.degrees(math.atan2(b[1], b[0] - G))
        t = xy(d, C + EXT, phi)
        links = [(o, a), (a, b), (d, t)]

        def visible(point):
            # All nominal mechanical outlines, including the enlarged grip at A.
            strip = min(distance_to_segment(point, *seg) for seg in links) > RADIUS
            hardware_radius = SPEC["joint"]["collar_od_envelope_max"] / 2
            return strip and all(math.dist(point, pin) > hardware_radius for pin in (o, a, b, d))

        def tick_clear(center, radius, angle):
            # Distance to any outline is 1-Lipschitz. Midpoint margin > half a
            # 4-mm tick proves the open segment misses every modeled body.
            for step in range(8):
                p = xy(center, radius + (step + 0.5) / 2, angle)
                if min(distance_to_segment(p, *seg) for seg in links) - RADIUS <= 0.25:
                    return False
                if (
                    min(math.dist(p, k) for k in (o, a, b, d)) - SPEC["joint"]["collar_od_envelope_max"] / 2
                    <= 0.25
                ):
                    return False
            return True

        assert 95 <= phi <= 166 and tick_clear(d, ROCKER_RING, phi)
        assert tick_clear(o, CRANK_RING, theta + 180 if extended else theta)
        if extended:
            assert not visible(xy(o, CRANK_RING, theta))
        print(
            f"Crank {theta:.2f}°, rocker {phi:.2f}°: rocker tick clear; use "
            f'{"opposite" if extended else "direct"} crank sightline.'
        )


def main():
    check_visible()
    # The physical 195 x 250-mm page fits the printable region of A4 and US Letter.
    pdf = canvas.Canvas(str(HERE / "angle-gauge.pdf"), pagesize=(195 * mm, 250 * mm), invariant=1)
    pdf.setTitle("Linkwright | Full-size crank and rocker reversal gauge")
    pdf.scale(mm, mm)

    def line(a, b, color=INK, width=0.18):
        pdf.setStrokeColor(colors.HexColor(color))
        pdf.setLineWidth(width)
        pdf.line(*a, *b)

    def text(x, y, s, size=2.8, color=INK, bold=False, center=False):
        pdf.setFont("Helvetica-Bold" if bold else "Helvetica", size)
        pdf.setFillColor(colors.HexColor(color))
        (pdf.drawCentredString if center else pdf.drawString)(x, y, s)

    def circle(origin, r, color):
        pdf.setStrokeColor(colors.HexColor(color))
        pdf.setLineWidth(0.18)
        pdf.circle(*origin, r)

    def marks(origin, r, angles, color, full):
        for angle in angles:
            if full and 123 < angle < 147:
                continue  # The output scale crosses here. Its opposite input graduations remain visible.
            length = 3.4 if angle % 10 == 0 else 2.5 if angle % 5 == 0 else 1.25
            line(xy(origin, r, angle), xy(origin, r + length, angle), color, 0.24 if angle % 5 == 0 else 0.13)
            if angle % 10 == 0:
                # Near the left edge, put the full-circle 180 label inside its ticks.
                offset = -3.7 if full and (angle >= 160 and angle <= 200) else (5.4 if full else 4.9)
                xx, yy = xy(origin, r + offset, angle)
                text(xx, yy - 0.8, str(angle), 2.45, color, True, True)
        if full:
            pdf.setStrokeColor(colors.HexColor(color))
            pdf.setLineWidth(0.18)
            for start, sweep in ((0, 123), (147, 213)):
                pdf.arc(origin[0] - r, origin[1] - r, origin[0] + r, origin[1] + r, start, sweep)
        else:
            pdf.setStrokeColor(colors.HexColor(color))
            pdf.setLineWidth(0.2)
            pdf.arc(
                origin[0] - r,
                origin[1] - r,
                origin[0] + r,
                origin[1] + r,
                angles.start,
                angles.stop - 1 - angles.start,
            )
        circle(origin, 8, FAINT)
        line((origin[0] - 2, origin[1]), (origin[0] + 2, origin[1]))
        line((origin[0], origin[1] - 2), (origin[0], origin[1] + 2))

    left, bottom = 10.5, 18
    o = (left - BASE["xmin"], bottom - BASE["ymin"])
    d = (o[0] + G, o[1])

    def at(x, y):
        return x + o[0], y + o[1]

    pdf.setLineWidth(0.2)
    pdf.setStrokeColor(colors.HexColor(FAINT))
    pdf.rect(left, bottom, BASE["xmax"] - BASE["xmin"], BASE["ymax"] - BASE["ymin"])
    text(left, 240, "Linkwright  /  reversal-angle gauge", 5.1, bold=True)
    text(
        left,
        232.8,
        "LW-159-R1    /    Top view   /    Degrees counterclockwise   /    Print at actual size",
        2.95,
    )
    line(o, (o[0] + 46, o[1]), FAINT)
    line((o[0] + 65, o[1]), d, FAINT)
    marks(o, CRANK_RING, range(360), BLUE, True)
    marks(d, ROCKER_RING, range(95, 167), OCHRE, False)
    text(o[0], o[1] - 12, "O   crank", 3.2, BLUE, True, True)
    text(d[0] - 3, d[1] - 12, "D   rocker", 3.2, OCHRE, True, True)
    text(*at(-56, 136), "Look straight down.", 3.55, bold=True)
    text(*at(-56, 129), "Extend D-B to the orange arc.", 3.1, OCHRE, True)
    text(*at(-56, 122), "Extend O-A to the blue circle.", 3.1, BLUE, True)
    text(*at(-56, 115), "For hidden ticks or the blank sector,", 2.9)
    text(*at(-56, 110), "read backwards along A-O. Add 180°", 2.9)
    text(*at(-56, 105), "and wrap back into 0-360°.", 2.9)
    # Both 50 mm scales are in the vacant lower-right corner of the base.
    x, y = at(104, -53)
    line((x, y), (x - 50, y), INK, 0.32)
    line((x, y), (x, y + 50), INK, 0.32)
    for k in (0, 50):
        line((x - k, y - 1.2), (x - k, y + 1.2), INK, 0.24)
        line((x - 1.2, y + k), (x + 1.2, y + k), INK, 0.24)
    text(x - 25, y + 2.4, "50.00 mm X", 2.7, center=True)
    text(x - 18.5, y + 25, "50.00 mm Y", 2.7)
    text(
        left,
        11,
        "Measure both 50 mm rules on this print. Include reading error in the 49.75-50.25 mm limits.",
        2.83,
    )
    text(
        left,
        6.1,
        "Trim to the board outline. Open the two 16 mm pivot keep-outs. Align both centres before use.",
        2.83,
    )
    pdf.showPage()
    pdf.save()


if __name__ == "__main__":
    main()
