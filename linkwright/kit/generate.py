"""Regenerate the Linkwright reference drawings. Needs reportlab and ezdxf."""

from pathlib import Path
import json
from math import sqrt
from xml.sax.saxutils import escape

from reportlab.lib.pagesizes import A4, letter
from reportlab.lib.units import mm
from reportlab.lib.colors import HexColor
from reportlab.pdfgen.canvas import Canvas
import ezdxf
from ezdxf import units

# Keep release checksums stable when regenerating manufacturing files.
ezdxf.options.write_fixed_meta_data_for_testing = True

ROOT = Path(__file__).resolve().parent
SPEC = json.loads((ROOT / "design.json").read_text())
S = SPEC["spacings"]
LINK = SPEC["links"]
BASE = SPEC["base"]
J = SPEC["joint"]
INK, MUTED, BORDER = "#192E38", "#425963", "#AABAC0"
TEAL, BLUE, ORANGE = "#00776B", "#29536F", "#AE461D"
WHITE, TINT, SKIN = "#FFFFFF", "#F4F7F5", "#D3E9E4"
REV = SPEC["revision"]
NAMES = [
    ("Crank", "crank", "O", "A"),
    ("Coupler", "coupler", "A", "B"),
    ("Rocker + pointer", "rocker", "D", "B"),
]


class Drawing:
    """Native vectors in physical millimetres; origin at top left."""

    def __init__(self, width, height, pdf=None, title="", desc=""):
        self.w, self.h, self.pdf = width, height, pdf
        self.svg = [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}mm" height="{height}mm" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">',
            f'<title id="title">{escape(title)}</title><desc id="desc">{escape(desc)}</desc>',
        ]

    def save_svg(self, name):
        (ROOT / name).write_text("\n".join(self.svg + ["</svg>"]) + "\n", encoding="utf-8")

    def style(self, stroke=None, fill=None, width=0.25):
        if self.pdf:
            self.pdf.setLineWidth(width * mm)
            if stroke:
                self.pdf.setStrokeColor(HexColor(stroke))
            if fill:
                self.pdf.setFillColor(HexColor(fill))

    def line(self, x, y, xx, yy, color=INK, width=0.3):
        self.svg.append(
            f'<line x1="{x:.5f}" y1="{y:.5f}" x2="{xx:.5f}" y2="{yy:.5f}" stroke="{color}" stroke-width="{width}"/>'
        )
        if self.pdf:
            self.style(color, width=width)
            self.pdf.line(x * mm, (self.h - y) * mm, xx * mm, (self.h - yy) * mm)

    def rect(self, x, y, w, h, radius=0, stroke=INK, fill=None, width=0.3):
        self.svg.append(
            f'<rect x="{x:.5f}" y="{y:.5f}" width="{w:.5f}" height="{h:.5f}" rx="{radius:.5f}" stroke="{stroke or "none"}" fill="{fill or "none"}" stroke-width="{width}"/>'
        )
        if self.pdf:
            self.style(stroke, fill, width)
            self.pdf.roundRect(
                x * mm,
                (self.h - y - h) * mm,
                w * mm,
                h * mm,
                radius * mm,
                stroke=int(bool(stroke)),
                fill=int(bool(fill)),
            )

    def circ(self, x, y, r, stroke=INK, fill=None, width=0.3):
        self.svg.append(
            f'<circle cx="{x:.5f}" cy="{y:.5f}" r="{r:.5f}" stroke="{stroke or "none"}" fill="{fill or "none"}" stroke-width="{width}"/>'
        )
        if self.pdf:
            self.style(stroke, fill, width)
            self.pdf.circle(x * mm, (self.h - y) * mm, r * mm, stroke=int(bool(stroke)), fill=int(bool(fill)))

    def text(self, x, y, s, size=3.6, color=INK, bold=False, align="start"):
        self.svg.append(
            f'<text x="{x:.4f}" y="{y:.4f}" font-family="Arial,Helvetica,sans-serif" text-anchor="{align}" fill="{color}" font-size="{size}" font-weight="{700 if bold else 400}">{escape(s)}</text>'
        )
        if self.pdf:
            self.style(fill=color)
            self.pdf.setFont("Helvetica-Bold" if bold else "Helvetica", size * mm)
            draw = {
                "start": self.pdf.drawString,
                "middle": self.pdf.drawCentredString,
                "end": self.pdf.drawRightString,
            }[align]
            draw(x * mm, (self.h - y) * mm, s)

    def words(self, x, y, lines, size=3.8, lead=5.3, color=INK):
        for i, s in enumerate(lines):
            self.text(x, y + i * lead, s, size, color)

    def bar(self, a, b, width, color):
        self.line(*a, *b, color, width)
        for p in (a, b):
            self.circ(*p, width / 2, None, color)

    def cross(self, x, y, radius=2.3):
        self.line(x - radius, y, x + radius, y, MUTED, 0.15)
        self.line(x, y - radius, x, y + radius, MUTED, 0.15)


def plan(d, x, y, scale):
    """An accurately scaled positive-branch assembly at theta=0 degrees."""
    g, a, b, c = [S[k] for k in ["ground", "crank", "coupler", "rocker"]]
    bx = a + (b * b + (g - a) ** 2 - c * c) / (2 * (g - a))
    by = sqrt(b * b - (bx - a) ** 2)

    def p(xx, yy):
        return x + scale * xx, y - scale * yy

    o, aa, bb, dd = p(0, 0), p(a, 0), p(bx, by), p(g, 0)
    extra = 1 + LINK["pointer_extension"] / c
    tip = p(g + (bx - g) * extra, by * extra)
    w = LINK["width"] * scale
    d.rect(
        *p(BASE["xmin"], BASE["ymax"]),
        scale * (BASE["xmax"] - BASE["xmin"]),
        scale * (BASE["ymax"] - BASE["ymin"]),
        0,
        BORDER,
        TINT,
        0.38,
    )
    d.line(*o, *dd, BORDER, 0.4)
    d.bar(o, aa, w, BLUE)
    d.bar(dd, tip, w, TEAL)
    d.bar(aa, bb, w, ORANGE)
    for point, lab, dx, dy in [(o, "O", -7, 8), (aa, "A", 5, 7), (bb, "B", 5, -1), (dd, "D", 3, 9)]:
        d.circ(*point, 3.5 * scale, INK, WHITE, 0.35)
        d.circ(*point, 1.5 * scale, INK, None, 0.25)
        d.text(point[0] + dx, point[1] + dy, lab, 4.2, INK, True)
    d.circ(*tip, 1.7 * scale, TEAL, WHITE, 0.3)
    return tip


def assembly_svg():
    d = Drawing(
        212,
        160,
        title="Linkwright 159 assembly",
        desc="Plan of the four-bar on its base at zero crank angle. Orange is the raised coupler. Dimensioned cut sheets and side views are in the build packet.",
    )
    d.rect(0, 0, 212, 160, 0, None, WHITE)
    d.text(7, 11, "Linkwright", 8, INK, True)
    d.text(7, 18, "REFERENCE ASSEMBLY  /  LW-159-R1", 3.5, MUTED)
    plan(d, 48, 111, 0.57)
    d.text(119, 39, "Four joints, two levels", 5.6, INK, True)
    yy = 52
    for k, a, b in [
        (BLUE, "O to A", "46.6 mm input crank"),
        (ORANGE, "A to B", "84.1 mm raised coupler"),
        (TEAL, "D to B", "111.3 mm rocker"),
    ]:
        d.line(119, yy, 128, yy, k, 1.7)
        d.text(132, yy + 1, a, 4, INK, True)
        d.text(132, yy + 7, b, 3.65, MUTED)
        yy += 19
    d.words(
        119,
        113,
        ["Ground O to D, 93.8 mm.", "Pointer extends 18 mm past B.", "Turn the grip on A by hand."],
        3.7,
        6,
    )
    d.text(7, 155, "TOP VIEW   /   +y points up   /   θ = 0°   /   drawing reduced", 3.4, MUTED)
    d.save_svg("assembly.svg")


def header(d, num, title, subtitle):
    d.text(15, 14, "LINKWRIGHT   /   MECHANICAL BUILD   /   " + REV, 3.0, TEAL, True)
    d.text(15, 27, title, 7, INK, True)
    d.text(15, 35, subtitle, 3.5, MUTED)
    d.line(15, 40, d.w - 15, 40, BORDER, 0.27)
    d.line(15, d.h - 15, d.w - 15, d.h - 15, BORDER, 0.25)
    d.text(15, d.h - 10, "UNLOADED HAND-OPERATED REFERENCE   /   DIMENSIONS IN mm", 2.85, MUTED, True)
    d.text(d.w - 15, d.h - 10, f"{num} / 6", 3, MUTED, True, "end")


def note(d, x, y, title, body):
    d.text(x, y, title, 4.3, TEAL, True)
    d.words(x, y + 6, body, 3.65, 5.1)


def introduction(d):
    header(
        d, 1, "Build the 60-degree wing", "One slow hand turn drives a bare rocker pointer back and forth."
    )
    plan(d, 51, 134, 0.58)
    d.text(121, 55, "At the alignment mark", 4.75, INK, True)
    for yy, color, label, detail in [
        (67, BLUE, "Crank O-A", "Push the grip at A."),
        (88, ORANGE, "Coupler A-B", "Runs above both links."),
        (109, TEAL, "Rocker D-B", "Carries the integral pointer."),
    ]:
        d.line(121, yy, 128, yy, color, 1.3)
        d.text(131, yy + 1, label, 4.15, INK, True)
        d.text(121, yy + 8, detail, 3.6)
    d.words(121, 135, ["O and D are fixed.", "At crank zero, place B", "above the line O-D."], 3.6, 5.2)
    note(
        d,
        15,
        182,
        "What to build",
        [
            "The two holes in the base are 93.8 apart. The three moving",
            "pin spacings are 46.6, 84.1 and 111.3. Use the 18 mm pointer",
            "extension as the visible wing. Keep the mechanism horizontal.",
        ],
    )
    note(
        d,
        15,
        213,
        "Which pages reach the bench?",
        [
            "2  Parts and measured fits     3  Full-size base drill sheet",
            "4  Full-size moving links       5  Joint stacks and fitting",
            "6  Assembly, first motion and acceptance",
        ],
    )
    d.text(
        15,
        246,
        "Print at actual size. Check both 50 mm axes on pages 3, 4 and the separate angle gauge.",
        3.65,
        INK,
        True,
    )


def stock(d):
    header(
        d,
        2,
        "Stock and hardware",
        "Select parts by the finished dimensions below. Inspect them before assembly.",
    )
    rows = [
        ("1", "Flat birch plywood base", "174 × 207 × 6.0 ±0.2; board flatness 0.30 max."),
        ("3", "Birch plywood links", "3.0 ±0.1 thick; profile width 14 ±0.4; flatness 0.20 max."),
        ("4", "Smooth straight steel shafts", "Diameter 2.98-3.00. Cut from one 150 mm or longer rod."),
        ("10", "Set-screw shaft collars", "3 mm bore; 5.0 ±0.1 thick. DIN 705 A, nominal OD 7."),
        ("12", "Thrust washers", "3.1-3.3 bore, OD 6.5-8; nominal OD 7. Thickness 0.50 ±0.05."),
        ("2 + 2", "Tube spacers", "Two 3.0 and two 10.0 long, ±0.1. Bore 3.05-3.30, OD 5-8."),
        ("1", "Finger grip and adhesive", "10 ±0.3 diameter, 15 ±0.5 long; close-grain wood or polymer."),
        ("4", "Adhesive feet", "12 ±0.5 installed height, width 15 each; one under each corner."),
    ]
    d.text(15, 49, "QTY", 3.1, MUTED, True)
    d.text(35, 49, "COMPONENT / ACCEPTANCE", 3.1, MUTED, True)
    for i, (qty, label, detail) in enumerate(rows):
        y = 59 + i * 16.7
        d.text(15, y, qty, 4.3, TEAL, True)
        d.text(35, y, label, 4.08, INK, True)
        d.text(35, y + 5.5, detail, 3.45)
        d.line(15, y + 9, d.w - 15, y + 9, BORDER, 0.17)
    note(
        d,
        15,
        199,
        "Before drilling",
        [
            "Use a square drill guide or drill press, backing block, fine blade, pin gauges,",
            "calipers, scale and feeler gauges. Use the collar's matching hex key.",
            "Check a coupon first. Moving holes accept 3.10-3.20; O and D accept 3.00-3.05.",
        ],
    )
    note(
        d,
        15,
        231,
        "Shafts and collars",
        [
            "Rough-cut blanks at 32, 32, 32 and 43, then trim to the measured stacks.",
            "Indicative finished lengths O/D 29.0, B 28.9, A 39.3. Remove burrs.",
            "Enclose all set screws within the 10.2 diameter hardware envelope.",
        ],
    )


def scale_l(d, x, y):
    d.line(x, y, x + 50, y, INK, 0.34)
    d.line(x, y, x, y - 50, INK, 0.34)
    for n in range(0, 51, 10):
        size = 1.4 if n in (0, 50) else 0.85
        d.line(x + n, y - size, x + n, y + size, INK, 0.25)
        d.line(x - size, y - n, x + size, y - n, INK, 0.25)
    d.text(x + 25, y + 5, "50.00 horizontal", 3.1, MUTED, False, "middle")
    d.text(x + 2, y - 44, "50.00", 3.1, MUTED)
    d.text(x + 2, y - 40, "vertical", 3.1, MUTED)


def links_layout():
    # Coordinates in the physical production layout, x-right/y-up from lower left.
    return [(name, key, first, last, 15, y) for (name, key, first, last), y in zip(NAMES, (85, 52, 19))]


def radius_ground():
    return sum(J["ground_bore"]) / 4


def radius_moving():
    return sum(J["moving_bore"]) / 4


def blank_width():
    return BASE["xmax"] - BASE["xmin"]


def blank_height():
    return BASE["ymax"] - BASE["ymin"]


def base_holes():
    return [(-BASE["xmin"], -BASE["ymin"]), (S["ground"] - BASE["xmin"], -BASE["ymin"])]


def base_page(d):
    header(
        d,
        3,
        "Base and drill centres  /  1:1",
        "Print at 100%. Each reference bar must measure 50.00 ±0.25 before transfer.",
    )
    x = (d.w - blank_width()) / 2
    y = 44
    w, h = blank_width(), blank_height()
    d.rect(x, y, w, h, 0, INK, None, 0.25)
    d.text(x + 21, y + 9, "BASE   /   6 mm stock   /   top face", 4.0, INK, True)
    d.text(x + 21, y + 16, "Outside 174 × 207. Hold each edge within ±0.5.", 3.3)
    for xx in [12, w - 12]:
        for yy in [12, h - 12]:
            # Foot references are on the underside, indicated with a small cross only.
            d.cross(x + xx, y + h - yy, 1.3)
    scale_l(d, x + 15, y + 82)
    d.text(x + 99, y + 32, "+y", 3.4, MUTED)
    d.line(x + 101, y + 46, x + 101, y + 35, MUTED, 0.3)
    d.text(x + 120, y + 47, "+x", 3.4, MUTED)
    d.line(x + 101, y + 46, x + 118, y + 46, MUTED, 0.3)
    for (xx, yy), label in zip(base_holes(), ("O", "D")):
        d.circ(x + xx, y + h - yy, radius_ground(), INK, None, 0.2)
        d.cross(x + xx, y + h - yy)
        d.text(x + xx, y + h - yy - 5, label, 4.5, INK, True, "middle")
    ox, oy = base_holes()[0]
    dx, _ = base_holes()[1]
    q = y + h - oy + 8
    d.line(x + ox, q, x + dx, q, MUTED, 0.21)
    for xx in [ox, dx]:
        d.line(x + xx, q - 1.5, x + xx, q + 1.5, MUTED, 0.21)
    d.text(x + (ox + dx) / 2, q + 5, "93.80 ±0.25", 3.55, TEAL, True, "middle")
    d.text(x + 7, y + h - 36, "O (62.00, 62.00)      D (155.80, 62.00)", 3.65, INK, True)
    d.text(x + 7, y + h - 30, "Measure x from left, y from bottom. Accept each coordinate within ±0.25.", 3.2)
    d.text(
        x + 7,
        y + h - 24,
        "Ream the two circles to 3.00-3.05. Crosses near corners mark underside feet only.",
        3.2,
    )
    d.text(
        x + 7,
        y + h - 18,
        "Feet attach with centre 12 ±1 from each adjacent edge. Keep the underside flat.",
        3.2,
    )


def links_page(d):
    header(
        d,
        4,
        "Moving profiles  /  1:1",
        "Cut on the centre of each outline and transfer the six pin centres before removing paper.",
    )
    x = (d.w - 174) / 2 + 1
    layout_top = 49
    for name, key, first, last, lx, ly in links_layout():
        xx = x + lx
        yy = layout_top + 104 - ly
        spacing = S[key]
        ext = LINK["pointer_extension"] if key == "rocker" else 0
        width = LINK["width"]
        rad = width / 2
        d.rect(xx - rad, yy - rad, spacing + ext + width, width, rad, INK, None, 0.25)
        for px, letter in [(xx, first), (xx + spacing, last)]:
            d.circ(px, yy, radius_moving(), INK, None, 0.21)
            d.cross(px, yy)
            d.text(px, yy + 12, letter, 3.3, MUTED, True, "middle")
        d.text(xx + spacing / 2, yy - 1, name, 3.75, INK, True, "middle")
        d.text(xx + spacing / 2, yy + 4, f"{spacing:.1f} ±0.25", 3.2, MUTED, False, "middle")
        if ext:
            d.line(xx + spacing + ext, yy - 1, xx + spacing + ext, yy + 1, TEAL, 0.25)
            d.text(xx + spacing + ext, yy + 12, "T", 3.3, TEAL, True, "middle")
    scale_l(d, 28, 213)
    note(
        d,
        95,
        169,
        "Finish the moving pieces",
        [
            "All six holes finish at 3.10-3.20.",
            "Strip width 14.0 ±0.4; ends R7.",
            "Keep profile within 0.2 of template.",
            "Grain of face plies runs along the bar.",
        ],
    )
    note(
        d,
        95,
        198,
        "Read the pointer",
        [
            "T is 18.0 ±0.3 beyond B on line D-B.",
            "Mark T in ink. Its end cap extends",
            "another 7 mm. T is a mark, not a hole.",
        ],
    )
    d.words(
        17,
        236,
        [
            "Hand work: drill/ream with scrap backing. Measure all bore spacings after finishing.",
            "CNC: use the separate layered DXFs or plain production SVGs; apply your own kerf offset.",
        ],
        3.55,
        5.5,
    )


def stacks(d):
    header(
        d,
        5,
        "How the four pins fit",
        "Side views enlarged. Both drawings use the same vertical scale and share the base-top datum.",
    )
    scale = 2.1
    base_y = 163

    def side(x, z, h, width, fill):
        d.rect(x - width / 2, base_y - (z + h) * scale, width, h * scale, 0, INK, fill, 0.19)

    def call(x, z, label, dy=0):
        y = base_y - z * scale
        d.line(x + 10, y, x + 14, y + dy, MUTED, 0.18)
        d.text(x + 15, y + dy + 1, label, 3.03, INK)

    for x in (42, 123):
        side(x, -6, 6, 38, "#E4DBCA")
        d.line(x - 18, base_y, x + 12, base_y, MUTED, 0.18)
    side(42, -11.3, 28.9, 4, "#AAB8BD")
    for z in [-11, 0, 12.25]:
        side(42, z, 5, 12, "#B7C3C7")
    side(42, 5, 3, 10, "#E2E9EA")
    for z in (8, 11.5):
        side(42, z, 0.5, 11, "#607682")
    side(42, 8.5, 3, 31, "#7DABBF")
    for z, label, dy in [
        (-8.6, "bottom collar", 0),
        (-3.0, "6 base", 0),
        (2.5, "clamp collar", 0),
        (6.5, "3 spacer", 1),
        (10, "lower link", -1),
        (15, "top collar", -1),
    ]:
        call(42, z, label, dy)
    side(123, 2.7, 39.3, 4, "#AAB8BD")
    for z in (3, 26.25):
        side(123, z, 5, 12, "#B7C3C7")
    for z in (8, 11.5, 22, 25.5):
        side(123, z, 0.5, 11, "#607682")
    side(123, 8.5, 3, 30, "#7DABBF")
    side(123, 12, 10, 10, "#E2E9EA")
    side(123, 22.5, 3, 30, "#E5A980")
    side(123, 32, 15, 14, "#D1B67A")
    d.line(123, base_y - 32 * scale, 123, base_y - 42 * scale, MUTED, 0.3)
    for z, label, dy in [
        (39, "A grip only", 2),
        (28.75, "top collar", 0),
        (24, "coupler", 0),
        (17, "10 spacer", 0),
        (10, "lower link", 0),
        (5.5, "bottom collar", 0),
        (-3, "6 base", 0),
    ]:
        call(123, z, label, dy)
    d.text(22, 52, "O + D  /  fixed axes", 4.3, TEAL, True)
    d.text(103, 52, "A + B  /  moving pins", 4.3, TEAL, True)
    note(
        d,
        15,
        197,
        "Stack order, from the underside upward",
        [
            "O and D: collar, base, clamp collar, 3 spacer, washer, lower link, washer, free play, top collar.",
            "A and B: collar, washer, lower link, washer, 10 spacer, washer, coupler, washer, free play, top collar.",
        ],
    )
    note(
        d,
        15,
        224,
        "Set each moving pack with a feeler gauge",
        [
            "Use 0.15-0.35 axial free play; .25 is the setting target. O and D clamp the board separately.",
            "Trim free shaft ends flush to +0.6. At A leave 0.5-1.0 gap above collar and bond the grip",
            "to 9.5-10.5 of shaft. Keep its adhesive off collars, washers, bores and other moving surfaces.",
        ],
    )


def commissioning(d):
    header(
        d,
        6,
        "Assemble and measure",
        "Enter the results in the commissioning sheet. A measured pass belongs to that physical build.",
    )
    groups = [
        (
            "1   Check the pieces",
            [
                "Check all six printed rules: two per template and two on the dial. Include reading error.",
                "Gauge every bore, thickness and pin spacing. Check flatness: 0.20 max for a link, 0.30 for",
                "the base. Sight shaft perpendicularity in both directions; total tilt must be within 1 degree.",
            ],
        ),
        (
            "2   Set the ground and fit the moving pins",
            [
                "Fit feet and clamp the O/D shafts to the base without crushing plywood. Check the base sits",
                "flat with under-board hardware clear of the table. Install the ground stacks from page 5.",
                "Fit bottom collars to A/B, pass shafts up through the lower links, then add the middle stacks.",
            ],
        ),
        (
            "3   Build the correct branch",
            [
                "Point the crank from O toward D, so A=(46.6,0) and B is above O-D. B should be near",
                "(13.90,77.48) relative to O. Add coupler above the two lower links. Retain all four stacks",
                "with the specified endplay. Fit the A grip only after checking the dry assembly.",
            ],
        ),
        (
            "4   Check movement and hardware",
            [
                "After curing, mark collars and shafts. Turn ten full cycles each way. Stop on any rub or jam.",
                "At each close approach stop and slip in a blade freely. Require thickness minus error >= 0.50.",
                "Example: 0.55 ±0.02 passing proves at least 0.53 mm. A 0.50 ±0.02 pass gives only 0.48.",
            ],
        ),
        (
            "5   Measure three forward sweeps",
            [
                "Turn CCW from minimum to maximum rocker angle. Log crank and rocker readings at both.",
                "Wrap crank travel to 0-360 degrees: increasing travel / (360 minus increasing travel).",
                "Each swing must span 57-63 degrees; each increasing/return ratio 1.65-1.95, with uncertainty.",
            ],
        ),
        (
            "6   Inspect again and record",
            [
                "Check the ten collars against their witness marks, endplay, base feet and the bonded grip.",
                "Rework and repeat the complete sequence if a limit fails. Record measured numbers, date",
                "and build ID. The reference is qualified by calculation; confirm operation on the finished unit.",
            ],
        ),
    ]
    for i, (title, lines) in enumerate(groups):
        note(d, 15, 53 + i * 34, title, lines)
    d.text(
        15,
        248,
        "Nominal reversals occur near crank 56.61 and 288.08 degrees. Rocker extremes are 101.33 and 161.32.",
        3.45,
        INK,
    )
    d.text(
        15, 255, "At constant crank speed, the crank-angle ratio equals the stroke-time ratio.", 3.25, MUTED
    )


def drawing_files():
    """Profiles use capsule arcs and circles. Coordinates are in mm without tool compensation."""
    for typ in ("base", "links"):
        w, h = (blank_width(), blank_height()) if typ == "base" else (168, 104)
        d = Drawing(
            w + 2,
            h + 2,
            title=f"Linkwright {typ} manufacturing outlines",
            desc="Millimetres, 1 mm margin. Finished bore and outside outlines; nominal toolpath requires machine compensation.",
        )
        doc = ezdxf.new("R2010", setup=False)
        doc.units = units.MM
        for layer, color in (("PROFILE", 7), ("BORE", 3)):
            doc.layers.add(layer, color=color)
        m = doc.modelspace()

        def drill(px, py, r):
            d.circ(1 + px, h + 1 - py, r, "#006B62", None, 0.10)
            m.add_circle((px, py), r, dxfattribs={"layer": "BORE"})

        if typ == "base":
            d.rect(1, 1, w, h, 0, INK, None, 0.1)
            m.add_lwpolyline([(0, 0), (w, 0), (w, h), (0, h)], close=True, dxfattribs={"layer": "PROFILE"})
            for x, y in base_holes():
                drill(x, y, radius_ground())
        else:
            for _, key, first, last, x, y in links_layout():
                length = S[key] + (LINK["pointer_extension"] if key == "rocker" else 0)
                r = LINK["width"] / 2
                d.rect(1 + x - r, 1 + h - y - r, length + 2 * r, 2 * r, r, INK, None, 0.10)
                # Closed polyline with two true semicircles. Four tangencies define the capsule.
                m.add_lwpolyline(
                    [(x, y - r, 0), (x + length, y - r, 1), (x + length, y + r, 0), (x, y + r, 1)],
                    format="xyb",
                    close=True,
                    dxfattribs={"layer": "PROFILE"},
                )
                for xp in [x, x + S[key]]:
                    drill(xp, y, radius_moving())
        d.save_svg(f"cut-{typ}.svg")
        doc.saveas(ROOT / f"cut-{typ}.dxf")


def main():
    assembly_svg()
    drawing_files()
    for name, page in [("a4", A4), ("letter", letter)]:
        document = Canvas(str(ROOT / f"build-{name}.pdf"), pagesize=page, invariant=1)
        document.setTitle("Linkwright reference build | LW-159-R1 | " + name.upper())
        document.setAuthor("Crest-Root / Linkwright")
        for render in (introduction, stock, base_page, links_page, stacks, commissioning):
            render(Drawing(page[0] / mm, page[1] / mm, document))
            document.showPage()
        document.save()


if __name__ == "__main__":
    main()
