"""Draw the symmetric six-bar in its theta=0 pose, in exact proportion."""

from html import escape
from math import cos, hypot, radians, sin, sqrt
from pathlib import Path


def triangle(ground_angle, multiplier, elbow_sign):
    """Intersect coupler and rocker circles at the selected assembly position."""
    radius = 49 * multiplier
    crank = (radius, 0)
    ground = (radius * cos(ground_angle), radius * sin(ground_angle))
    dx, dy = ground[0] - crank[0], ground[1] - crank[1]
    distance = hypot(dx, dy)
    ux, uy = dx / distance, dy / distance
    coupler, rocker = 47 * multiplier, 24 * multiplier
    along = (distance**2 + coupler**2 - rocker**2) / (2 * distance)
    across = elbow_sign * sqrt(coupler**2 - along**2)
    elbow = (crank[0] + along * ux - across * uy, crank[1] + along * uy + across * ux)
    assert abs(hypot(elbow[0] - ground[0], elbow[1] - ground[1]) - rocker) < 1e-10
    assert abs(hypot(elbow[0] - crank[0], elbow[1] - crank[1]) - coupler) < 1e-10
    return ground, crank, elbow


def screen(point):
    x, y = point
    return (108 + 1.27 * x, 268 - 1.27 * y)


def xy(point):
    x, y = screen(point)
    return f"{x:.3f},{y:.3f}"


def line(start, end, color, width):
    x1, y1 = screen(start)
    x2, y2 = screen(end)
    return (
        f'<line x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" '
        f'y2="{y2:.3f}" stroke="{color}" stroke-width="{width}" '
        'stroke-linecap="round"/>'
    )


def label(x, y, word, style="label", anchor="start"):
    return f'<text x="{x}" y="{y}" class="{style}" text-anchor="{anchor}">{escape(word)}</text>'


def build():
    origin = (0, 0)
    upper, upper_pin, upper_elbow = triangle(radians(60), 3, -1)
    lower, lower_pin, lower_elbow = triangle(radians(-60), 5, 1)
    joints = [origin, upper, upper_elbow, upper_pin, lower, lower_elbow, lower_pin]
    frame = "#506477"
    crank = "#172f54"
    coupling = "#087e75"
    rocker = "#ac4723"
    parts = [
        """<svg xmlns="http://www.w3.org/2000/svg" width="810" height="623" viewBox="0 0 810 623" role="img" aria-labelledby="title desc">
<title id="title">Geometry of the degree-three six-bar example</title>
<desc id="desc">True-proportion planar diagram at crank angle zero. A rigid ternary crank shares an origin with the fixed triangular frame. Its two pins drive couplers of length 141 and 235 and rockers of length 72 and 120. The fixed pivot rays are 60 degrees above and below the crank.</desc>
<style>
text{font-family:system-ui,sans-serif;fill:#152d43}.title{font-size:23px;font-weight:730}.sub{font-size:14px;fill:#405871}.label{font-size:14px;font-weight:650;paint-order:stroke;stroke:#fff;stroke-width:4px;stroke-linejoin:round}.note{font-size:14px}.head{font-size:15px;font-weight:700}.mini{font-size:12px;fill:#466078}.dim{font-size:13px;fill:#304d64;paint-order:stroke;stroke:#fff;stroke-width:4px}
</style>
<rect x=".5" y=".5" width="809" height="622" rx="10" fill="white" stroke="#ccd5dc"/>"""
    ]
    parts.extend(
        [
            label(25, 36, "Two loops on one crank", "title"),
            label(
                25,
                60,
                "One real assembly at crank angle 0° · all centre distances use the same unit",
                "sub",
            ),
            f'<polygon points="{xy(origin)} {xy(upper)} {xy(lower)}" fill="#eef2f5" stroke="{frame}" stroke-width="3"/>',
        ]
    )
    parts.append(line(origin, lower_pin, crank, 8))
    for fixed, pin, elbow in [
        (upper, upper_pin, upper_elbow),
        (lower, lower_pin, lower_elbow),
    ]:
        parts.extend([line(pin, elbow, coupling, 7), line(fixed, elbow, rocker, 7)])
    for joint in joints:
        x, y = screen(joint)
        parts.append(
            f'<circle cx="{x:.3f}" cy="{y:.3f}" r="5.8" fill="white" stroke="#172f42" stroke-width="2.5"/>'
        )
    parts.extend(
        [
            label(91, 266, "O", anchor="end"),
            label(190, 98, "G₁", anchor="end"),
            label(255, 551, "G₂", anchor="end"),
            label(314, 197, "141"),
            label(253, 92, "72"),
            label(431, 422, "235"),
            label(329, 562, "120"),
            label(140, 169, "147", "dim"),
            label(177, 411, "245", "dim"),
            label(234, 367, "343", "dim"),
            label(296, 254, "P₁"),
            label(425, 255, "P₂"),
            label(310, 82, "J₁"),
            label(430, 579, "J₂"),
        ]
    )
    for x1, x2, caption in [
        (108, screen(upper_pin)[0], "147"),
        (screen(upper_pin)[0], screen(lower_pin)[0], "98"),
    ]:
        parts.append(
            f'<path d="M{x1:.2f},291v8 m0,-4H{x2:.2f} m0,-4v8" fill="none" stroke="#708293"/>'
        )
        parts.append(label((x1 + x2) / 2, 313, caption, "dim", "middle"))
    parts.append('<line x1="489" y1="88" x2="489" y2="580" stroke="#d5dfe6"/>')
    parts.extend([label(515, 111, "Six rigid pieces", "head")])
    legend = [
        (frame, "Ground frame", "O, G₁ and G₂ are fixed."),
        (crank, "One crank", "O, P₁ and P₂ turn together."),
        (coupling, "Two couplers", "P₁–J₁ and P₂–J₂."),
        (rocker, "Two rockers", "J₁–G₁ and J₂–G₂."),
    ]
    for i, (color, title, caption) in enumerate(legend):
        y = 146 + 57 * i
        parts.extend(
            [
                f'<line x1="516" y1="{y-4}" x2="549" y2="{y-4}" stroke="{color}" stroke-width="5" stroke-linecap="round"/>',
                label(560, y, title, "note"),
                label(560, y + 20, caption, "mini"),
            ]
        )
    parts.extend(
        [
            label(515, 410, "Ground pivot directions", "head"),
            label(515, 437, "OG₁ is 60° above the crank.", "note"),
            label(515, 460, "OG₂ is 60° below the crank.", "note"),
            label(515, 511, "Length scale", "head"),
            '<path d="M517,538v9m0,-4.5h127m0,-4.5v9" fill="none" stroke="#28435c" stroke-width="2"/>',
            label(580.5, 566, "100 units", "mini", "middle"),
            "</svg>",
        ]
    )
    return "\n".join(parts) + "\n"


if __name__ == "__main__":
    Path(__file__).with_name("six-bar.svg").write_text(build())
