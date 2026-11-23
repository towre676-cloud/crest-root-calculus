"""Regenerate the README illustration with Python and Matplotlib."""

from math import atan, cos, hypot, pi, radians, sin, sqrt
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Arc
from matplotlib.ticker import FuncFormatter

PAPER = "#fcfbf7"
INK = "#192e39"
MUTED = "#475e68"
TEAL = "#00706d"
ORANGE = "#b84c20"
GRID = "#d7ddda"
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 12,
        "text.color": INK,
        "axes.labelcolor": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "svg.fonttype": "none",
        "svg.hashsalt": "crest-root-four-bar",
    }
)


def label(ax, point, name, offset, color=INK):
    ax.annotate(
        name,
        point,
        xytext=offset,
        textcoords="offset points",
        color=color,
        fontweight="bold",
        fontsize=12,
        ha="center",
    )


def assemblies(ax):
    ground, crank, coupler, rocker = 59, 29.3, 52.9, 70
    origin = (0, 0)
    pivot = (ground, 0)
    theta = radians(120)
    pin = (crank * cos(theta), crank * sin(theta))
    dx, dy = ground - pin[0], -pin[1]
    distance = hypot(dx, dy)
    ex, ey = dx / distance, dy / distance
    along = (coupler**2 + distance**2 - rocker**2) / (2 * distance)
    height = sqrt(coupler**2 - along**2)
    ax.plot(*zip(pin, pivot), color="#8e9ba0", ls=(0, (3, 3)), lw=1.5, zorder=1)
    ax.plot(*zip(origin, pivot), color=INK, lw=3, zorder=2)
    ax.plot(*zip(origin, pin), color=INK, lw=4, zorder=3)
    for sign, color, style, name, offset in [
        (1, TEAL, "-", "B₊", (13, 4)),
        (-1, ORANGE, "--", "B₋", (-15, -5)),
    ]:
        joint = (
            pin[0] + along * ex - sign * height * ey,
            pin[1] + along * ey + sign * height * ex,
        )
        ax.plot(*zip(pin, joint, pivot), ls=style, lw=2.6, color=color)
        ax.scatter(*joint, s=52, facecolor=PAPER, edgecolor=color, lw=2, zorder=5)
        label(ax, joint, name, offset, color)
    for position in (origin, pivot, pin):
        ax.scatter(*position, s=54, facecolor=PAPER, edgecolor=INK, lw=2.1, zorder=5)
    for position, name, offset in [
        (origin, "O", (0, -13)),
        (pivot, "D", (12, 2)),
        (pin, "A", (-11, 6)),
    ]:
        label(ax, position, name, offset)
    ax.add_patch(Arc(origin, 18, 18, theta1=0, theta2=120, lw=1, color=MUTED))
    ax.text(7.5, 7, "120°", fontsize=10)
    ax.text(31, -4.7, "59 mm", ha="center", color=MUTED, fontsize=10)
    ax.set(xlim=(-38, 78), ylim=(-31, 77), aspect="equal")
    ax.axis("off")


def closure_window(ax):
    ground, crank, coupler, rocker = 100, 30, 75, 54.99994
    center = 180.25

    def margin(angle):
        distance = hypot(
            ground - crank * cos(radians(angle - 0.25)),
            crank * sin(radians(angle - 0.25)),
        )
        return 1000 * (coupler + rocker - distance)

    angles = [center + j / 1000 for j in range(-310, 311)]
    values = [margin(angle) for angle in angles]
    alpha0 = (coupler + rocker) ** 2 - (ground - crank) ** 2
    alpha1 = (coupler + rocker) ** 2 - (ground + crank) ** 2
    half_width = 2 * atan(sqrt(-alpha1 / alpha0)) * 180 / pi
    left, right = center - half_width, center + half_width
    ax.axvspan(left, right, fc="#f9e5d7", ec="none", zorder=0)
    ax.axhline(0, color=INK, lw=1.2)
    ax.plot(angles, values, lw=2.6, color=INK)
    for angle in (180, 180.5):
        ax.scatter(
            angle, margin(angle), s=60, facecolor=TEAL, edgecolor=PAPER, zorder=6
        )
        ax.annotate(
            "sample",
            (angle, margin(angle)),
            xytext=((18 if angle == 180 else -18), 11),
            textcoords="offset points",
            color=TEAL,
            fontsize=10,
            weight="bold",
            ha="center",
        )
    ax.text(
        center,
        0.21,
        "Cannot close",
        ha="center",
        color=ORANGE,
        weight="bold",
        fontsize=11,
    )
    ax.text(center, 0.169, "0.2613° wide", ha="center", color=ORANGE, fontsize=10)
    ax.scatter(center, margin(center), s=32, color=ORANGE, zorder=4)
    ax.annotate(
        "−0.060 µm",
        (center, margin(center)),
        xytext=(0, -16),
        textcoords="offset points",
        ha="center",
        weight="bold",
        fontsize=10,
        color=ORANGE,
    )
    ax.set(
        xlim=(179.945, 180.555),
        ylim=(-0.13, 0.285),
        yticks=[-0.1, 0, 0.1, 0.2],
        xticks=[180, 180.125, 180.25, 180.375, 180.5],
        xlabel="Crank angle in the fixed frame (degrees)",
        ylabel="Outer closure margin (µm)",
    )
    ax.set_xticklabels(["180.0", "", "180.25", "", "180.5"])
    ax.yaxis.set_major_formatter(
        FuncFormatter(lambda y, _: f"{y:.1f}".replace("-", "−"))
    )
    ax.grid(axis="y", lw=0.7, color=GRID, zorder=0)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    for side in ["bottom", "left"]:
        ax.spines[side].set_color(GRID)
    ax.tick_params(labelsize=10, color="#83969b")
    ax.xaxis.label.set_size(11)
    ax.yaxis.label.set_size(11)


if __name__ == "__main__":
    figure = plt.figure(figsize=(11, 4.25), facecolor=PAPER)
    for x, title, subtitle in [
        (0.055, "One crank angle, two poses", "Linkwright wing at 100% scale"),
        (0.515, "A failure between samples", "A second linkage, checked every 0.5°"),
    ]:
        figure.text(x, 0.906, title, fontsize=17, weight="bold")
        figure.text(x, 0.836, subtitle, fontsize=11, color=MUTED)
    assemblies(figure.add_axes([0.045, 0.145, 0.38, 0.63], facecolor=PAPER))
    closure_window(figure.add_axes([0.579, 0.243, 0.37, 0.515], facecolor=PAPER))
    figure.text(
        0.055,
        0.083,
        "O and D are fixed. Holding A still lets B reflect",
        fontsize=10,
        color=INK,
    )
    figure.text(
        0.055, 0.035, "across the dotted line from A to D.", fontsize=10, color=INK
    )
    figure.text(
        0.515,
        0.083,
        "Above zero, the bars reach. Below zero, a gap opens.",
        fontsize=10,
        color=INK,
    )
    figure.text(
        0.515,
        0.035,
        "Both sampled poses pass, even though the crank jams between them.",
        fontsize=10,
        color=INK,
    )
    target = Path(__file__).with_name("four-bar.svg")
    figure.savefig(
        target,
        metadata={
            "Date": None,
            "Title": "Two four-bar poses and a failure hidden between sampled angles",
        },
    )
    lines = target.read_text(encoding="utf-8").splitlines()
    target.write_text("\n".join(line.rstrip() for line in lines) + "\n", encoding="utf-8")
