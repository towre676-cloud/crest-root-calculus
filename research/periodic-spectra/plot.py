"""Plot every input, comparing each baseline with its own study's Crest anchor."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter

HERE = Path(__file__).resolve().parent


def main():
    original = json.loads((HERE / "summary.json").read_text())
    original["cases"] = {key: value["methods"] for key, value in original["cases"].items()}
    native = json.loads((HERE / "native-summary.json").read_text())
    panel = json.loads((HERE / "panel.json").read_text())
    keys = [case["id"] for case in panel["cases"]]
    labels = []
    for key in keys:
        if key.startswith("spring-"):
            _, cell, length = key.split("-")
            labels.append(f"Spring · m = {cell[1:]} · L = {length[1:]}")
        else:
            stem, length = key.rsplit("-L", 1)
            name = {"hidden-period":"Hidden period", "large-coefficients":"Large integers",
                    "prime-repeat":"Prime repeats", "single-defect":"Single defect", "aperiodic":"Aperiodic"}[stem]
            labels.append(f"{name} · L = {length}")
    comparisons = [(original, "pari-dense", "General PARI/GP matrix routine", "#B56A1F", (0.7, 3600), [1,10,100,1000]),
                   (native, "pari-native-doubling", "Compiled PARI with the same reduction", "#087F8C", (0.025, 4), [.03,.1,.3,1,3])]
    plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":10, "svg.hashsalt":"crest-periodic-2026-11-23"})
    fig, axes = plt.subplots(1, 2, figsize=(13.8, 11.2), sharey=True)
    fig.patch.set_facecolor("#FBFCFE")
    for ax, (data, method, title, color, limits, ticks) in zip(axes, comparisons):
        ratios = [data["cases"][key][method]["kernel_median_seconds"] /
                  data["cases"][key]["crest-cpp"]["kernel_median_seconds"] for key in keys]
        ax.set_facecolor("#FBFCFE")
        ax.axhspan(11.5, 14.5, color="#EAF0F5", zorder=0)
        ax.axhspan(14.5, 19.5, color="#F4F1ED", zorder=0)
        ax.axvline(1, color="#4A5568", linewidth=1.3, linestyle=(0,(4,3)), zorder=1)
        ax.hlines(range(len(keys)), 1, ratios, color=color, linewidth=1.8, alpha=.55, zorder=2)
        ax.scatter(ratios, range(len(keys)), color=color, s=46, zorder=3)
        for i, value in enumerate(ratios):
            ax.annotate(f"{value:.2f}×" if value < 10 else f"{value:.1f}×", (value, i),
                        xytext=(7,0), textcoords="offset points", va="center", fontsize=9, color="#253047")
        ax.set_xscale("log")
        ax.set_xlim(*limits)
        ax.set_ylim(19.8,-.8)
        ax.xaxis.set_major_locator(FixedLocator(ticks))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x,pos:f"{x:g}×"))
        ax.tick_params(axis="x", which="minor", bottom=False)
        ax.tick_params(axis="y", length=0)
        ax.grid(axis="x", which="major", color="#DDE3EC", linewidth=.65, zorder=0)
        ax.set_title(title, loc="left", pad=25, fontsize=12, fontweight="bold", color="#14243A")
        wins = sum(v>1 for v in ratios)
        ax.text(0, 1.015, f"Crest faster: {wins}/20     PARI faster: {20-wins}/20", transform=ax.transAxes,
                fontsize=10, color="#536275")
        for spine in ax.spines.values(): spine.set_visible(False)
        ax.set_xlabel("PARI CPU time ÷ Crest CPU time (log scale)", labelpad=12, color="#33445B")
    axes[0].set_yticks(range(len(keys)), labels=labels)
    axes[0].tick_params(axis="y", pad=12)
    fig.suptitle("Recognizing structure helps. Native PARI usually wins.", x=.065, y=.98,
                 ha="left", fontsize=20, fontweight="bold", color="#14243A")
    fig.text(.065,.944,"Exact complete characteristic polynomials · 20 fixed inputs · medians of three CPU measurements",
             fontsize=11, color="#536275")
    fig.text(.065,.073,"Above 1× favors Crest; below 1× favors PARI. Each panel uses its own contemporaneous Crest measurements.",
             fontsize=10, color="#33445B")
    fig.text(.065,.050,"Shading: hidden / large-integer / prime-repeat controls; then broken-period controls. All 420 scheduled rows were correct.",
             fontsize=9.5, color="#536275")
    fig.text(.065,.027,"PARI 2.19.1 · The compiled control was designed after the original study; every original input was retained. No independent confirmation claim.",
             fontsize=9.2, color="#536275")
    fig.subplots_adjust(left=.22,right=.97,top=.875,bottom=.15,wspace=.14)
    fig.savefig(HERE / "comparison.png", dpi=160, facecolor=fig.get_facecolor())
    fig.savefig(HERE / "comparison.svg", metadata={"Date":"2026-11-23"}, facecolor=fig.get_facecolor())
    svg = HERE / "comparison.svg"
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")
    plt.close(fig)


if __name__ == "__main__":
    main()
