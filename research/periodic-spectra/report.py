"""Render all frozen observations; ratios > 1 mean Crest was faster."""
from collections import Counter
import json
import math
from pathlib import Path
import statistics

from audit import audit
from model import coefficient_period, expand_case

HERE=Path(__file__).resolve().parent


def summarize():
    verification=audit(controls=False)
    panel=json.loads((HERE/"panel.json").read_text())
    raw=json.loads((HERE/"results.json").read_text())
    protocol=json.loads((HERE/"results-protocol.json").read_text())
    refs=json.loads((HERE/"results-reference.json").read_text())
    methods=protocol["methods"]
    cases={}
    for case in panel["cases"]:
        a,b=expand_case(case)
        row={"cohort":case["cohort"],"sites":len(a),"actual_period":coefficient_period(a,b),
             "max_coefficient_bits":max(abs(int(v)).bit_length() for v in refs[case["id"]]),"methods":{}}
        for method in methods:
            runs=[r for r in raw["runs"] if r["case"]==case["id"] and r["method"]==method]
            item={"statuses":[r["status"] for r in runs]}
            if len(runs)==protocol["replicates"] and all(r["status"]=="complete" for r in runs):
                values=[r["kernel_seconds_per_solve"] for r in runs]
                item.update(kernel_median_seconds=statistics.median(values),
                            kernel_min_seconds=min(values),kernel_max_seconds=max(values),
                            cold_wall_median_seconds=statistics.median(r["cold"]["process_wall_seconds"] for r in runs),
                            cold_cpu_median_seconds=statistics.median(r["cold"]["process_cpu_seconds"] for r in runs))
            row["methods"][method]=item
        cases[case["id"]]=row
    groups={"all":list(cases),**panel["cohorts"]}
    comparisons={}
    for group,ids in groups.items():
        comparisons[group]={}
        for method in methods[1:]:
            ratios=[];cold=[];paired=[]
            for key in ids:
                a,b=cases[key]["methods"]["crest-cpp"],cases[key]["methods"][method]
                if "kernel_median_seconds" in a and "kernel_median_seconds" in b:
                    ratios.append(b["kernel_median_seconds"]/a["kernel_median_seconds"])
                    cold.append(b["cold_wall_median_seconds"]/a["cold_wall_median_seconds"])
                    paired.append(key)
            comparisons[group][method]={"paired_inputs":len(ratios),"paired_ids":paired,
                "median_cpu_ratio":statistics.median(ratios) if ratios else None,
                "geomean_cpu_ratio":math.exp(statistics.mean(math.log(v) for v in ratios)) if ratios else None,
                "median_cold_wall_ratio":statistics.median(cold) if cold else None,
                "crest_cpu_wins":sum(r>1 for r in ratios),"comparator_cpu_wins":sum(r<1 for r in ratios)}
    return {"scope":protocol["scope"],"pari_version":protocol["pari_version"],
            "replicates":protocol["replicates"],"verification":verification,
            "method_statuses":{m:dict(Counter(r["status"] for r in raw["runs"] if r["method"]==m)) for m in methods},
            "cases":cases,"comparisons":comparisons}


def render(summary):
    text=["# Exact periodic spectra: measured results", "",
          f"PARI/GP **{summary['pari_version']}** (GMP kernel) and the new Crest C++/GMP module were run on the same CPU. "
          f"The frozen panel has **{len(summary['cases'])} inputs**, four implementations and **{summary['replicates']} repetitions**, giving "
          f"**{summary['verification']['planned_rows']} planned rows**. Every returned integer coefficient was compared with an independent full-ring determinant recurrence. "
          "The complete schedule, calibration batches, software hashes and all outcomes are retained.", "",
          "## What is being compared", "",
          "The common input is a cyclic Jacobi operator supplied as integer diagonal and bond arrays. `crest-cpp` detects the least joint period and uses the cell transfer polynomial with a linear Chebyshev recurrence. "
          "`pari-dense` builds the full matrix and calls the recommended default `charpoly`, which uses modular arithmetic for integer matrices. "
          "`pari-periodic-linear` independently implements the same periodic reduction in GP. `pari-periodic-doubling` also uses binary powering identities, exploiting PARI's polynomial multiplication. "
          "These are separate implementation comparisons; the classical transfer/Floquet reduction is not a new theorem of this project.","",
          "Compute timing includes validation, period discovery, construction and complete polynomial computation on every solve; it excludes text parsing and serialization for every method. "
          "The dense baseline's matrix construction is included because it starts from the same arrays. GP's millisecond timer is used only for batches lasting at least 100 ms; the native solver uses its process CPU clock. "
          "Cold wall timing separately measures a single solve in a fresh process, including parsing and output. Compilation and independent verification are excluded. Three repetitions characterize this host and do not establish a population-level confidence interval.","",
          "## Completion", "", "| Method | Correct rows | Timeouts | Incorrect outputs | Process errors |", "|---|---:|---:|---:|---:|"]
    for method,status in summary["method_statuses"].items():
        text.append(f"| {method} | {status.get('complete',0)} | {status.get('timeout',0)} | {status.get('wrong_answer',0)} | {status.get('error',0)} |")
    text += ["", "## Paired comparisons", "", "Ratios divide each comparator's median time by Crest's median time on the same input. **Above 1 means Crest is faster; below 1 means PARI is faster.** "
             "Only inputs correctly completed in every repetition by both implementations enter a finite ratio. The geometric mean is reported alongside the median. "
             "Broken-period cases still have cyclic nearest-neighbor structure; they are not generic dense matrices.", "",
             "| Cohort | Comparator | Paired inputs | Median CPU ratio | Geometric mean CPU ratio | Median cold wall ratio | Crest CPU wins |", "|---|---|---:|---:|---:|---:|---:|"]
    for group,methods in summary["comparisons"].items():
        for method,row in methods.items():
            if row["paired_inputs"]:
                text.append(f"| {group} | {method} | {row['paired_inputs']} | {row['median_cpu_ratio']:.3f} | {row['geomean_cpu_ratio']:.3f} | {row['median_cold_wall_ratio']:.3f} | {row['crest_cpu_wins']}/{row['paired_inputs']} |")
    text += ["", "## Every input", "", "Times are median **milliseconds per complete polynomial** from the batched CPU observations. The full JSON also records the minimum and maximum of the three repetitions, and cold wall/CPU times. No case is omitted.", "",
             "| Input | Sites | Found period | Maximum coefficient bits | Crest C++ | PARI dense | PARI periodic linear | PARI periodic doubling |", "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for key,row in summary["cases"].items():
        cells=[]
        for method,item in row["methods"].items():
            cells.append(f"{1000*item['kernel_median_seconds']:.4f}" if "kernel_median_seconds" in item else "/".join(item["statuses"]))
        text.append(f"| {key} | {row['sites']} | {row['actual_period']} | {row['max_coefficient_bits']} | " + " | ".join(cells) + " |")
    text += ["", "## Interpretation and limits", "",
             "This panel measures synthetic one-dimensional harmonic/signed cyclic operators and exact characteristic polynomials. It does not benchmark full crystallographic space-group recognition, LLL, DFT, measured phonon data, molecular-isomer enumeration or drug discovery. "
             "The native implementation's comparison with GP also includes language/dispatch and data-representation differences; it is not a comparison with every possible hand-written libpari or gp2c implementation. "
             "The structured GP comparators are essential: a gain over a general matrix routine alone does not establish a superior periodic algorithm. "
             "Coefficient equality establishes the same finite spectral polynomial; the physical model's accuracy remains a separate question.","",
             "[Frozen panel](panel.json) · [Protocol](results-protocol.json) · [Raw observations](results.json) · [Exact reference polynomials](results-reference.json) · [Summary JSON](summary.json) · [Proof](proof.md) · [Field comparison and next challenge](opportunities.md)", ""]
    return "\n".join(text)


if __name__=="__main__":
    summary=summarize()
    (HERE/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    (HERE/"results.md").write_text(render(summary))
    print(json.dumps({"comparisons":summary["comparisons"],"statuses":summary["method_statuses"]},indent=2))
