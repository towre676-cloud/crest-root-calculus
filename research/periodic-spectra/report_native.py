"""Summarize the compiled-library supplement without changing the original study."""
from collections import Counter
import json
import math
from pathlib import Path
import statistics

from audit_native import audit_native

HERE=Path(__file__).resolve().parent


def summary():
    verified=audit_native()
    raw=json.loads((HERE/"native-results.json").read_text())
    protocol=json.loads((HERE/"native-results-protocol.json").read_text())
    panel=json.loads((HERE/"panel.json").read_text())
    cases={}
    for case in panel["cases"]:
        cases[case["id"]]={}
        for method in protocol["methods"]:
            rows=[r for r in raw["runs"] if r["case"]==case["id"] and r["method"]==method]
            item={"statuses":[r["status"] for r in rows]}
            if len(rows)==protocol["replicates"] and all(r["status"]=="complete" for r in rows):
                item.update(kernel_median_seconds=statistics.median(r["kernel_seconds_per_solve"] for r in rows),
                            kernel_min_seconds=min(r["kernel_seconds_per_solve"] for r in rows),
                            kernel_max_seconds=max(r["kernel_seconds_per_solve"] for r in rows),
                            cold_wall_median_seconds=statistics.median(r["cold"]["process_wall_seconds"] for r in rows))
            cases[case["id"]][method]=item
    comparisons={}
    for group,ids in {"all":list(cases),**panel["cohorts"]}.items():
        comparisons[group]={}
        for method in protocol["methods"][1:]:
            ratios=[];cold=[]
            for key in ids:
                own,other=cases[key]["crest-cpp"],cases[key][method]
                if "kernel_median_seconds" in own and "kernel_median_seconds" in other:
                    ratios.append(other["kernel_median_seconds"]/own["kernel_median_seconds"])
                    cold.append(other["cold_wall_median_seconds"]/own["cold_wall_median_seconds"])
            comparisons[group][method]={"paired_inputs":len(ratios),
                    "median_cpu_ratio":statistics.median(ratios) if ratios else None,
                    "geomean_cpu_ratio":math.exp(statistics.mean(math.log(v) for v in ratios)) if ratios else None,
                    "median_cold_wall_ratio":statistics.median(cold) if cold else None,
                    "crest_cpu_wins":sum(v>1 for v in ratios),"pari_cpu_wins":sum(v<1 for v in ratios)}
    return {"design_status":protocol["design_status"],"verification":verified,
            "statuses":{m:dict(Counter(r["status"] for r in raw["runs"] if r["method"]==m)) for m in protocol["methods"]},
            "cases":cases,"comparisons":comparisons}


def render(data):
    doubling=data["comparisons"]["all"]["pari-native-doubling"]
    periodic=data["comparisons"]["periodic-springs"]["pari-native-doubling"]
    crest_wins=[key for key,methods in data["cases"].items()
                if "kernel_median_seconds" in methods["crest-cpp"]
                and "kernel_median_seconds" in methods["pari-native-doubling"]
                and methods["crest-cpp"]["kernel_median_seconds"]
                    < methods["pari-native-doubling"]["kernel_median_seconds"]]
    text=["# Compiled PARI control", "",
          "After observing the original GP comparison, we tested the concrete language/dispatch confound: both the project module and the PARI comparators now run as optimized native executables. "
          "The PARI programs link the same official PARI 2.19.1 / GMP build and implement the same validated period reduction, with linear or binary trace powering. "
          "All 20 original inputs are retained. A contemporaneous Crest rerun is the anchor. "
          "The separate protocol was frozen before its first timing; this is an explicitly post-study control, not a new independently selected confirmation panel.","",
          "## Main finding", "",
          f"**Compiled PARI with binary trace powering is faster on {doubling['pari_cpu_wins']} of {doubling['paired_inputs']} paired inputs**, "
          f"including {periodic['pari_cpu_wins']} of {periodic['paired_inputs']} periodic spring inputs. "
          f"Its median PARI/Crest CPU ratio is {doubling['median_cpu_ratio']:.3f}, and its geometric mean ratio is {doubling['geomean_cpu_ratio']:.3f}. "
          "These are ratios in the stated PARI/Crest direction; no reciprocal-median equivalence is assumed. "
          + ("Crest's compute wins against that comparator are: "+", ".join(f"`{key}`" for key in crest_wins)+". " if crest_wins else "Crest has no compute wins against that comparator. ")
          + "The periodic reduction therefore provides useful structure to exploit, but this native control does not support a general arithmetic-speed advantage over PARI.", "",
          f"The median cold-process wall ratio is {doubling['median_cold_wall_ratio']:.3f}. "
          "Cold execution includes initialization and output, so it answers a different question from repeated polynomial computation. The tables retain both costs.", "",
          "All three implementations receive identical expanded integer arrays. Every timed solve repeats validation, period discovery and complete characteristic-polynomial computation. "
          "Process CPU batches last at least 100 ms; parsing, library initialization and output are included separately in cold process wall time. "
          "The native PARI source and compile command are included. No original observation has been replaced.","",
          "## Completion", "", "| Method | Correct rows | Timeouts | Incorrect outputs | Errors |", "|---|---:|---:|---:|---:|"]
    for method,row in data["statuses"].items():
        text.append(f"| {method} | {row.get('complete',0)} | {row.get('timeout',0)} | {row.get('wrong_answer',0)} | {row.get('error',0)} |")
    text += ["", "## All paired comparisons", "", "Ratios are PARI time divided by Crest time. Values above 1 favor Crest; values below 1 favor PARI. Every input must complete correctly in all three repetitions to be paired.", "",
             "| Cohort | Compiled comparator | Inputs | Median CPU ratio | Geometric mean CPU ratio | Median cold wall ratio | Crest wins | PARI wins |", "|---|---|---:|---:|---:|---:|---:|---:|"]
    for group,methods in data["comparisons"].items():
        for method,row in methods.items():
            if row["paired_inputs"]:
                text.append(f"| {group} | {method} | {row['paired_inputs']} | {row['median_cpu_ratio']:.3f} | {row['geomean_cpu_ratio']:.3f} | {row['median_cold_wall_ratio']:.3f} | {row['crest_cpu_wins']}/{row['paired_inputs']} | {row['pari_cpu_wins']}/{row['paired_inputs']} |")
    text += ["", "## Every input", "", "Median milliseconds per full characteristic polynomial. Per-input ranges and all calibration batches remain in the JSON.", "",
             "| Input | Crest C++ | Native PARI linear | Native PARI doubling |", "|---|---:|---:|---:|"]
    for key,methods in data["cases"].items():
        values=[f"{1000*row['kernel_median_seconds']:.4f}" if "kernel_median_seconds" in row else "/".join(row["statuses"]) for row in methods.values()]
        text.append("| "+key+" | "+" | ".join(values)+" |")
    text += ["", "This compiled comparison removes the main GP-language/native-language mismatch. It still compares these explicit implementations, polynomial representations, algorithms and this machine. "
             "It does not prove that either implementation is optimal or establish a result for general crystallography or chemistry.","",
             "[Protocol](native-results-protocol.json) · [Raw results](native-results.json) · [Summary](native-summary.json) · [Native source](pari_native.c) · [Build provenance](pari-native-build.json) · [Original GP study](results.md)",""]
    return "\n".join(text)


if __name__=="__main__":
    data=summary()
    (HERE/"native-summary.json").write_text(json.dumps(data,indent=2)+"\n")
    (HERE/"native-results.md").write_text(render(data))
    print(json.dumps({"comparisons":data["comparisons"],"statuses":data["statuses"]},indent=2))
