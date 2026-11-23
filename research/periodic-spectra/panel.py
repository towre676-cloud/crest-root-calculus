"""Deterministic, physically interpretable and adversarial input panel.

Freeze panel.json before timing. No instance is selected by observed speed.
"""
import json
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
SEED = 76126


def build():
    rng = random.Random(SEED)
    cases = []
    cells = [([1], [0]), ([1, 3], [1, 2]),
             ([1, 4, 2], [0, 2, 1]), ([2, 5, 1, 3], [1, 0, 3, 2])]
    for stiffness, pinning in cells:
        m = len(stiffness)
        for length in (24, 96, 384):
            cases.append({"id": f"spring-m{m}-L{length}", "cohort": "periodic-springs",
                          "diagonal": [stiffness[i-1]+stiffness[i]+pinning[i] for i in range(m)],
                          "bonds": [-k for k in stiffness], "repeats": length//m,
                          "physical_model": {"stiffness": stiffness, "pinning": pinning,
                                             "mass": 1}})
    # A full-cell declaration deliberately hides the smaller period from metadata.
    cases.append({"id": "hidden-period-L192", "cohort": "structure-controls",
                  "diagonal": [5, 8, 6]*64, "bonds": [-1, -3, -2]*64, "repeats": 1})
    # Integer coefficient growth, a prime number of cells and a signed-bond control.
    cases.append({"id": "large-coefficients-L96", "cohort": "structure-controls",
                  "diagonal": [2000004, 5000013, 4000011],
                  "bonds": [-1000001, -2000003, -1000007], "repeats": 32})
    cases.append({"id": "prime-repeat-L183", "cohort": "structure-controls",
                  "diagonal": [3, -2, 5], "bonds": [1, -2, 3], "repeats": 61})
    for length in (24, 96, 384):
        cases.append({"id": f"single-defect-L{length}", "cohort": "broken-period",
                      "diagonal": [3, 5, 6], "bonds": [-1, -2, -3], "repeats": length//3,
                      "defect": {"index": length//2, "diagonal_delta": 1}})
    for length in (24, 96):
        cases.append({"id": f"aperiodic-L{length}", "cohort": "broken-period",
                      "diagonal": [rng.randint(-7, 11) for _ in range(length)],
                      "bonds": [rng.choice([-5,-3,-2,-1,1,2,3,5]) for _ in range(length)],
                      "repeats": 1})
    return {"schema": 1, "seed": SEED, "scope": "Full integer characteristic polynomial of a cyclic Jacobi matrix",
            "cohorts": {name: [c["id"] for c in cases if c["cohort"] == name]
                        for name in ("periodic-springs", "structure-controls", "broken-period")},
            "cases": cases}


if __name__ == "__main__":
    target = HERE / "panel.json"
    if target.exists():
        raise SystemExit("Refusing to replace a frozen panel; choose a new study directory.")
    target.write_text(json.dumps(build(), indent=2) + "\n")
    print(f"Wrote {len(build()['cases'])} inputs to {target}")
