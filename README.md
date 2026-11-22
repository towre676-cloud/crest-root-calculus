# Crest–Root Calculus

**Mathematics you can actually check.**

Crest–Root Calculus is an open research project about turning difficult mathematical questions into **exact calculations, explicit answers, and evidence that can be checked again**.

This repository contains Python programs, mathematical papers, worked examples, and verification files. The project ranges across algebra and number theory, with recent work focused on finding hidden elliptic-curve structures inside plane quartic curves.

**New to mathematics or programming? You're welcome here.** You don't have to understand the equations to run a demonstration.

## What does it do?

Imagine being handed a complicated equation describing a geometric shape. Sometimes that shape contains a surprising connection to a simpler object called an **elliptic curve**. Finding the connection is hard; proving you've found *all* of them can be harder.

For certain questions, Crest–Root Calculus can:

- **Discover** candidate mathematical structures using calculations over small finite fields.
- **Recover** exact answers instead of relying only on decimal approximations.
- **Verify** candidate answers by substituting them into the original equations.
- **Certify completeness** in supported cases: show not just that an answer works, but that no others were missed.
- **Reproduce research results** from saved, checkable computational evidence.

It is research software, **not** a general-purpose solver for every equation, and not every calculation is fast.

## A real example

One recent question is: *Does a smooth plane quartic curve have a degree-three map to a genus-one curve? How many distinct maps are there, up to the chosen equivalence?*

The Part XXX decision program can analyze an example called `C*` and certify **one geometric degree-three quotient class**. On a different example, it certifies **zero**. This distinction matters: a finite-field scan can show convincing-looking candidates that do not correspond to genuine solutions in characteristic zero.

Here, **degree three** describes the algebraic map, not the power of the original quartic equation. A **certificate** means the result is accompanied by mathematical conditions that the code checks; it does not mean that all of the repository has been formally proved in a proof assistant.

> **Repository setup:** The research code and monograph are distributed in the separately packaged integrated Part XXX archive and have not yet been uploaded to this GitHub repository. The commands below apply to that archive. A GitHub download will support them once the research source tree is added here.

## Try it on Windows (beginner instructions)

### 1. Get the files

On this GitHub repository page, click **Code → Download ZIP**, then extract the ZIP. Alternatively, if you use Git, clone the repository after it has been published.

The commands below assume the extracted folder contains the `crest` directory, which in turn contains `scan_completeness_xxx`.

### 2. Install Python

Install **Python 3.12** from [python.org](https://www.python.org/downloads/). During installation on Windows, enable the option to add Python to PATH if it is offered.

Open **PowerShell** or **Git Bash** in the extracted project folder. Check Python:

```bash
python --version
```

If `python` is not recognized on Windows, try `py -3.12` instead of `python` in the steps below.

### 3. Install the needed packages

The Part XXX release was tested with **SymPy 1.14**, plus NumPy for independent point-counting checks. Install them with:

```bash
python -m pip install "sympy==1.14.0" numpy
```

For an isolated setup, you can use a Python virtual environment, but it is not necessary for this first example.

### 4. Run your first exact calculation

From the extracted project folder:

```bash
cd crest/scan_completeness_xxx
python -B decision/decide.py scan/cstar.json --prime 5
```

The result is JSON (structured text). Look for:

```json
"complete": true,
"lower": 1,
"upper": 1
```

The matching lower and upper bounds certify that this example has exactly **one** geometric degree-three quotient class.

Want to see an example with none?

```bash
python -B decision/decide.py decision/nonlift_quartic.json --prime 13
```

In this example, the result has `"complete": true`, `"lower": 0`, and `"upper": 0`: no geometric degree-three quotient classes survive the exact lifting check.

### 5. Check the packaged research results

From `crest/scan_completeness_xxx`, run:

```bash
python -B verify_release.py --jobs 3
```

This reruns **18 Part XXX verification groups** and compares their exact outputs against the archived receipt files. Passing means these specified calculations reproduce the stored results on your machine; it does not certify every broader mathematical claim without reading its proof.

If you receive an import error, verify that the packages from step 3 installed into the same Python interpreter you're using here.

## What's in the repository?

The integrated research archive is unusually large because it keeps historical work and reproducible evidence instead of only the newest scripts.

| Location | What you'll find |
|---|---|
| `crest/index.html` | The integrated illustrated mathematical monograph, Parts I–XXX, §§1–148. Open it in a browser. |
| `crest/scan_completeness_xxx/` | The most recent equation-only discovery, lifting, decision, and verification programs. |
| `crest/README_PART_XXX.md` | Detailed commands for the newest release. |
| `crest/pointfree_xxix/` | Earlier point-free discovery experiments and checks. |
| `crest/exhaustion_xxviii/` | Earlier geometric classification and exclusion experiments. |
| `crest/scan_completeness_xxx/part_xxx.pdf` | Standalone research paper for §§138–148. |
| Archived JSON receipt files | Expected exact outputs used to check reproducibility. |

> **Note:** In the present archive, the Part XXX quick-start file is `crest/README_PART_XXX.md`. Historical directories contain additional specialized READMEs.

## For curious readers: what's the mathematics?

A **plane quartic** is a curve given by a homogeneous polynomial of degree four in three coordinates, such as

\[
X^4 + Y^4 + Z^4 = 0.
\]

An **elliptic curve** is a smooth genus-one curve with a selected base point. One familiar form looks like

\[
y^2=x^3+ax+b.
\]

The recent program searches for algebraic maps from quartics to genus-one curves, including targets that may not themselves possess a rational point. To avoid overlooking solutions, it combines geometry, finite-field searches, exact reconstruction, local lifting, and global completeness arguments. These steps have different logical roles: finding candidates is not the same thing as proving exhaustiveness.

The monograph also develops other strands of the broader Crest–Root program. Start with its introduction; you do not need to read all 148 sections to run a first example.

## What is proved, and what remains difficult?

The Part XXX work provides a **theoretical terminating decision procedure for geometric degree-three quotient classes of smooth plane quartics over the rationals**, together with concrete certificates for worked examples. However, the finite fields and lifting precision required by the general procedure can become impractically large. A prior twelve-chart elimination procedure also exceeded 1,600 seconds on representative inputs. The fast examples above are not a promise of fast performance on arbitrary curves.

There is also an important distinction between a genus-one quotient **defined over** the rationals and one that has a **rational point**. The latter question is not automatically resolved by discovering its quotient class.

## How can I help?

Beginners can help by reproducing the examples on different computers, filing clear bug reports, improving the installation guide, or making the outputs easier to understand. More experienced contributors can work on runtime improvements, independent implementations, algebraic proof checks, test coverage, and new examples.

When reporting a problem, include the command you ran, the Python version, your operating system, and the exact error message. Avoid posting claims of new mathematical theorems without precise statements and reproducible verification steps.

## Research and citation

The primary technical reference is `crest/index.html`, the integrated monograph through **Part XXX (§148)**. For claims from a particular experiment, consult its corresponding section, script, and receipt. A formal bibliographic citation, contributor list, release DOI, and software license should be added before a public archival release; none are assumed here.

**Project status:** Active research software. Results and evidence have different levels of verification; consult the proof and receipt associated with each claim.