# Branch reconciliation — 23 November 2026

All six research branches have their content integrated into main. The publication branch contains the same monograph already preserved in the research tree, plus an empty deployment marker. This audit records the state immediately before the periodic spectral challenge; later work can advance main without changing these frozen facts.

## Frozen baseline

Repository: [`towre676-cloud/crest-root-calculus`](https://github.com/towre676-cloud/crest-root-calculus). Main commit: [`479938206f125ca62e67897d376a6f5e6b7f4671`](https://github.com/towre676-cloud/crest-root-calculus/commit/479938206f125ca62e67897d376a6f5e6b7f4671). Its complete tree is `c3eea471e925725275ad3c6fb8c1a3d1c6019190`, containing **506 files**. The [machine-readable receipt](branch-reconciliation-2026-11-23.json) records every observed remote ref, all tree hashes, all differing blob hashes and explicit ancestry witnesses.

GitHub branch, pull-request, commit and recursive-tree responses were checked against the fetched local objects. No research branch has a file absent from the baseline main tree.

## Research branches

| Branch | Head | PR | Files | Identical files at the same path in main | Absent from main |
|---|---|---:|---:|---:|---:|
| `research/row-23640-degree-spectra` | `d802bbad561e18edcfe4095a520ff229de55dc12` | [#1](https://github.com/towre676-cloud/crest-root-calculus/pull/1) | 263 | 261 | 0 |
| `research/primitive-maps-and-tower-factors` | `da416135f74bf2f8312a9e7b09642f5b716aa6f2` | [#2](https://github.com/towre676-cloud/crest-root-calculus/pull/2) | 304 | 297 | 0 |
| `research/quartic-elimination-tower-factors` | `0d883e2424c77fb64a4697a9760ae39930be2f32` | [#3](https://github.com/towre676-cloud/crest-root-calculus/pull/3) | 324 | 316 | 0 |
| `research/positive-quartics-geometric-tower` | `aa354c29a6c7f364ec34f6598d12567cbe5d3be9` | [#4](https://github.com/towre676-cloud/crest-root-calculus/pull/4) | 364 | 359 | 0 |
| `research/priority-certificates-and-publication` | `4d36b805cf4a9a9cc815cd187f22710f6d5bf802` | [#5](https://github.com/towre676-cloud/crest-root-calculus/pull/5) | 439 | 433 | 0 |
| `research/introduction-and-finite-map-reconstruction` | `0f9b8ecfae3732b4f97d1bc7577d88eb31454117` | [#6](https://github.com/towre676-cloud/crest-root-calculus/pull/6) | 506 | 506 | 0 |

For PRs **#1–5**, the merge commit records the research head as its second parent, and its entire tree equals that branch head’s tree. Those merge commits lie on the recorded first-parent chain leading to the frozen main commit. Thus both their results and their original historical revisions are integrated.

PR **#6** used a squash commit: its source head is not a recorded merge parent. Nevertheless, its **entire source tree equals the frozen main tree**, including every path, file mode and blob. A branch ancestry indicator therefore does not establish missing content for this final research branch. The previously closed duplicate PR #7 does not require another merge.

### Earlier differences are superseded revisions

The smaller identical-file counts on older branches arise only from later revisions of the main README, research README, dated backlog, a few research-note READMEs, the aggregate research ZIP and `research/replay.py`. The JSON receipt lists each such path and both blob hashes. No differing mathematical implementation module, input panel or certificate remains stranded on these branches. The newer documentation and replay runner incorporate subsequent studies; the historical versions remain available through the integrated merge history.

| PR | Merge commit | Branch and merge tree SHA |
|---|---|---|
| #1 | `813336a1c82b8b63473aebf63fb78818a1df7963` | `974bd28433f963c700f2ec881df9b76303512adb` |
| #2 | `c4c78f13048e673a6d25ec50e0312cbb23d26c48` | `4cbebac45e30a0bf513d9b35fde1ac71e531bc02` |
| #3 | `2c568e2b3a93b591e51ede04129a7124acc8c9d9` | `ddf38db23077e3568c235bd57317a6dfb011c4d8` |
| #4 | `743b23b8c736d8207fbfe384b293bfac3fb8af68` | `83aa81b12ff8e9f71f4074eb6a361551023be38a` |
| #5 | `e9ab97224070ad7e6264c3a1e64c52fe5e9198b2` | `1ce257a3be01ce849e31197598edab0b1b06792a` |
| #6 | `479938206f125ca62e67897d376a6f5e6b7f4671` | `c3eea471e925725275ad3c6fb8c1a3d1c6019190` |

## Publication branch mapping

The `gh-pages` head is `488abf904a1165955526184f002114c47dd24b83`, with tree `b08765434159415b74e83d8da382f60f819f70b6`. It has exactly two files.

| Deployment path | Source preserved on main | Verification |
|---|---|---|
| `index.html` | `research/monograph-publication/index.html` | Identical Git blob `0c54959e797ffd4567bb1e7adc09838e2db42026`; independently compared as raw bytes |
| `.nojekyll` | `research/monograph-publication/.nojekyll` | Empty file copied in this reconciliation; identical empty blob `e69de29bb2d1d6434b8b29ae775ad8c2e48c5391` |

The monograph has **1,031,075 bytes** and SHA-256 `978a72fe6d31b629aa18ca72cd4aa3649fe8a6191e5d8656693861aa047961a0`. The empty marker disables Jekyll processing when this source directory is published at a deployment root. It contributes no mathematical content.

The publication commit intentionally has a two-file tree. Merging that branch wholesale into the research branch would import deletion changes for nonpublication files. The lossless integration is the path mapping above; the monograph itself was already present, and this change adds its empty marker beside the canonical source. Existing branches are retained.

## Reproduce the checks

In a clone containing the frozen objects, inspect the complete branch inventory and final content equality:

```sh
git for-each-ref --format="%(refname:short) %(objectname)" refs/remotes/origin
git rev-parse 479938206f125ca62e67897d376a6f5e6b7f4671^{tree}
git rev-parse 0f9b8ecfae3732b4f97d1bc7577d88eb31454117^{tree}
git diff --exit-code 479938206f125ca62e67897d376a6f5e6b7f4671 0f9b8ecfae3732b4f97d1bc7577d88eb31454117 --
```

Both tree commands must print `c3eea471e925725275ad3c6fb8c1a3d1c6019190`; the diff must be empty and exit successfully. For every earlier row, run `git rev-parse <head-sha>^{tree}` and `git rev-parse <merge-sha>^{tree}` and compare them with the corresponding table entry. `git ls-tree -r --full-tree <head-sha>` enumerates the exact path, mode and blob records used for the per-file comparison.

The audit checkout was shallow. Git ancestry queries can stop at a shallow boundary even when the raw commit objects contain recorded parents. This audit therefore read each raw commit with:

```sh
git cat-file commit 479938206f125ca62e67897d376a6f5e6b7f4671
git cat-file commit e9ab97224070ad7e6264c3a1e64c52fe5e9198b2
```

Continue through the `main_first_parent_lineage_through_pr_1` list in the JSON. Each next commit must occur in the preceding commit’s `parent` header. At merges #1–5 the research branch head must also occur as a parent. With full history, the ordinary equivalent is `git merge-base --is-ancestor <merge-sha> <frozen-main-sha>`; each of the six merge commits must succeed.

Check the publication source mapping independently:

```sh
git rev-parse 488abf904a1165955526184f002114c47dd24b83:index.html
git rev-parse 479938206f125ca62e67897d376a6f5e6b7f4671:research/monograph-publication/index.html
git show 488abf904a1165955526184f002114c47dd24b83:index.html | sha256sum
git show 479938206f125ca62e67897d376a6f5e6b7f4671:research/monograph-publication/index.html | sha256sum
```

The blob commands and the SHA-256 commands must agree respectively. Use the frozen commits in this document when reproducing the audit after main advances. This audit made no remote writes and deleted no branches. The new periodic challenge is a subsequent change based on this verified baseline.
