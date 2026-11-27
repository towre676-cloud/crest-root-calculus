/* Run-level paired audit of operator-adjudicated trial records. */
(function (root) {
  "use strict";
  function evaluate(records, engine) {
    if (!records.length) throw Error("Select at least one saved trial record.");
    const seen = new Set();
    const rows = [],
      counts = {
        total: 0,
        publishedBench: 0,
        synthetic: 0,
        unlabeled: 0,
        incomplete: 0,
        paired: 0,
        faults: 0,
        clean: 0,
        bothCatch: 0,
        trialOnly: 0,
        baselineOnly: 0,
        bothMiss: 0,
        bothFlagClean: 0,
        trialOnlyFlagClean: 0,
        baselineOnlyFlagClean: 0,
        bothClear: 0,
      };
    for (const record of records) {
      const result = engine.replay(record);
      const key = JSON.stringify([result.setup, result.run]);
      if (seen.has(key))
        throw Error(
          "Duplicate setup and run: " + result.setup + " / " + result.run + ".",
        );
      seen.add(key);
      const base = record.comparison?.baseline,
        truth = record.comparison?.observed;
      const trial = result.status,
        caseRow = {
          setup: result.setup,
          run: result.run,
          source: result.source,
          truth: truth || "Unrecorded",
          baseline: base || "Unrecorded",
          trial,
        };
      counts.total++;
      let category;
      if (result.source === "synthetic") category = "synthetic";
      else if (result.source === "published-bench") category = "publishedBench";
      else if (
        !["Completed within limits", "Confirmed deviation"].includes(truth)
      )
        category = "unlabeled";
      else if (
        !["Passed", "Flagged"].includes(base) ||
        !["violation", "meets"].includes(trial)
      )
        category = "incomplete";
      else {
        counts.paired++;
        const flag = trial === "violation",
          previous = base === "Flagged";
        if (truth === "Confirmed deviation") {
          counts.faults++;
          category = flag
            ? previous
              ? "bothCatch"
              : "trialOnly"
            : previous
              ? "baselineOnly"
              : "bothMiss";
        } else {
          counts.clean++;
          category = flag
            ? previous
              ? "bothFlagClean"
              : "trialOnlyFlagClean"
            : previous
              ? "baselineOnlyFlagClean"
              : "bothClear";
        }
      }
      counts[category]++;
      rows.push({ ...caseRow, category });
    }
    return { unit: "one unique setup and run per record", counts, rows };
  }
  root.TrialAudit = { evaluate };
})(globalThis);
