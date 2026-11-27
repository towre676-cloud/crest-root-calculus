"use strict";
const fs = require("node:fs");
const e = require("./engine.js");
const args = process.argv.slice(2);
try {
  if (args.length !== 2)
    throw Error(
      "Usage: node cli.cjs setup.json trace.csv OR node cli.cjs --replay record.json",
    );
  if (args[0] === "--replay") {
    const record = JSON.parse(fs.readFileSync(args[1], "utf8"));
    const r = e.replay(record);
    console.log(
      "Verified",
      r.setup,
      r.run,
      r.status,
      r.cycles.length + " cycles",
    );
  } else {
    const setup = JSON.parse(fs.readFileSync(args[0], "utf8"));
    const csv = fs.readFileSync(args[1], "utf8");
    const result = e.summary(e.analyze(setup, csv));
    console.log(
      JSON.stringify({ engine: e.VERSION, setup, csv, result }, null, 2),
    );
  }
} catch (error) {
  console.error(error.message);
  process.exitCode = 1;
}
