"use strict";
const { test } = require("node:test");
const a = require("node:assert/strict"),
  fs = require("node:fs"),
  path = require("node:path");
const { analyze, summary, replay } = require("./engine");
require("./examples");
const header = "setup,run,cycle,phase_deg,channel,value\n";
function setup(changes = {}) {
  return {
    schema: 1,
    setup: "R1",
    name: "Bench",
    source: "machine",
    loops: [],
    events: [],
    handoffs: [],
    channels: [
      {
        id: "x",
        label: "Gap",
        unit: "mm",
        sense: "min",
        limit: 0,
        error: 0,
        rate: 1,
        window: [0, 10],
      },
    ],
    ...changes,
  };
}
function calculate(p, points) {
  const csv =
    header +
    points
      .map(([x, v, ch = "x"]) => "R1,run,1," + x + "," + ch + "," + v)
      .join("\n");
  return analyze(p, csv);
}
function check(p, points) {
  return calculate(p, points).cycles[0].checks[0];
}
test("packaged fixtures match the files, validate and replay", () => {
  for (const [key, { csv, setup: p }] of Object.entries(trialExamples)) {
    a.deepEqual(
      p,
      JSON.parse(
        fs.readFileSync(
          path.join(__dirname, "examples", key + ".json"),
          "utf8",
        ),
      ),
    );
    a.equal(
      csv,
      fs.readFileSync(path.join(__dirname, "examples", key + ".csv"), "utf8"),
    );
    const r = analyze(p, csv);
    const result = summary(r);
    a.deepEqual(replay({ engine: r.version, setup: p, csv, result }), r);
    result.status = "forged";
    a.throws(
      () => replay({ engine: r.version, setup: p, csv, result }),
      /disagrees/,
    );
  }
  const evaluate = (k) => analyze(trialExamples[k].setup, trialExamples[k].csv);
  a.equal(evaluate("opener").cycles[0].checks[0].lower, "-4");
  a.equal(evaluate("opener").cycles[1].checks[0].status, "meets");
  a.equal(evaluate("folder").cycles[0].checks[1].upper, "-0.25");
  a.equal(evaluate("reach").loops[0].margin, "-0.00006");
  a.equal(evaluate("reach").loops[0].phase, "180.25");
});
test("threshold equality, strict below threshold, sensor uncertainty, maximum limit, clipped windows", () => {
  const p = setup();
  a.equal(
    check(p, [
      [0, 5],
      [10, 5],
    ]).lower,
    "0",
  );
  a.equal(
    check(p, [
      [0, 5],
      [10, 5],
    ]).status,
    "meets",
  );
  a.equal(
    check(p, [
      [0, 4.999999],
      [10, 5],
    ]).lower,
    "-0.0000005",
  );
  Object.assign(p.channels[0], {
    window: [2, 4],
    sense: "max",
    limit: 10,
    error: 1,
  });
  const r = check(p, [
    [0, 4],
    [10, 4],
  ]);
  a.equal(r.lower, "1");
  a.equal(r.upper, "9");
  Object.assign(p.channels[0], { window: [0, 10], rate: 0 });
  a.equal(
    check(p, [
      [0, 12],
      [10, 12],
    ]).status,
    "violation",
  );
});
test("all observations constrain a single trajectory, including nonadjacent samples", () => {
  const p = setup();
  Object.assign(p.channels[0], { rate: 0, error: 1 });
  a.equal(
    check(p, [
      [0, 0],
      [5, 2],
      [10, 4],
    ]).status,
    "inconsistent",
  );
  a.equal(
    check(p, [
      [0, 0],
      [5, 2],
      [10, 2],
    ]).lower,
    "1",
  );
});
test("missing coverage, missing rate, wrap, re-ordered samples and measured failure", () => {
  const p = setup();
  a.equal(
    check(p, [
      [1, 4],
      [10, 4],
    ]).status,
    "missing",
  );
  a.equal(
    check(p, [
      [1, -4],
      [10, 4],
    ]).status,
    "violation",
  );
  p.channels[0].rate = null;
  a.equal(
    check(p, [
      [0, 9],
      [10, 9],
    ]).status,
    "unresolved",
  );
  a.equal(
    check(p, [
      [0, -1],
      [10, 9],
    ]).status,
    "violation",
  );
  Object.assign(p.channels[0], { rate: 0, window: [350, 10] });
  a.equal(
    check(p, [
      [360, 0],
      [0, 0],
    ]).status,
    "meets",
  );
  a.equal(
    check(p, [
      [350, 0],
      [0, 0],
    ]).status,
    "missing",
  );
});
test("schema and trace identity, missing/sparse events, wide CSV and quoted fields", () => {
  const p = setup({
    events: [
      { id: "on", errorDeg: 1 },
      { id: "off", errorDeg: 1 },
    ],
    handoffs: [{ label: "Handoff", first: "on", second: "off", minDeg: 5 }],
  });
  const wide =
    'setup,run,cycle,phase_deg,x,on,off\r\nR1,run,1,0,"5",1,\r\nR1,run,1,10,5,,1\n';
  const r = analyze(p, wide);
  a.equal(r.cycles[0].checks[1].lower, "3");
  a.deepEqual(
    r,
    calculate(p, [
      [0, 5],
      [0, "", "on"],
      [10, 5],
      [10, "", "off"],
    ]),
  );
  a.throws(() => analyze(p, wide.replace("R1", "R2")), /Setup ID/);
  a.throws(() => analyze(p, wide.replace("x,on,off", "x,on,on")), /duplicate/);
  a.throws(() => analyze(p, wide + "R1,run,1,11,,1,\n"), /Repeated/);
  a.throws(() => analyze(p, wide + "R1,run,1,10,5,,\n"), /Duplicate phase/);
  a.throws(
    () => analyze(p, wide.replace("10,5,,1", "10,5,0,1")),
    /once per event/,
  );
  a.equal(
    calculate(p, [
      [0, 5],
      [10, 5],
    ]).cycles[0].checks[1].status,
    "missing",
  );
  a.equal(
    calculate(p, [
      [0, 5],
      [10, 5],
      [3, "", "on"],
      [7, "", "off"],
    ]).cycles[0].checks[1].status,
    "unresolved",
  );
  a.equal(
    calculate(p, [
      [0, 5],
      [10, 5],
      [7, "", "on"],
      [8, "", "off"],
    ]).cycles[0].checks[1].status,
    "violation",
  );
});
test("bounded-stroke minima include an interior projection on a dimensional edge", () => {
  const l = {
    name: "Test",
    ground: [8, 12],
    crank: [5, 15],
    coupler: [5, 5],
    rocker: [6, 6],
    strokeCos: [0.8, 0.8],
  };
  const g = calculate(setup({ loops: [l] }), [
    [0, 9],
    [10, 9],
  ]).loops[0];
  a.equal(g.folded, "22.04");
  a.equal(g.witness.ground, "8");
  a.equal(g.witness.crank, "6.4");
  for (const scope of [
    [-1, 1],
    [0, 0],
    [-0.63, 0.72],
    [1, 1],
    [-1, -1],
  ]) {
    const r = calculate(setup({ loops: [{ ...l, strokeCos: scope }] }), [
      [0, 9],
      [10, 9],
    ]).loops[0];
    for (let i = 8; i <= 12; i += 0.25)
      for (let j = 5; j <= 15; j += 0.25)
        for (const k of scope) {
          const d = i * i + j * j - 2 * i * j * k;
          a.ok(d - 1 >= Number(r.folded) - 1e-10);
          a.ok(121 - d >= Number(r.extended) - 1e-10);
        }
    const w = r.witness,
      k = Number(r.cosine),
      x = +w.ground,
      y = +w.crank,
      d = x * x + y * y - 2 * x * y * k;
    a.ok(
      Math.abs(
        Number(r.margin) - (r.constraint === "folded" ? d - 1 : 121 - d),
      ) < 1e-9,
    );
  }
});
test("exact bounds match an independent all-observation envelope on 240 seeded cases", () => {
  let seed = 403;
  function random(n) {
    seed = (seed * 48271) % 2147483647;
    return seed % n;
  }
  for (let t = 0; t < 240; t++) {
    const p = setup(),
      k = 1 + random(4),
      error = random(3),
      l = random(8),
      h = l + 1 + random(12);
    Object.assign(p.channels[0], { rate: k, error, window: [l, h] });
    let obs = [[0, random(8)]];
    for (const x of [4, 9, 14, 20]) {
      const [px, y] = obs.at(-1),
        step = Math.min((x - px) * k, 4);
      obs.push([x, y + random(2 * step + 1) - step]);
    }
    const phases = new Set([
      l,
      h,
      ...obs.map((v) => v[0]).filter((x) => l <= x && x <= h),
    ]);
    // A pair of opposite-slope lines meets at a rational multiple of 1/(2k).
    for (let j = l * 2 * k; j <= h * 2 * k; j++) phases.add(j / (2 * k));
    let low = Infinity,
      high = Infinity;
    for (const x of phases) {
      low = Math.min(
        low,
        Math.max(...obs.map(([v, y]) => y - error - k * Math.abs(x - v))),
      );
      high = Math.min(
        high,
        Math.min(...obs.map(([v, y]) => y + error + k * Math.abs(x - v))),
      );
    }
    const r = check(p, obs);
    a.ok(Math.abs(Number(r.lower) - low) < 1e-9);
    a.ok(Math.abs(Number(r.upper) - high) < 1e-9);
  }
});
test("toggle equality, unsupported numeric input, and wrapped phase precision", () => {
  const p = setup({
    loops: [
      {
        name: "toggle",
        ground: [10, 10],
        crank: [5, 5],
        coupler: [8, 8],
        rocker: [7, 7],
        groundDeg: 300,
      },
    ],
  });
  a.equal(
    calculate(p, [
      [0, 9],
      [10, 9],
    ]).loops[0].status,
    "violation",
  );
  a.equal(
    calculate(p, [
      [0, 9],
      [10, 9],
    ]).loops[0].phase,
    "120",
  );
  a.throws(
    () =>
      check(p, [
        [0, 1.1234567],
        [10, 1],
      ]),
    /at most six/,
  );
  a.throws(
    () =>
      check(p, [
        [0, "NaN"],
        [10, 1],
      ]),
    /decimal/,
  );
});
