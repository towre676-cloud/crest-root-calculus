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
require("./field-data");
require("./bench");
require("./audit");
test("twenty measured windows reproduce the withheld overshoot; endpoints alone remain unresolved", () => {
  for (const [key, b] of Object.entries(fieldData.bench)) {
    const both = [false, true].map((sparse) => {
      const { setup: p, csv } = BenchTrial.fixture(key, sparse);
      const r = analyze(p, csv);
      a.equal(r.axis, "time_s");
      a.equal(r.loops[0].folded, "21");
      a.equal(r.loops[0].extended, "39");
      return r.cycles[0].checks[0];
    });
    a.equal(both[0].status, "violation");
    a.equal(both[1].status, "unresolved");
    a.ok(Number(both[1].sampled) > 0);
    a.ok(!Object.hasOwn(both[0], "phase"));
    a.ok(
      b.example.some(
        ([t, y]) => t / 1e6 === Number(both[0].time_s) && y === b.peak,
      ),
    );
  }
});
test("time coordinates retain precision, long recordings and edge units; historical records still replay", () => {
  const p = setup({
    schema: 2,
    axis: "time_s",
    events: [
      { id: "on", errorS: "0.001" },
      { id: "off", errorS: "0.001" },
    ],
    handoffs: [{ label: "Transfer", first: "on", second: "off", minS: 0.01 }],
  });
  p.channels[0].window = [400, 401];
  const csv =
    "setup,run,cycle,time_s,x,on,off\nR1,r,c,400,4,1,\nR1,r,c,400.009,, ,1\nR1,r,c,401,4,,\n";
  const r = analyze(p, csv);
  a.equal(r.cycles[0].checks[0].lowerWitness.seconds, 400.5);
  a.equal(r.cycles[0].checks[1].unit, "s");
  a.equal(r.cycles[0].checks[1].lower, "-0.003");
  a.throws(() => analyze(p, csv.replace("time_s", "phase_deg")), /time_s/);
  p.channels[0].window = [401, 400];
  a.throws(() => analyze(p, csv), /increase/);
  const ex = trialExamples.opener,
    result = summary(analyze(ex.setup, ex.csv));
  result.version = "1.0.0";
  a.equal(replay({ ...ex, engine: "1.0.0", result }).version, "2.0.0");
  result.status = "changed";
  a.throws(() => replay({ ...ex, engine: "1.0.0", result }), /disagrees/);
});
test("native import uses selected column and original increasing timestamps", () => {
  const s =
    "% Time [s], Goal Current [mA], Pres. Current [mA], Pres. Voltage [V], Pres. Position [rad], Pres. Velocity [rad/s]\n0.001,1,9,12.1,1,-3\n0.005,1,11,12.1,0,-4";
  const ex = BenchTrial.fromNative(
    s,
    "native.txt",
    "current",
    10,
    "max",
    0,
    null,
  );
  a.equal(analyze(ex.setup, ex.csv).cycles[0].checks[0].sampled, "-1");
  a.throws(
    () =>
      BenchTrial.fromNative(
        s.replace("0.005", "0.001"),
        "f",
        "position",
        0,
        "min",
        0,
        null,
      ),
    /increase/,
  );
  a.throws(
    () =>
      BenchTrial.fromNative(
        s.replace("Time", "Phase"),
        "f",
        "current",
        10,
        "max",
        0,
        null,
      ),
    /header/,
  );
});
test("data totals preserve source counts, holdout membership and excluded readings", () => {
  const b = Object.values(fieldData.bench),
    total = (xs, f) => xs.reduce((s, x) => s + f(x), 0);
  a.equal(b.length, 20);
  a.equal(
    total(b, (x) => total(x.runs, (r) => r.rows)),
    905885,
  );
  a.equal(
    total(b, (x) => total(x.runs, (r) => r.flagged)),
    19,
  );
  a.equal(
    total(b, (x) => x.raw.hidden),
    226452,
  );
  a.equal(
    total(b, (x) => x.hidden),
    226431,
  );
  a.equal(
    total(b, (x) => x.raw.outside),
    5,
  );
  a.equal(
    total(b, (x) => x.outside),
    4,
  );
  a.equal(
    total(Object.values(fieldData.piade), (x) => x.intervals),
    429394,
  );
  a.equal(
    total(Object.values(fieldData.comopi), (x) => x.targetRows),
    41,
  );
  a.equal(
    total(Object.values(fieldData.alpi), (x) => x.rows),
    444834,
  );
});
test("paired audit handles abstentions, clean runs, duplicate identities, source and tampering", () => {
  const engine = require("./engine");
  function make(id, reading, base, truth, source = "machine", rate = 0) {
    const p = setup({ source });
    p.channels[0].rate = rate;
    const csv =
      header +
      "R1," +
      id +
      ",1,0,x," +
      reading +
      "\nR1," +
      id +
      ",1,10,x," +
      reading +
      "\n";
    return {
      engine: engine.VERSION,
      setup: p,
      csv,
      result: summary(analyze(p, csv)),
      comparison: { baseline: base, observed: truth },
    };
  }
  const yes = "Confirmed deviation",
    no = "Completed within limits";
  const records = [
    make("a", -1, "Passed", yes),
    make("b", 1, "Flagged", yes),
    make("c", -1, "Passed", no),
    make("d", 0, "Passed", no),
    make("e", 1, "Passed", no, "machine", null),
    make("f", 1, "Passed", "Unrecorded"),
    make("g", -1, "Passed", yes, "synthetic"),
  ];
  const b = BenchTrial.fixture("empty/Sine_0", false);
  records.push({
    ...b,
    engine: engine.VERSION,
    result: summary(analyze(b.setup, b.csv)),
  });
  const c = TrialAudit.evaluate(records, engine).counts;
  for (const key of [
    "trialOnly",
    "baselineOnly",
    "trialOnlyFlagClean",
    "bothClear",
    "incomplete",
    "unlabeled",
    "publishedBench",
    "synthetic",
  ])
    a.equal(c[key], 1);
  a.equal(c.paired, 4);
  a.throws(
    () => TrialAudit.evaluate([records[0], records[0]], engine),
    /Duplicate/,
  );
  records[0].result.status = "meets";
  a.throws(() => TrialAudit.evaluate(records, engine), /disagrees/);
});
