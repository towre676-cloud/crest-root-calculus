"use strict";
const el = (id) => document.getElementById(id);
const fmt = (x, digits = 0) =>
  x.toLocaleString("en-US", { maximumFractionDigits: digits });
const groups = Object.entries(fieldData.bench);
function cellRow(parent, values) {
  const tr = document.createElement("tr");
  for (const value of values) {
    const td = document.createElement("td");
    if (value instanceof Node) td.append(value);
    else td.textContent = value;
    tr.append(td);
  }
  parent.append(tr);
}
function chart(b) {
  const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
  svg.setAttribute("viewBox", "0 0 650 258");
  svg.setAttribute("role", "img");
  svg.setAttribute(
    "aria-label",
    "Recorded current in milliamps, with retained endpoints and withheld peak over elapsed milliseconds.",
  );
  const data = b.example.filter(([t]) => t >= b.window[0] && t <= b.window[1]);
  const current = data.map((p) => p[1] / 1e6);
  const low = Math.min(...current),
    high = Math.max(...current);
  const pad = Math.max((high - low) * 0.16, 0.1);
  const X = (t) => 66 + (554 * (t - b.window[0])) / (b.window[1] - b.window[0]);
  const Y = (v) => 199 - (160 * (v - low + pad)) / (high - low + 2 * pad);
  function draw(tag, attrs, label) {
    const n = document.createElementNS(svg.namespaceURI, tag);
    for (const [k, v] of Object.entries(attrs)) n.setAttribute(k, v);
    if (label !== undefined) n.textContent = label;
    svg.append(n);
  }
  const label = (x, y, txt, anchor = "middle") =>
    draw(
      "text",
      { x, y, fill: "#dceee5", "font-size": 14, "text-anchor": anchor },
      txt,
    );
  label(12, 22, "mA", "start");
  for (let j = 0; j < 4; j++) {
    const v = low + ((high - low) * j) / 3;
    draw("line", { x1: 63, x2: 626, y1: Y(v), y2: Y(v), stroke: "#40645d" });
    label(54, Y(v) + 5, fmt(v, 2), "end");
  }
  const threshold = (b.peak + b.endpointMax) / 2e6;
  draw("line", {
    x1: 64,
    x2: 624,
    y1: Y(threshold),
    y2: Y(threshold),
    stroke: "#ffb58e",
    "stroke-width": 2,
    "stroke-dasharray": "7 5",
  });
  draw("polyline", {
    points: data.map(([t, v]) => X(t) + "," + Y(v / 1e6)).join(" "),
    fill: "none",
    stroke: "#f2d793",
    "stroke-width": 2,
  });
  for (const [t, value] of data) {
    const x = X(t),
      y = Y(value / 1e6);
    if (b.window.includes(t)) {
      draw("path", {
        d: "M " + x + " " + (y - 6) + " l 6 6 -6 6 -6 -6 Z",
        fill: "#103c3b",
        stroke: "#eef7ef",
        "stroke-width": 2.5,
      });
    } else draw("circle", { cx: x, cy: y, r: 4.6, fill: "#ffdf8b" });
    label(x, 224, fmt((t - b.window[0]) / 1000, 2));
  }
  label(337, 249, "Time since first endpoint (ms)");
  el("wave").replaceChildren(svg);
}
function renderBench() {
  const key = el("group").value,
    b = fieldData.bench[key];
  chart(b);
  el("peak").textContent = fmt(b.peak / 1e6, 2) + " mA";
  el("detail").textContent =
    "The higher endpoint read " +
    fmt(b.endpointMax / 1e6, 2) +
    " mA. The original trace reached " +
    fmt(b.peak / 1e6, 2) +
    " mA at " +
    fmt(b.peakTime / 1e6, 6) +
    " s. The endpoints were " +
    fmt((b.window[1] - b.window[0]) / 1000, 3) +
    " ms apart.";
  const url = "index.html?bench=" + encodeURIComponent(key);
  el("allRead").href = url;
  el("sparseRead").href = url + "&samples=ends";
  el("holdout").replaceChildren();
  for (const [label, v] of [
    ["All readings", b.raw],
    ["Screened intervals", b],
  ])
    cellRow(el("holdout"), [
      label,
      fmt(v.rate),
      fmt(v.hidden),
      fmt(v.outside),
      fmt(v.inconsistent) + " / " + fmt(v.brackets),
    ]);
  const total = (key, raw = false) =>
    groups.reduce((s, [, x]) => s + (raw ? x.raw : x)[key], 0);
  el("pooled").textContent =
    "Across all 20 held-out runs the unscreened calculation excluded " +
    fmt(total("outside", true)) +
    " of " +
    fmt(total("hidden", true)) +
    " withheld readings from its bounds. After screening, " +
    fmt(total("outside")) +
    " of " +
    fmt(total("hidden")) +
    " were outside. These are empirical tests of the rate estimate. They do not establish a physical rate bound.";
  el("quality").textContent =
    "Nineteen raw readings across five recordings exceed the bench-specific 10,000 mA review threshold. The screened calculation omits adjacent training pairs and entire held-out brackets touching those readings. The published values remain intact.";
  el("runs").replaceChildren();
  for (const [i, r] of b.runs.entries()) {
    const link = document.createElement("a");
    link.href = r.url;
    link.textContent = r.file;
    cellRow(el("runs"), [
      link,
      i === 2 ? "Held out" : "Rate estimate",
      fmt(r.rows),
      r.current.map((v) => fmt(v / 1e6, 2)).join(" to "),
      fmt(r.flagged),
    ]);
  }
}
const states = [
  "production",
  "downtime",
  "scheduled_downtime",
  "performance_loss",
  "idle",
];
for (const [k, p] of Object.entries(fieldData.piade)) {
  cellRow(el("piade"), [
    k,
    fmt(p.intervals),
    ...states.map((s) => fmt((p.duration_ms[s] || 0) / 36e5, 1)),
  ]);
  el("machine").add(new Option(k, k));
}
for (const [k, v] of Object.entries(fieldData.comopi))
  cellRow(el("comopi"), [
    k,
    ...[
      "sensorRows",
      "alarmRows",
      "warningRows",
      "targetRows",
      "alarms53",
      "alarms54",
    ].map((s) => fmt(v[s])),
  ]);
function scenario() {
  const inputs = ["share", "value", "cost"].map((id) => el(id));
  if (inputs.some((input) => input.value === "" || !input.validity.valid)) {
    el("economics").textContent =
      "Enter a recoverable share, hourly value and trial cost to calculate a scenario.";
    return;
  }
  const [share, value, cost] = inputs.map((i) => i.valueAsNumber);
  const hours =
    ((fieldData.piade[el("machine").value].duration_ms.downtime / 36e5) *
      share) /
    100;
  const net = hours * value - cost;
  el("economics").textContent = Number.isFinite(net)
    ? fmt(hours, 2) +
      " recovered hours × $" +
      fmt(value, 2) +
      "/h − $" +
      fmt(cost, 2) +
      " trial cost = $" +
      fmt(net, 2) +
      " net for this observation period."
    : "Enter finite scenario values.";
}
for (const id of ["share", "value", "cost", "machine"])
  el(id).addEventListener("input", scenario);
for (const [k] of groups)
  el("group").add(
    new Option(
      k
        .replace("empty/", "Unloaded · ")
        .replace("honey/", "Honey stirrer · ")
        .replace("_", " "),
      k,
    ),
  );
el("group").value = "empty/Sine_0";
el("group").onchange = renderBench;
renderBench();
scenario();

const alarms = Object.entries(fieldData.alpi);
const totalAlarms = alarms.reduce((n, [, x]) => n + x.rows, 0);
el("alarmSummary").textContent =
  fmt(totalAlarms) +
  " logged alarm events and " +
  fieldData.alpiCodes +
  " distinct codes across " +
  alarms.length +
  " machines. Each event has a timestamp, machine and alarm code. Repeated events can belong to the same stop.";
for (const [serial, v] of alarms)
  cellRow(el("alpi"), [
    serial,
    fmt(v.rows),
    v.distinct,
    v.top,
    fmt((100 * v.topCount) / v.rows, 1) + "%",
    v.first.slice(0, 10),
    v.last.slice(0, 10),
  ]);
let audit,
  auditGeneration = 0;
el("auditFiles").onchange = async (event) => {
  const generation = ++auditGeneration,
    files = [...event.target.files];
  audit = null;
  el("auditResult").hidden = true;
  el("auditStatus").textContent = "Replaying selected records…";
  try {
    if (files.reduce((n, f) => n + f.size, 0) > 64000000)
      throw Error("Choose a batch up to 64 MB.");
    const records = [];
    for (const f of files) records.push(JSON.parse(await f.text()));
    if (generation !== auditGeneration) return;
    audit = TrialAudit.evaluate(records, Linkwright);
    const c = audit.counts;
    el("auditStatus").textContent =
      c.paired +
      " paired runs of " +
      c.total +
      " submitted · " +
      c.synthetic +
      " synthetic · " +
      c.publishedBench +
      " published bench · " +
      c.unlabeled +
      " lack a definitive bench outcome · " +
      c.incomplete +
      " lack two decisive inspection results.";
    el("auditTable").replaceChildren();
    cellRow(el("auditTable"), [
      "Confirmed deviation (" + c.faults + ")",
      c.bothCatch,
      c.trialOnly,
      c.baselineOnly,
      c.bothMiss,
    ]);
    cellRow(el("auditTable"), [
      "Completed within limits (" + c.clean + ")",
      c.bothFlagClean,
      c.trialOnlyFlagClean,
      c.baselineOnlyFlagClean,
      c.bothClear,
    ]);
    el("auditResult").hidden = false;
  } catch (err) {
    if (generation === auditGeneration)
      el("auditStatus").textContent = err.message;
  } finally {
    if (generation === auditGeneration) el("auditFiles").value = "";
  }
};
el("auditSave").onclick = () => {
  if (!audit) return;
  const url = URL.createObjectURL(
    new Blob([JSON.stringify(audit, null, 2)], { type: "application/json" }),
  );
  const a = document.createElement("a");
  a.href = url;
  a.download = "linkwright-paired-audit.json";
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 5000);
};
