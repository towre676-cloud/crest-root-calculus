"use strict";
const $ = (id) => document.getElementById(id);
const statuses = {
  meets: "Meets requirement",
  unresolved: "More evidence needed",
  violation: "Requirement violated",
  missing: "Coverage missing",
  inconsistent: "Trace contradicts bounds",
};
const explain = {
  meets:
    "Every compatible trajectory meets the threshold throughout the required interval.",
  violation:
    "At least one phase must violate the threshold under the declared error and rate bounds.",
  unresolved:
    "All readings can fit the model, but compatible paths allow both outcomes. A reading near the limiting phase can narrow the result.",
  missing:
    "The recorded phase range does not cover the required interval. Add the missing readings.",
  inconsistent:
    "The declared sensor error and rate limit cannot explain these readings together. Check phase alignment, units and calibration.",
};
let active = null;
function node(tag, text, cls) {
  const e = document.createElement(tag);
  if (text !== undefined) e.textContent = text;
  if (cls) e.className = cls;
  return e;
}
function badge(status) {
  return node("span", statuses[status], "pill " + status);
}
function signed(x) {
  if (x == null) return "—";
  return Number(x) > 0 ? "+" + x : String(x);
}
function margins(x) {
  return x.lower === null || x.lower === undefined
    ? "—"
    : signed(x.lower) + " to " + signed(x.upper);
}
function options(el, pairs) {
  el.replaceChildren(
    ...pairs.map(([value, label]) => new Option(label, value)),
  );
}
function save(name, type, body) {
  const url = URL.createObjectURL(new Blob([body], { type }));
  const a = node("a");
  a.href = url;
  a.download = name;
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 5000);
}
function dirty(msg) {
  $("report").hidden = true;
  $("export").disabled = true;
  $("src").textContent = "Unanalyzed inputs";
  $("notice").textContent =
    msg || "Inputs changed. Analyze them to refresh the report.";
}
function activate(setup, csv, selector) {
  const result = Linkwright.analyze(setup, csv);
  $("notes").value = "";
  $("baseline").selectedIndex = 0;
  $("observed").selectedIndex = 0;
  active = { setup, csv, result, selector };
  $("setupText").value = JSON.stringify(setup, null, 2);
  $("traceText").value = csv;
  $("case").value = selector;
  $("report").hidden = false;
  $("export").disabled = false;
  $("notice").textContent = "";
  $("src").textContent =
    setup.source === "synthetic"
      ? "Synthetic reference data"
      : setup.source === "published-bench"
        ? "Published four-bar measurements"
        : "Machine measurements";
  $("axisLabel").textContent =
    setup.axis === "time_s" ? "Elapsed time" : "Encoder phase";
  $("job").textContent = setup.name;
  $("basis").textContent =
    typeof setup.provenance?.basis === "string" ? setup.provenance.basis : "";
  $("overall").replaceWith(
    Object.assign(badge(result.status), { id: "overall" }),
  );
  $("meta").textContent = [
    result.setup,
    "run " + result.run,
    result.cycles.length +
      " recorded cycle" +
      (result.cycles.length === 1 ? "" : "s"),
  ].join("  /  ");
  options(
    $("cycle"),
    result.cycles.map((c, i) => [i, c.cycle]),
  );
  options(
    $("channel"),
    setup.channels.map((c, i) => [i, c.label]),
  );
  const tbl = $("all");
  tbl.replaceChildren();
  for (const cycle of result.cycles)
    for (const x of cycle.checks) {
      const row = node("tr");
      for (const v of [
        cycle.cycle,
        x.label,
        statuses[x.status],
        margins(x),
        x.unit,
      ])
        row.append(node("td", v));
      tbl.append(row);
    }
  const loops = $("loops");
  loops.replaceChildren();
  if (!result.loops.length)
    loops.append(node("p", "This setup defines process checks only."));
  for (const l of result.loops) {
    const item = node("div", undefined, "check");
    item.append(
      node("strong", l.label),
      badge(l.status),
      node(
        "p",
        "Worst " +
          l.constraint +
          " margin " +
          signed(l.margin) +
          " " +
          l.unit +
          (l.cosine != null
            ? " at crank cosine " + l.cosine
            : " at " + l.phase + "°") +
          ". Folded " +
          signed(l.folded) +
          " · extended " +
          signed(l.extended) +
          " " +
          l.unit +
          ".",
      ),
      node(
        "p",
        "Attaining dimensions (mm): " +
          Object.entries(l.witness)
            .map(([k, v]) => k + " " + v)
            .join(" / "),
      ),
    );
    loops.append(item);
  }
  renderCycle();
}
function renderCycle() {
  const c = active.result.cycles[Number($("cycle").value)];
  $("checks").replaceChildren();
  for (const x of c.checks) {
    const div = node("div", undefined, "check");
    const location = x.time_s ?? x.phase;
    const phase =
      location == null
        ? ""
        : " · " +
          Number(Number(location).toFixed(6)) +
          (active.setup.axis === "time_s" ? " s" : "°");
    div.append(
      node("strong", x.label),
      badge(x.status),
      node("p", "Minimum margin: " + margins(x) + " " + x.unit + phase),
    );
    if (x.lower == null && x.upper != null)
      div.append(
        node("p", "Available upper bound " + signed(x.upper) + " " + x.unit),
      );
    $("checks").append(div);
  }
  renderGraph();
}
function renderGraph() {
  const idx = Number($("channel").value),
    spec = active.setup.channels[idx];
  const result = active.result.cycles[Number($("cycle").value)].checks[idx];
  $("graph").replaceChildren();
  if (!spec) {
    $("chartTitle").textContent = "Event timing";
    $("explain").textContent =
      "Timestamp bounds for this cycle appear under Stroke checks.";
    return;
  }
  $("chartTitle").textContent = spec.label;
  $("explain").textContent =
    (active.setup.axis === "time_s"
      ? explain[result.status].replaceAll("phase", "time")
      : explain[result.status]) +
    (result.sampled != null
      ? " Lowest recorded margin: " +
        signed(result.sampled) +
        " " +
        spec.unit +
        "."
      : "");
  if (result.unbounded)
    $("explain").textContent =
      "No rate bound supplied. Measurements can identify a violation at a recorded point. Between-sample margins remain open.";
  if (!result.points?.length) return;
  const [start, end] = spec.window.map(Number),
    full = $("view").value === "all" || start > end;
  const timed = active.setup.axis === "time_s";
  let from = full ? (timed ? result.points[0][0] : 0) : start,
    to = full ? (timed ? result.points.at(-1)[0] : 360) : end;
  if (from === to) to = from + (timed ? 0.001 : 1);
  const width = 690,
    height = 250,
    left = 65,
    right = 674,
    top = 22,
    bottom = 203;
  const k = Number(spec.rate),
    data = result.points;
  // Intersect the two endpoint cones. Each side changes slope at most once per segment.
  const upper = [],
    lower = [];
  function append(which, x, y) {
    which.push([x, y]);
  }
  for (let i = 1; !result.unbounded && i < data.length; i++) {
    const [x, , lo, hi] = data[i - 1],
      [z, , lo2, hi2] = data[i];
    const a = Math.max(x, from),
      b = Math.min(z, to);
    if (a > b) continue;
    const loAt = (t) => Math.max(lo - k * (t - x), lo2 - k * (z - t));
    const hiAt = (t) => Math.min(hi + k * (t - x), hi2 + k * (z - t));
    for (const [list, cross, evalAt] of [
      [lower, k ? ((lo - lo2) / k + x + z) / 2 : a, loAt],
      [upper, k ? ((hi2 - hi) / k + x + z) / 2 : a, hiAt],
    ]) {
      append(list, a, evalAt(a));
      if (cross > a && cross < b) append(list, cross, evalAt(cross));
      append(list, b, evalAt(b));
    }
  }
  const dots = data.filter(([x]) => x >= from && x <= to);
  let lo = 0,
    hi = 0;
  for (const list of [lower, upper, dots])
    for (const [, v] of list) {
      lo = Math.min(lo, v);
      hi = Math.max(hi, v);
    }
  const pad = Math.max((hi - lo) * 0.1, 1e-7);
  lo -= pad;
  hi += pad;
  const X = (t) => left + ((t - from) * (right - left)) / (to - from),
    Y = (v) => bottom - ((v - lo) * (bottom - top)) / (hi - lo);
  const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
  svg.setAttribute("viewBox", "0 0 " + width + " " + height);
  svg.setAttribute("role", "img");
  svg.setAttribute(
    "aria-label",
    spec.label +
      " margin in " +
      spec.unit +
      (timed ? " versus elapsed seconds" : " versus encoder phase"),
  );
  function draw(tag, attrs, text) {
    const e = document.createElementNS(svg.namespaceURI, tag);
    for (const [key, v] of Object.entries(attrs)) e.setAttribute(key, v);
    if (text !== undefined) e.textContent = text;
    svg.append(e);
  }
  for (const [l, r] of start > end
    ? [
        [0, end],
        [start, 360],
      ]
    : [[start, end]]) {
    draw("rect", {
      x: X(Math.max(from, l)),
      y: top,
      width: Math.max(0, X(Math.min(to, r)) - X(Math.max(from, l))),
      height: bottom - top,
      fill: "#194a46",
    });
  }
  for (let i = 0; i <= 4; i++) {
    const x = from + (i * (to - from)) / 4,
      y = lo + (i * (hi - lo)) / 4;
    draw("line", {
      x1: left,
      x2: right,
      y1: Y(y),
      y2: Y(y),
      stroke: "#38605b",
    });
    draw(
      "text",
      {
        x: left - 8,
        y: Y(y) + 4,
        fill: "#c6ded4",
        "text-anchor": "end",
        "font-size": 13,
      },
      Number(y.toPrecision(3)),
    );
    draw(
      "text",
      {
        x: X(x),
        y: bottom + 24,
        fill: "#c6ded4",
        "text-anchor": "middle",
        "font-size": 14,
      },
      Number(x.toFixed(timed ? 6 : 3)) + (timed ? " s" : "°"),
    );
  }
  // Large traces keep the exact numerical result and show extrema per screen column.
  function compress(arr, low) {
    const out = [];
    for (const p of arr) {
      const prev = out.at(-1);
      if (prev && Math.floor(X(prev[0])) === Math.floor(X(p[0]))) {
        if (low ? p[1] < prev[1] : p[1] > prev[1]) out[out.length - 1] = p;
      } else out.push(p);
    }
    return out;
  }
  const polygon = compress(lower, true).concat(
    compress(upper, false).reverse(),
  );
  draw("polygon", {
    points: polygon
      .map(([x, y]) => X(x).toFixed(3) + "," + Y(y).toFixed(3))
      .join(" "),
    fill: "#4e9b88",
    opacity: 0.7,
  });
  draw("line", {
    x1: left,
    x2: right,
    y1: Y(0),
    y2: Y(0),
    stroke: "#ffb58e",
    "stroke-width": 2,
    "stroke-dasharray": "6 4",
  });
  const stride = Math.max(1, Math.ceil(dots.length / 1000));
  dots.forEach(([x, y], i) => {
    if (i % stride === 0)
      draw("circle", { cx: X(x), cy: Y(y), r: 3.6, fill: "#ffe8a5" });
  });
  draw(
    "text",
    { x: left, y: 14, fill: "#daf1e7", "font-size": 13 },
    "Margin (" + spec.unit + ")",
  );
  $("graph").append(svg);
}
function loadCase(k) {
  if (k === "bench") {
    const example = BenchTrial.fixture(benchGroup, benchEveryFourth);
    activate(example.setup, example.csv, k);
    $("editor").hidden = true;
    return;
  }
  const e = trialExamples[k];
  activate(e.setup, e.csv, k);
  $("editor").hidden = true;
}
$("case").onchange = (e) => loadCase(e.target.value);
$("inputs").onclick = () => {
  $("editor").hidden = !$("editor").hidden;
};
for (const id of ["setupText", "traceText"]) $(id).oninput = () => dirty();
let importToken = { setupFile: 0, traceFile: 0 };
for (const [input, area] of [
  ["setupFile", "setupText"],
  ["traceFile", "traceText"],
])
  $(input).onchange = async (e) => {
    const f = e.target.files[0];
    if (!f) return;
    const token = ++importToken[input];
    const selected = active;
    dirty("Reading " + f.name + "…");
    try {
      if (f.size > 16000000) throw Error("File exceeds 16 MB.");
      const data = await f.text();
      if (token !== importToken[input] || active !== selected) return;
      $(area).value = data;
      dirty();
    } catch (err) {
      dirty(err.message);
    } finally {
      $(input).value = "";
    }
  };
$("run").onclick = () => {
  try {
    activate(JSON.parse($("setupText").value), $("traceText").value, "custom");
    $("editor").hidden = true;
  } catch (e) {
    dirty(e.message);
  }
};
$("reset").onclick = () => activate(active.setup, active.csv, active.selector);
$("cycle").onchange = renderCycle;
$("channel").onchange = renderGraph;
$("view").onchange = renderGraph;
$("saveSetup").onclick = () =>
  save("setup.json", "application/json", $("setupText").value);
$("saveTrace").onclick = () =>
  save("trace.csv", "text/csv", $("traceText").value);
$("export").onclick = () => {
  const { setup, csv, result } = active;
  save(
    result.setup + "-" + result.run + "-record.json",
    "application/json",
    JSON.stringify(
      {
        engine: Linkwright.VERSION,
        exportedUTC: new Date().toISOString(),
        setup,
        csv,
        result: Linkwright.summary(result),
        comparison: {
          baseline: $("baseline").value,
          observed: $("observed").value,
          notes: $("notes").value,
        },
      },
      null,
      2,
    ),
  );
};
$("reopen").onclick = () => $("recordFile").click();
$("recordFile").onchange = async (e) => {
  try {
    const file = e.target.files[0];
    if (!file) return;
    if (file.size > 48000000) throw Error("Record exceeds 48 MB.");
    const record = JSON.parse(await file.text());
    Linkwright.replay(record);
    activate(record.setup, record.csv, "custom");
    $("editor").hidden = true;
    for (const key of ["baseline", "observed", "notes"])
      if (typeof record.comparison?.[key] === "string")
        $(key).value = record.comparison[key];
    $("notice").textContent = "Saved calculations match the replay.";
  } catch (err) {
    dirty(err.message);
  } finally {
    $("recordFile").value = "";
  }
};
$("nativeFile").onchange = async (event) => {
  const file = event.target.files[0];
  if (!file) return;
  const selected = active;
  try {
    if (file.size > 16000000) throw Error("File exceeds 16 MB.");
    if ($("nativeLimit").value === "" || $("nativeError").value === "")
      throw Error(
        "Enter a limit and absolute sensor error before opening a recording.",
      );
    const values = [
      "nativeSignal",
      "nativeLimit",
      "nativeSense",
      "nativeError",
      "nativeRate",
    ].map((k) => $(k).value);
    const body = await file.text();
    if (active !== selected) return;
    const e = BenchTrial.fromNative(
      body,
      file.name,
      ...values.slice(0, 4),
      values[4] || null,
    );
    activate(e.setup, e.csv, "custom");
    $("nativeMessage").textContent = "";
    $("editor").hidden = true;
  } catch (err) {
    dirty(err.message);
    $("nativeMessage").textContent = err.message;
  } finally {
    $("nativeFile").value = "";
  }
};
const query = new URLSearchParams(location.search);
const benchGroup = Object.hasOwn(fieldData.bench, query.get("bench"))
  ? query.get("bench")
  : "empty/Sine_0";
const benchEveryFourth = query.get("samples") === "ends";
loadCase("bench");
