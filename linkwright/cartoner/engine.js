/* Linkwright | Exact comparisons for planar closure, phase traces and event order. */
(function (root) {
  "use strict";
  const VERSION = "2.0.0";
  const S = 1000000n;
  const min = (a, b) => (a < b ? a : b);
  const max = (a, b) => (a > b ? a : b);
  const abs = (a) => (a < 0n ? -a : a);
  function requireThat(condition, message) {
    if (!condition) throw new Error(message);
  }
  function decimal(
    input,
    name,
    lower = -100000000n * S,
    upper = 100000000n * S,
  ) {
    const raw = String(input);
    requireThat(
      /^-?\d{1,9}(\.\d{1,6})?$/.test(raw),
      name + " needs a decimal with at most six places.",
    );
    const [whole, fraction = ""] = raw.replace("-", "").split(".");
    const v =
      (BigInt(whole) * S + BigInt(fraction.padEnd(6, "0"))) *
      (raw[0] === "-" ? -1n : 1n);
    requireThat(
      v >= lower && v <= upper,
      name + " is outside its allowed range.",
    );
    return v;
  }
  // All returned decimal margins are exact. Phase of an interior witness is also returned as a ratio.
  function format(v, places = 6) {
    const n = abs(v)
      .toString()
      .padStart(places + 1, "0");
    const part = n.slice(-places).replace(/0+$/, "");
    return (v < 0n ? "-" : "") + n.slice(0, -places) + (part ? "." + part : "");
  }
  const twiceMargin = (v) => format(v * 5n, 13);
  const angle = (x, label) => decimal(x, label, 0n, 360n * S);
  function name(x, label) {
    requireThat(
      typeof x === "string" &&
        x.trim() === x &&
        x.length > 0 &&
        x.length <= 120 &&
        !/[\x00-\x1f]/.test(x),
      "Invalid " + label,
    );
    return x;
  }
  function id(x, label) {
    requireThat(
      typeof x === "string" && /^[a-zA-Z0-9][a-zA-Z0-9_.-]{0,79}$/.test(x),
      "Invalid " + label,
    );
    return x;
  }
  function bounds(x, label) {
    requireThat(
      Array.isArray(x) && x.length === 2,
      label + " needs [low, high] in mm.",
    );
    const a = decimal(x[0], label, 1n),
      b = decimal(x[1], label, a);
    return [a, b];
  }
  // For a restricted sweep d²=g²+a²−2ga cos(theta). In a positive length
  // rectangle its maximum is at a corner. Its minimum is at a corner or
  // the orthogonal projection onto an edge (including the g=a diagonal).
  function restrictedLoop(spec, lengths) {
    const [[g, G], [a, A], [b, B], [c, C]] = lengths;
    requireThat(
      Array.isArray(spec.strokeCos) && spec.strokeCos.length === 2,
      "strokeCos needs ordered cosine bounds.",
    );
    const left = decimal(spec.strokeCos[0], "Lower cosine", -S, S);
    const right = decimal(spec.strokeCos[1], "Upper cosine", left, S);
    const clamp = (x, l, h) => min(h, max(l, x));
    function extremes(cosine) {
      const points = [];
      for (const x of [g, G]) {
        for (const y of [a, A]) points.push([x * S, y * S]);
        points.push([x * S, clamp(x * cosine, a * S, A * S)]);
      }
      for (const y of [a, A])
        points.push([clamp(y * cosine, g * S, G * S), y * S]);
      return points.map(([x, y]) => ({
        v: (x * x + y * y) * S - 2n * x * y * cosine,
        x,
        y,
      }));
    }
    const closest = extremes(right).reduce((x, y) => (x.v <= y.v ? x : y));
    const furthest = extremes(left).reduce((x, y) => (x.v >= y.v ? x : y));
    const farBC = abs(B - c) >= abs(b - C) ? [B, c] : [b, C];
    const squaredScale = S * S * S;
    const folded = closest.v - (farBC[0] - farBC[1]) ** 2n * squaredScale;
    const extended = (b + c) ** 2n * squaredScale - furthest.v;
    const inner = folded <= extended,
      w = inner ? closest : furthest;
    const worst = min(folded, extended);
    return {
      kind: "geometry",
      label: spec.name,
      status: worst > 0n ? "meets" : "violation",
      unit: "mm²",
      margin: format(worst, 30),
      folded: format(folded, 30),
      extended: format(extended, 30),
      constraint: inner ? "folded" : "extended",
      cosine: format(inner ? right : left),
      witness: {
        ground: format(w.x, 12),
        crank: format(w.y, 12),
        coupler: format(inner ? farBC[0] : b),
        rocker: format(inner ? farBC[1] : c),
      },
    };
  }
  function loop(spec) {
    const names = ["ground", "crank", "coupler", "rocker"];
    name(spec.name, "loop name");
    const lengths = names.map((k) => bounds(spec[k], k));
    if (spec.strokeCos !== undefined) return restrictedLoop(spec, lengths);
    const [[g, G], [a, A], [b, B], [c, C]] = lengths;
    let nearG, nearA;
    if (g > A) {
      nearG = g;
      nearA = A;
    } else if (a > G) {
      nearG = G;
      nearA = a;
    } else {
      nearG = nearA = max(g, a);
    }
    const farBC = abs(B - c) >= abs(b - C) ? [B, c] : [b, C];
    const folded = abs(nearG - nearA) - abs(farBC[0] - farBC[1]);
    const extended = b + c - G - A;
    const limiting = folded <= extended;
    const worst = min(folded, extended);
    const witness = limiting ? [nearG, nearA, ...farBC] : [G, A, b, c];
    const datum = angle(spec.groundDeg, "Ground direction");
    return {
      kind: "geometry",
      label: name(spec.name, "loop name"),
      status: worst > 0n ? "meets" : "violation",
      margin: format(worst),
      unit: "mm",
      folded: format(folded),
      extended: format(extended),
      constraint: limiting ? "folded" : "extended",
      phase: format((datum + (limiting ? 0n : 180n * S)) % (360n * S)),
      witness: Object.fromEntries(names.map((k, i) => [k, format(witness[i])])),
    };
  }
  function project(p) {
    requireThat(
      p && (p.schema === 1 || p.schema === 2),
      "Use setup schema 1 or 2.",
    );
    const time = p.schema === 2 && p.axis === "time_s";
    requireThat(
      p.schema === 1 ? p.axis === undefined : time || p.axis === "phase_deg",
      "Schema 2 requires axis time_s or phase_deg. Use schema 2 to specify an axis.",
    );
    const coordinate = time ? (x, key) => decimal(x, key, 0n) : angle;
    id(p.setup, "setup ID");
    name(p.name, "setup name");
    requireThat(
      ["synthetic", "machine"].includes(p.source) ||
        (p.schema === 2 && p.source === "published-bench"),
      "Use machine, synthetic, or published-bench (schema 2) for source.",
    );
    for (const key of ["loops", "channels", "events", "handoffs"])
      requireThat(
        Array.isArray(p[key]) && p[key].length <= 32,
        "Invalid " + key + " list (up to 32).",
      );
    requireThat(
      p.channels.length + p.handoffs.length > 0,
      "Declare at least one trace requirement.",
    );
    const seen = new Set();
    function register(k) {
      id(k, "channel ID");
      requireThat(!seen.has(k), "Duplicate ID " + k);
      seen.add(k);
    }
    const channels = p.channels.map((c) => {
      register(c.id);
      name(c.label, "channel name");
      name(c.unit, "unit");
      requireThat(
        ["min", "max"].includes(c.sense),
        "Channel sense must be min or max.",
      );
      requireThat(
        Array.isArray(c.window) && c.window.length === 2,
        "Set both phase window endpoints.",
      );
      const w = c.window.map((x) => coordinate(x, "Window"));
      if (time) requireThat(w[1] > w[0], "Elapsed-time windows must increase.");
      requireThat(w[0] !== w[1], "Use [0,360] for a complete cycle.");
      return {
        ...c,
        time,
        w,
        threshold: decimal(c.limit, "Limit"),
        e: decimal(c.error, "Sensor error", 0n),
        rate:
          c.rate === null
            ? null
            : decimal(c.rate, "Maximum rate per axis unit", 0n),
      };
    });
    const events = new Map(
      p.events.map((e) => {
        register(e.id);
        return [e.id, coordinate(time ? e.errorS : e.errorDeg, "Event error")];
      }),
    );
    const handoffs = p.handoffs.map((h) => {
      name(h.label, "event check");
      requireThat(
        events.has(h.first) && events.has(h.second) && h.first !== h.second,
        "Unknown or identical handoff events.",
      );
      return {
        ...h,
        time,
        dwell: coordinate(
          time ? h.minS : h.minDeg,
          "Required event separation",
        ),
      };
    });
    return {
      loops: p.loops.map(loop),
      channels,
      events,
      handoffs,
      seen,
      coordinate,
      time,
    };
  }
  function csvRows(data, axis) {
    requireThat(
      typeof data === "string" && data.length < 16000000,
      "Trace file exceeds 16 MB.",
    );
    let rows = [],
      row = [],
      cell = "",
      quoted = false,
      closed = false;
    for (const ch of data.replace(/^\uFEFF/, "").replace(/\r\n?/g, "\n") +
      "\n") {
      if (quoted) {
        if (ch === '"') {
          quoted = false;
          closed = true;
        } else cell += ch;
      } else if (ch === '"' && (!cell || closed)) {
        if (closed) cell += '"';
        quoted = true;
        closed = false;
      } else if (ch === "," || ch === "\n") {
        row.push(cell);
        cell = "";
        closed = false;
        if (ch === "\n") {
          if (row.length > 1 || row[0]) rows.push(row);
          row = [];
        }
      } else {
        requireThat(!closed && ch !== '"', "Malformed quoted CSV field.");
        cell += ch;
      }
    }
    requireThat(!quoted, "Unclosed CSV quote.");
    const header = rows.shift() || [];
    requireThat(
      header.slice(0, 4).join(",") === "setup,run,cycle," + axis,
      "Expected setup,run,cycle," + axis + " as the first four columns.",
    );
    requireThat(
      rows.length > 0 && rows.length <= 100000,
      "Supply between 1 and 100,000 trace rows.",
    );
    return { header, rows };
  }
  function readTrace(p, spec, csv) {
    const cycles = new Map();
    let run = null;
    const { header, rows } = csvRows(csv, spec.time ? "time_s" : "phase_deg");
    const long = header.slice(4).join(",") === "channel,value";
    const keys = header.slice(4);
    if (!long)
      requireThat(
        keys.length > 0 &&
          new Set(keys).size === keys.length &&
          keys.every((k) => spec.seen.has(k)),
        "Wide CSV contains an unknown or duplicate channel.",
      );
    let count = 0;
    for (const [i, raw] of rows.entries()) {
      try {
        requireThat(
          raw.length === header.length,
          "CSV row length differs from header.",
        );
        const readings = long
          ? [raw]
          : keys.flatMap((k, j) => {
              const v = raw[j + 4];
              if (!v.trim()) return [];
              if (spec.events.has(k) && v !== "1")
                throw Error("Use 1 once per event, blank elsewhere.");
              return [[...raw.slice(0, 4), k, spec.events.has(k) ? "" : v]];
            });
        requireThat(readings.length > 0, "Empty measurement row.");
        for (const row of readings) {
          requireThat(
            ++count <= 100000,
            "One file supports up to 100,000 observations.",
          );
          const [setup, r, cycle, phase, ch, value] = row;
          requireThat(
            setup === p.setup,
            "Setup ID differs from selected revision.",
          );
          id(r, "run ID");
          id(cycle, "cycle ID");
          if (run === null) run = r;
          requireThat(run === r, "Import one run per file.");
          requireThat(spec.seen.has(ch), "Unknown channel " + ch);
          const x = spec.coordinate(
            phase,
            spec.time ? "Elapsed time" : "Encoder phase",
          );
          if (!cycles.has(cycle)) cycles.set(cycle, new Map());
          const signals = cycles.get(cycle);
          if (!signals.has(ch)) signals.set(ch, []);
          const event = spec.events.has(ch);
          requireThat(
            !event || value === "",
            "Timestamp events require an empty value.",
          );
          signals
            .get(ch)
            .push({ x, v: event ? null : decimal(value, "Reading") });
        }
      } catch (e) {
        throw new Error("CSV row " + (i + 2) + ": " + e.message);
      }
    }
    requireThat(cycles.size <= 2000, "One import supports up to 2,000 cycles.");
    for (const [cycle, signals] of cycles)
      for (const [ch, points] of signals) {
        points.sort((a, b) => (a.x < b.x ? -1 : a.x > b.x ? 1 : 0));
        requireThat(
          !spec.events.has(ch) || points.length === 1,
          "Repeated timestamp " + ch + " in cycle " + cycle,
        );
        requireThat(
          points.every((p, i) => i === 0 || p.x !== points[i - 1].x),
          "Duplicate phase for " + ch + " in cycle " + cycle,
        );
      }
    return { cycles, run };
  }
  function channel(c, points) {
    const result = {
      kind: "trace",
      id: c.id,
      label: c.label,
      unit: c.unit,
      window: c.window,
      samples: points.length,
    };
    const windows =
      c.w[0] < c.w[1]
        ? [c.w]
        : [
            [0n, c.w[1]],
            [c.w[0], 360n * S],
          ];
    if (points.length === 0) return { ...result, status: "missing" };
    const sign = c.sense === "min" ? 1n : -1n;
    const q = points.map(({ x, v }) => {
      const y = (v - c.threshold) * sign * S;
      return { x, y, lo: y - c.e * S, hi: y + c.e * S };
    });
    if (c.rate === null) {
      const inside = q.filter((p) =>
        windows.some(([a, b]) => a <= p.x && p.x <= b),
      );
      const best = inside.reduce((v, p) => (!v || p.hi < v.hi ? p : v), null);
      return {
        ...result,
        status:
          best && best.hi < 0n ? "violation" : best ? "unresolved" : "missing",
        lower: null,
        upper: best ? format(best.hi, 12) : null,
        sampled: best ? format(best.y, 12) : null,
        phase: best ? format(best.x) : null,
        unbounded: true,
        points: q.map((p) => [
          Number(p.x) / 1e6,
          Number(p.y) / 1e12,
          Number(p.lo) / 1e12,
          Number(p.hi) / 1e12,
        ]),
      };
    }
    // Tighten each sensor interval against all other observations in linear time.
    function tighten(i, j) {
      const a = q[i],
        b = q[j],
        step = c.rate * abs(a.x - b.x);
      a.lo = max(a.lo, b.lo - step);
      a.hi = min(a.hi, b.hi + step);
    }
    for (let i = 1; i < q.length; i++) tighten(i, i - 1);
    for (let i = q.length - 2; i >= 0; i--) tighten(i, i + 1);
    const invalid = q.find((p) => p.lo > p.hi);
    if (invalid)
      return { ...result, status: "inconsistent", phase: format(invalid.x) };
    let low = null,
      high = null,
      witness = null,
      upperAt = null,
      sampleMin = null;
    for (const [start, end] of windows) {
      for (const p of q)
        if (p.x >= start && p.x <= end) {
          sampleMin = sampleMin === null ? p.y : min(sampleMin, p.y);
          if (high === null || 2n * p.hi < high) {
            high = 2n * p.hi;
            upperAt = p.x;
          }
        }
      for (let i = 1; i < q.length; i++) {
        const a = q[i - 1],
          b = q[i],
          l = max(a.x, start),
          r = min(b.x, end);
        if (l > r) continue;
        const k = c.rate;
        const lower = max(
          max(2n * (a.lo - k * (r - a.x)), 2n * (b.lo - k * (b.x - l))),
          a.lo + b.lo - k * (b.x - a.x),
        );
        if (low === null || lower < low) {
          low = lower;
          const den = k === 0n ? 1n : 2n * k;
          const num =
            k === 0n
              ? l
              : min(r * den, max(l * den, a.lo - b.lo + k * (a.x + b.x)));
          witness = {
            numerator: String(num),
            denominator: String(den * S),
            [c.time ? "seconds" : "degrees"]: Number(num) / Number(den * S),
          };
        }
        for (const t of [l, r]) {
          const top = 2n * min(a.hi + k * (t - a.x), b.hi + k * (b.x - t));
          if (high === null || top < high) {
            high = top;
            upperAt = t;
          }
        }
      }
    }
    const covered = windows.every(([l, r]) => q[0].x <= l && q.at(-1).x >= r);
    const status =
      high !== null && high < 0n
        ? "violation"
        : !covered || low === null
          ? "missing"
          : low >= 0n
            ? "meets"
            : "unresolved";
    return {
      ...result,
      status,
      covered,
      lower: covered && low !== null ? twiceMargin(low) : null,
      upper: high === null ? null : twiceMargin(high),
      sampled: sampleMin === null ? null : format(sampleMin, 12),
      phase:
        status === "violation"
          ? format(upperAt)
          : witness
            ? String(c.time ? witness.seconds : witness.degrees)
            : null,
      lowerWitness: witness,
      points: q.map((p) => [
        Number(p.x) / 1e6,
        Number(p.y) / 1e12,
        Number(p.lo) / 1e12,
        Number(p.hi) / 1e12,
      ]),
    };
  }
  function handoff(h, events, data) {
    const first = data.get(h.first)?.[0],
      second = data.get(h.second)?.[0];
    const base = {
      kind: "handoff",
      label: h.label,
      unit: h.time ? "s" : "deg",
      first: h.first,
      second: h.second,
    };
    if (!first || !second) return { ...base, status: "missing" };
    const margin = second.x - first.x - h.dwell;
    const error = events.get(h.first) + events.get(h.second);
    return {
      ...base,
      status:
        margin - error >= 0n
          ? "meets"
          : margin + error < 0n
            ? "violation"
            : "unresolved",
      lower: format(margin - error),
      upper: format(margin + error),
      phase: format(second.x),
    };
  }
  function analyze(p, csv) {
    const spec = project(p),
      { cycles, run } = readTrace(p, spec, csv);
    const records = Array.from(cycles, ([cycle, data]) => ({
      cycle,
      checks: [
        ...spec.channels.map((c) => channel(c, data.get(c.id) || [])),
        ...spec.handoffs.map((h) => handoff(h, spec.events, data)),
      ],
    }));
    const statuses = new Set(
      [...spec.loops, ...records.flatMap((r) => r.checks)].map((x) => x.status),
    );
    const status = [
      "violation",
      "inconsistent",
      "missing",
      "unresolved",
      "meets",
    ].find((s) => statuses.has(s));
    if (spec.time) {
      for (const cycle of records)
        for (const check of cycle.checks) {
          if (Object.hasOwn(check, "phase")) {
            check.time_s = check.phase;
            delete check.phase;
          }
        }
    }
    return {
      ...(p.schema === 2 ? { axis: p.axis } : {}),
      version: VERSION,
      setup: p.setup,
      source: p.source,
      run,
      status,
      loops: spec.loops,
      cycles: records,
    };
  }
  function summary(r) {
    return {
      ...r,
      cycles: r.cycles.map((c) => ({
        ...c,
        checks: c.checks.map(({ points, ...x }) => x),
      })),
    };
  }
  function replay(record) {
    requireThat(
      record.engine === VERSION ||
        (record.engine === "1.0.0" &&
          record.setup?.schema === 1 &&
          !("axis" in record.setup)),
      "This record needs engine " + record.engine,
    );
    const fresh = analyze(record.setup, record.csv);
    function stable(x) {
      if (Array.isArray(x)) return x.map(stable);
      if (x && typeof x === "object")
        return Object.fromEntries(
          Object.keys(x)
            .sort()
            .map((k) => [k, stable(x[k])]),
        );
      return x;
    }
    requireThat(
      JSON.stringify(stable({ ...summary(fresh), version: record.engine })) ===
        JSON.stringify(stable(record.result)),
      "Stored result disagrees with replay.",
    );
    return fresh;
  }
  const api = { VERSION, analyze, project, summary, replay };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.Linkwright = api;
})(globalThis);
