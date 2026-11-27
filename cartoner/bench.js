/* Physical four-bar replay and the source-native DaRUS import adapter. */
(function (root) {
  "use strict";
  const signals = {
    current: [2, "Motor current", "mA"],
    command: [1, "Command current", "mA"],
    voltage: [3, "Supply voltage", "V"],
    position: [4, "Crank position", "rad"],
    velocity: [5, "Crank velocity", "rad/s"],
  };
  function geometry() {
    return [
      {
        name: "Stuttgart four-bar · nominal pivot dimensions",
        groundDeg: 0,
        ground: [95, 95],
        crank: [30, 30],
        coupler: [104, 104],
        rocker: [60, 60],
      },
    ];
  }
  const header = "setup,run,cycle,time_s,channel,value";
  function setup(name, id, channel, window, limit, sense, error, rate) {
    const [, label, unit] = signals[channel];
    return {
      schema: 2,
      axis: "time_s",
      setup: id,
      name,
      source: "machine",
      channels: [
        { id: channel, label, unit, window, limit, sense, error, rate },
      ],
      loops: [],
      events: [],
      handoffs: [],
    };
  }
  function fixture(key, endpoints) {
    const b = fieldData.bench[key];
    if (!b) throw Error("Unknown published experiment.");
    const id = "darus-" + key.replace("/", "-");
    const limit = Math.round((b.endpointMax + b.peak) / 2) / 1e6;
    const p = setup(
      "Measured current · " + key,
      id,
      "current",
      b.window.map((x) => x / 1e6),
      limit,
      "max",
      0,
      null,
    );
    p.source = "published-bench";
    p.loops = geometry();
    p.provenance = {
      doi: "10.18419/darus-4152",
      file: b.runs[2].file,
      condition: key.split("/")[0],
      selection: endpoints
        ? "Two retained readings"
        : "Original readings in and around the selected window",
      basis:
        "Illustrative midpoint threshold selected after seeing the readings. Zero error refers to recorded values. Physical sensor error and rate are uncalibrated.",
    };
    const points = endpoints
      ? b.example.filter(([t]) => b.window.includes(t))
      : b.example;
    return {
      setup: p,
      csv:
        [
          header,
          ...points.map(([t, i]) =>
            [
              id,
              "repeat-02",
              "excerpt",
              (t / 1e6).toFixed(6),
              "current",
              (i / 1e6).toFixed(6),
            ].join(","),
          ),
        ].join("\n") + "\n",
    };
  }
  function fromNative(text, file, channel, limit, sense, error, rate) {
    if (!Object.hasOwn(signals, channel))
      throw Error("Choose a recorded signal.");
    const lines = text
      .trim()
      .replace(/^\uFEFF/, "")
      .split(/\r?\n/);
    const columns =
      "% Time [s], Goal Current [mA], Pres. Current [mA], Pres. Voltage [V], Pres. Position [rad], Pres. Velocity [rad/s]";
    if (lines.shift()?.trim() !== columns)
      throw Error("Expected the original six-column DaRUS servo header.");
    if (lines.length < 2 || lines.length > 100000)
      throw Error("Supply 2 to 100,000 native readings.");
    const id = "native-servo",
      col = signals[channel][0];
    const output = [header];
    let first, last;
    for (let j = 0; j < lines.length; j++) {
      const row = lines[j].split(",").map((s) => s.trim());
      if (row.length !== 6 || row.some((s) => !/^-?\d+(\.\d{1,6})?$/.test(s)))
        throw Error("Invalid native reading at line " + (j + 2));
      const t = Number(row[0]);
      if (t < 0 || (last !== undefined && t <= last))
        throw Error("Native timestamps must increase.");
      if (j === 0) first = row[0];
      last = t;
      output.push(
        [id, "native", "recording", row[0], channel, row[col]].join(","),
      );
    }
    const p = setup(
      file.slice(0, 110) || "Native servo recording",
      id,
      channel,
      [first, last.toFixed(6)],
      limit,
      sense,
      error,
      rate,
    );
    p.provenance = { file, format: "DaRUS six-column servo export" };
    return { setup: p, csv: output.join("\n") + "\n" };
  }
  root.BenchTrial = { fixture, fromNative };
})(globalThis);
