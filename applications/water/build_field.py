"""Rebuild the blinded-location BKTown summaries from the authors' pinned ZIP and EPANET 2.2."""

import argparse
import csv
import ctypes as ct
from datetime import datetime, timedelta
from hashlib import sha256
import io
import json
import math
from pathlib import Path
from statistics import median
from tempfile import TemporaryDirectory
from zipfile import ZipFile

ARCHIVE_SHA256 = "51474a77a35414dfdc0c2748450729fc4ea78fc2bf8c7b2f4cc3e56b6e3f7374"
PREFIX = "data/input/BKTown/HPC_run_08/"
RAW = "data/raw/BKTown/"
BASE = datetime(2020, 12, 8)


def stream_rows(z, path):
    with z.open(path) as f:
        yield from csv.DictReader(io.TextIOWrapper(f, encoding="utf-8-sig"))


class Snapshot:
    """Use the EPANET project API, with a separate solver initialization at each candidate."""

    def __init__(self, library, network, report):
        self.lib = ct.CDLL(str(library.resolve()))
        self.h = ct.c_void_p()
        self.warnings = set()
        self.check(self.lib.EN_createproject(ct.byref(self.h)))
        self.check(
            self.lib.EN_open(self.h, str(network).encode(), str(report).encode(), b"")
        )

    def check(self, code):
        if code >= 100:
            raise RuntimeError("EPANET error " + str(code))
        if code:
            self.warnings.add(code)

    def index(self, n):
        out = ct.c_int()
        self.check(self.lib.EN_getnodeindex(self.h, n.encode(), ct.byref(out)))
        return out.value

    def read(self, time, sensor_nodes, leak_node=None, flow=None):
        lib, h = self.lib, self.h
        self.check(
            lib.EN_settimeparam(h, 4, ct.c_long(int((time - BASE).total_seconds())))
        )
        if leak_node is not None:
            self.check(lib.EN_adddemand(h, leak_node, ct.c_double(flow), b"", b"crest"))
        self.check(lib.EN_openH(h))
        self.check(lib.EN_initH(h, 0))
        t = ct.c_long()
        self.check(lib.EN_runH(h, ct.byref(t)))
        result = []
        for i in sensor_nodes:
            v = ct.c_double()
            self.check(lib.EN_getnodevalue(h, i, 11, ct.byref(v)))
            result.append(v.value * 0.0980665)  # m of water head to bar
        self.check(lib.EN_closeH(h))
        if leak_node is not None:
            n = ct.c_int()
            self.check(lib.EN_getnumdemands(h, leak_node, ct.byref(n)))
            self.check(lib.EN_deletedemand(h, leak_node, n.value))
        return result

    def close(self):
        self.check(self.lib.EN_close(self.h))
        self.check(self.lib.EN_deleteproject(self.h))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("archive", type=Path)
    p.add_argument("epanet_library", type=Path)
    p.add_argument(
        "--output", type=Path, default=Path(__file__).with_name("field.json")
    )
    args = p.parse_args()
    h = sha256()
    with args.archive.open("rb") as stream:
        for b in iter(lambda: stream.read(1048576), b""):
            h.update(b)
    if h.hexdigest() != ARCHIVE_SHA256:
        raise SystemExit("Unexpected source archive")
    with ZipFile(args.archive) as z, TemporaryDirectory() as work:
        net = Path(work) / "field.inp"
        net.write_bytes(
            z.read(PREFIX + "BKTown_original_WDN_8Dec2020_17Dec2020_uncalibrated.inp")
        )
        sim = Snapshot(args.epanet_library, net, Path(work) / "hydraulics.rpt")
        events = []
        for row in stream_rows(z, RAW + "experiments/infield_logging_leaks_MC_3.csv"):
            a, b = [
                datetime.strptime(row["Datum"] + " " + row[k], "%d/%m/%Y %H:%M")
                for k in ("Beginuur", "Einduur")
            ]
            events.append(
                (
                    a,
                    b,
                    row["3GE FID"],
                    float(row["Opmerkingen"].split()[0])
                    * 3600
                    / (b - a).total_seconds(),
                )
            )
        sensor_ids = [
            r["asset_id"]
            for r in stream_rows(z, PREFIX + "pressure_loggers_asset_id.csv")
        ]
        present = {}
        for r in stream_rows(
            z,
            RAW + "pressure_meas/mobile_pressure_data-08122020_15022021-UTCp01h00.csv",
        ):
            if "2020-12-09" <= r["datetime"][:10] <= "2020-12-17":
                present[datetime.fromisoformat(r["datetime"])] = [
                    float(r[s]) if r[s] and math.isfinite(float(r[s])) else None
                    for s in sensor_ids
                ]

        def drop(start, i):
            def sample(offsets):
                values = [
                    present[start + timedelta(minutes=t)][i]
                    for t in offsets
                    if start + timedelta(minutes=t) in present
                ]
                return (
                    median(values)
                    if len(values) == len(offsets) and None not in values
                    else None
                )

            pre, active = sample(range(-8, -2)), sample(range(10, 30))
            return pre - active if pre is not None and active is not None else None

        runs = []
        omitted = []
        for ix, (start, end, site, flow) in enumerate(events):
            if any(
                other_start < start and other_end > start - timedelta(minutes=8)
                for other_start, other_end, _, _ in events
            ):
                omitted.append(
                    {
                        "date": start.isoformat(),
                        "reason": "baseline overlaps previous controlled discharge",
                    }
                )
                continue
            history = [
                [drop(start.replace(day=day), i) for day in range(9, 16)]
                for i in range(len(sensor_ids))
            ]
            measured = [drop(start, i) for i in range(len(sensor_ids))]
            keep = [
                i
                for i, hs in enumerate(history)
                if None not in hs and measured[i] is not None
            ]
            aliases = ["S%02d" % (i + 1) for i in keep]
            nodes = [sim.index(sensor_ids[i]) for i in keep]
            when = start + timedelta(minutes=20)
            ref = sim.read(when, nodes)
            if not all(0 < v < 16 for v in ref):
                raise ValueError("Out-of-range modeled sensor pressure")
            # Only compact deltas and anonymous logger labels are redistributed.
            predicted = {}
            for k, (_, _, possible, _) in enumerate(events):
                result = sim.read(when, nodes, sim.index(possible), flow)
                predicted["Z%02d" % (k + 1)] = [
                    round(a - b, 6) for a, b in zip(ref, result)
                ]
            runs.append(
                {
                    "date": start.isoformat(),
                    "true_zone": "Z%02d" % (ix + 1),
                    "discharge_m3_h": round(flow, 4),
                    "sensors": aliases,
                    "drop_bar": [round(measured[i], 6) for i in keep],
                    "background_differential_bar": [
                        [round(min(t), 6), round(max(t), 6)]
                        for t in (
                            [x - y for x, y in zip(history[i], history[keep[0]])]
                            for i in keep
                        )
                    ],
                    "hydraulic_drop_bar": predicted,
                }
            )
        sim.close()
        out = {
            "source_doi": "10.5281/zenodo.7255403",
            "archive_sha256": ARCHIVE_SHA256,
            "solver_warning_codes": sorted(sim.warnings),
            "excluded_events": omitted,
            "events": runs,
        }
        args.output.write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
