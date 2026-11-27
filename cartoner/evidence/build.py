"""Build the published measurements with Python 3 and NumPy. See sources.md."""
import argparse
import csv
import hashlib
import json
import io
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from urllib.request import Request, urlopen
import numpy as np

HERE = Path(__file__).resolve().parent
MANIFEST = json.loads((HERE / 'manifest.json').read_text())


def verify(cache, download):
    for item in MANIFEST:
        path = cache / item['path']
        if download and not path.exists():
            path.parent.mkdir(parents=True, exist_ok=True)
            request = Request(item['url'], headers={'User-Agent': 'Linkwright/2.0 (data reproduction)'})
            with urlopen(request, timeout=120) as response:
                contents = response.read()
            if hashlib.sha256(contents).hexdigest() != item['sha256']:
                raise ValueError(f'Download differs from manifest: {path}')
            path.write_bytes(contents)
        if hashlib.sha256(path.read_bytes()).hexdigest() != item['sha256']:
            raise ValueError(f'Checksum mismatch: {path}')


def run_data(path):
    data = np.loadtxt(path, delimiter=',', comments='%', ndmin=2)
    if data.shape[1] != 6 or not np.isfinite(data).all():
        raise ValueError(path)
    t = np.rint(data[:, 0] * 1e6).astype('int64')
    i = np.rint(data[:, 2] * 1e6).astype('int64')
    if not np.allclose(i, data[:, 2] * 1e6, atol=1e-6, rtol=0):
        raise ValueError('Source current resolution has changed')
    if (np.diff(t) <= 0).any():
        raise ValueError('Timestamps must increase')
    summary = dict(file=path.name, rows=len(t), first=int(t[0]), last=int(t[-1]),
                   current=[int(i.min()), int(i.max())], flagged=int((abs(i)>10_000_000_000).sum()),
                   position=[float(data[:, 4].min()), float(data[:, 4].max())],
                   voltage=[float(data[:, 3].min()), float(data[:, 3].max())],
                   velocity=[float(data[:, 5].min()), float(data[:, 5].max())])
    return t, i, summary


def benchmark(paths):
    runs = [run_data(path) for path in paths]
    raw = compare(runs, False)
    screened = compare(runs, True)
    return dict(runs=[r[2] for r in runs], raw={k: raw[k] for k in ('outside','hidden','rate','inconsistent','brackets')}, **screened)


def compare(runs, screen):
    # The 10 A screen is specific to this servo corpus. No replacement or gap bridging.
    def admissible(i):
        return abs(i) <= 10_000_000_000 if screen else np.ones(len(i), dtype=bool)
    # Least integer mA/s at least as large as every measured training secant.
    rate = 0
    for t, i, _ in runs[:2]:
        dt = np.diff(t)
        eligible = admissible(i)[:-1] & admissible(i)[1:]
        rate = max(rate, int(((abs(np.diff(i)) + dt - 1) // dt)[eligible].max()))
    t, i, _ = runs[2]
    # Every fourth reading plus the last. Retain short final brackets too.
    keep = np.unique(np.append(np.arange(0, len(t), 4), len(t) - 1))
    hidden = np.setdiff1d(np.arange(len(t)), keep, assume_unique=True)
    right = np.searchsorted(keep, hidden)
    bad = np.concatenate(([0], np.cumsum(~admissible(i))))
    good_bracket = bad[keep[1:] + 1] == bad[keep[:-1]]
    eligible = good_bracket[right - 1]
    hidden, right = hidden[eligible], right[eligible]
    a, b = keep[right - 1], keep[right]
    left_allow, right_allow = rate * (t[hidden] - t[a]), rate * (t[b] - t[hidden])
    outside = (abs(i[hidden] - i[a]) > left_allow) | (abs(i[b] - i[hidden]) > right_allow)
    inconsistent = abs(np.diff(i[keep])) > rate * np.diff(t[keep])
    exceed = i[hidden] - np.maximum(i[a], i[b])
    k = int(exceed.argmax())
    left, right = int(a[k]), int(b[k])
    trace = [[int(t[j]), int(i[j])] for j in range(max(0,left-12), min(len(t), right+13))]
    return dict(rate=rate, hidden=len(hidden),
                outside=int(outside.sum()), brackets=int(good_bracket.sum()), inconsistent=int((inconsistent & good_bracket).sum()),
                overshoots=int((exceed > 0).sum()), peakLift=int(exceed[k]),
                window=[int(t[left]), int(t[right])], peakTime=int(t[hidden[k]]),
                example=trace, endpointMax=int(max(i[left],i[right])),
                peak=int(i[hidden[k]]))


def packaging(cache):
    piade = defaultdict(lambda: dict(intervals=0, duration_ms=Counter(), hoursRows=0))
    with (cache/'piade/raw_data.csv').open() as src:
        for r in csv.DictReader(src):
            group = piade[r['equipment_ID']]
            duration = int(r['elapsed'])
            if duration < 0 or abs(duration - (float(r['end'])-float(r['start']))*1000) > 2:
                raise ValueError('Unexpected elapsed unit or duration')
            group['intervals'] += 1
            group['duration_ms'][r['type']] += duration
    with (cache/'piade/sequences_1h_data.csv').open() as src:
        for r in csv.DictReader(src):
            piade[r['equipment_ID']]['hoursRows'] += 1
    comopi = defaultdict(lambda: dict(sensorRows=0, alarmRows=0, warningRows=0, targetRows=0, alarms53=0, alarms54=0))
    for name, key in [('sensors','sensorRows'), ('alarm','alarmRows'), ('warnings','warningRows')]:
        with (cache/f'comopi/industrial_dataset_{name}_10m_agg.csv').open() as src:
            reader = csv.reader(src)
            cols = next(reader)
            serial = cols.index('_serial')
            if name == 'alarm':
                idx53, idx54 = cols.index('AL_53'), cols.index('AL_54')
            for row in reader:
                group = comopi[row[serial]]
                group[key] += 1
                if name == 'alarm':
                    x, y = int(row[idx53]), int(row[idx54])
                    group['alarms53'] += x
                    group['alarms54'] += y
                    group['targetRows'] += int(x+y > 0)
    counts, first, last = defaultdict(Counter), {}, {}
    with zipfile.ZipFile(cache/'alpi/data.zip') as archive:
        with archive.open('alarms_log_data/raw/alarms.csv') as body:
            for r in csv.DictReader(io.TextIOWrapper(body, encoding='utf-8')):
                serial, alarm, stamp = r['serial'], r['alarm'], r['timestamp']
                counts[serial][alarm] += 1
                first[serial] = min(first.get(serial, stamp), stamp)
                last[serial] = max(last.get(serial, stamp), stamp)
    alpi = {}
    for serial, counter in sorted(counts.items()):
        code, freq = counter.most_common(1)[0]
        alpi[serial] = dict(rows=sum(counter.values()), distinct=len(counter), top=code,
                            topCount=freq, first=first[serial], last=last[serial])
    return dict(piade=dict(sorted(piade.items())), comopi=dict(sorted(comopi.items())), alpi=alpi,
                alpiCodes=len(set().union(*counts.values())))


def main():
    arg = argparse.ArgumentParser(description=__doc__)
    arg.add_argument('cache', type=Path)
    arg.add_argument('--download', action='store_true')
    arg.add_argument('--output', type=Path, default=HERE.parent/'field-data.js')
    args = arg.parse_args()
    verify(args.cache, args.download)
    groups = defaultdict(list)
    for item in MANIFEST:
        path = args.cache/item['path']
        if path.suffix == '.txt':
            group = path.parent.name.removeprefix('output_') + '/' + path.stem.split('_output_')[0]
            groups[group].append(path)
    urls = {str(args.cache / f["path"]):f["url"] for f in MANIFEST}
    benches = {}
    for group, paths in sorted(groups.items()):
        paths = sorted(paths)
        if [p.stem[-2:] for p in paths] != ['00','01','02']:
            raise ValueError('Incomplete repeated experiment')
        benches[group] = benchmark(paths)
        for run, path in zip(benches[group]["runs"], paths):
            run["url"] = urls[str(path)]
    result = dict(bench=benches, **packaging(args.cache))
    args.output.write_bytes(('globalThis.fieldData = '+json.dumps(result, separators=(',',':'),sort_keys=True)+';\n').encode())
    print(len(benches), 'matched groups;', sum(r['rows'] for b in benches.values() for r in b['runs']), 'servo measurements')


if __name__ == '__main__':
    main()
