"""Embed the canonical CSV/JSON examples and package the offline browser trial."""
import json
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED

root = Path(__file__).resolve().parent
examples = {}
for key in ('opener', 'folder', 'reach'):
    examples[key] = {
        'setup': json.loads((root / 'examples' / (key + '.json')).read_text()),
        'csv': (root / 'examples' / (key + '.csv')).read_text(),
    }
(root / 'examples.js').write_bytes(('globalThis.trialExamples = ' +
    json.dumps(examples, indent=2) + ';\n').encode())
files = ['index.html', 'engine.js', 'examples.js', 'ui.js', 'style.css',
         'guide.html', 'test.cjs', 'cli.cjs', 'bench.js', 'field.html', 'field.js',
         'field.css', 'field-data.js', 'audit.js']
files += sorted('evidence/' + f.name for f in (root / 'evidence').iterdir() if f.is_file())
files += sorted('examples/' + f.name for f in (root / 'examples').iterdir()
                if f.suffix in {'.csv', '.json'})
output = root.parents[1] / 'downloads' / 'linkwright-cartoner-trial.zip'
with TemporaryDirectory(dir=output.parent) as staging:
    partial = Path(staging) / output.name
    with ZipFile(partial, 'w') as archive:
        for name in sorted(files):
            info = ZipInfo(name, (2026, 1, 1, 0, 0, 0))
            info.external_attr = 0o100644 << 16
            archive.writestr(info, (root / name).read_bytes(),
                             compress_type=ZIP_DEFLATED, compresslevel=9)
    partial.replace(output)
print(output.name, sha256(output.read_bytes()).hexdigest())
