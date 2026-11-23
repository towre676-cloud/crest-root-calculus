"""Verify the delivered Linkwright build. Optional checks use the named CAD/drawing/browser tools."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parent

def run(args, folder=ROOT):
    subprocess.run(args, cwd=folder, check=True)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--drawings',action='store_true',help='check delivered PDF/DXF/SVG with PyMuPDF and ezdxf')
    parser.add_argument('--cad',action='store_true',help='replay the 39-solid STEP checks with CadQuery 2.7')
    parser.add_argument('--browser',type=Path,metavar='CHROMIUM',help='test the offline UI with Playwright')
    args=parser.parse_args()
    if not __debug__:
        raise SystemExit('Run without -O so that numerical assertions execute.')
    for line in (ROOT/'SHA256SUMS.txt').read_text(encoding='utf8').splitlines():
        expected, name=line.split('  ',1)
        if sha256((ROOT/name).read_bytes()).hexdigest()!=expected:
            raise SystemExit('Checksum mismatch: '+name)
    print('Verified release checksums.', flush=True)
    run([sys.executable,'check.py','--check'],ROOT/'verification')
    run([sys.executable,'-m','unittest','discover','-v'],ROOT/'verification')
    node=shutil.which('node')
    if not node:
        raise SystemExit('Install Node.js to run the builder and inspection tests.')
    run([node,'test_geometry.cjs'])
    run([node,'inspection/test_bench.cjs'])
    run([node,'cad/test_motion.cjs'])
    if args.drawings:
        run([sys.executable,'kit/test_drawings.py'])
    if args.cad:
        run([sys.executable,'cad/test_assembly.py'])
        with TemporaryDirectory() as folder:
            run([sys.executable,'cad/build.py','--output',folder])
            rebuilt=json.loads((Path(folder)/'receipt.json').read_text(encoding='utf8'))
            archived=json.loads((ROOT/'cad/receipt.json').read_text(encoding='utf8'))
            if rebuilt!=archived:
                raise SystemExit('Regenerated CAD properties differ from the archived model.')
    if args.browser:
        exe=str(args.browser.resolve())
        run([sys.executable,'test_page.py','--chromium',exe])
        run([sys.executable,'inspection/test_workshop.py','--chromium',exe])
    print('All requested digital checks passed.')

if __name__=='__main__':
    main()
