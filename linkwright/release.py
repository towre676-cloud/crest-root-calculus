"""Build the offline ZIP from the named release files, omitting local measurements and caches."""
import argparse
from hashlib import sha256
from html import escape
from pathlib import Path
import re
import subprocess
from tempfile import TemporaryDirectory
from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED

ROOT=Path(__file__).resolve().parent
GUIDES={
    'kit/README.md':('kit/assembly-guide.html','Fabrication and assembly'),
    'verification/README.md':('verification/engineering.html','Engineering verification'),
    'custom-design.md':('custom-guide.html','Design another motion')
}
def guides():
    targets={ROOT/source:ROOT/output for source,(output,_) in GUIDES.items()}
    targets[ROOT/'README.md']=ROOT/'index.html'
    for source,(output,title) in GUIDES.items():
        out=ROOT/output; source=ROOT/source
        body=subprocess.check_output(['pandoc',str(source),'-f','gfm+tex_math_dollars','-t','html5','--mathml','--wrap=none'],text=True,encoding='utf8')
        def link(m):
            href=m.group(1); path,sep,anchor=href.partition('#')
            if '://' in path or not path:return m.group(0)
            target=(out.parent/path).resolve()
            if target.is_dir():target=target/'README.md'
            if target not in targets:return m.group(0)
            import os
            new=os.path.relpath(targets[target],out.parent).replace(os.sep,'/')
            return 'href="'+new+(sep+anchor if sep else '')+'"'
        body=re.sub(r'href="([^"]+)"',link,body)
        prefix='' if out.parent==ROOT else '../'
        out.write_text(f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Linkwright · {escape(title)}</title><link rel="stylesheet" href="{prefix}inspection/bench.css"><style>main{{max-width:960px}} h2,h3{{margin-top:1.4em}} img{{max-width:100%;height:auto}} table{{border-collapse:collapse;width:100%;font-size:.94rem}} th,td{{text-align:left;vertical-align:top;border-bottom:1px solid #b8c9c6;padding:.5em .6em;overflow-wrap:anywhere}}th{{background:#edf2ef}} li{{margin:.45em 0}} pre{{overflow:auto;padding:1em;background:#eef2f2}} math{{overflow:auto}} @media(max-width:600px){{table{{display:block;overflow:auto}}}}</style></head><body><header><div><a href="{prefix}index.html">Linkwright</a><span>LW-159-R1</span></div></header><main><nav><a href="{prefix}index.html">Workshop</a><a href="{prefix}assembly-view.html">3D assembly</a><a href="{prefix}commissioning.html">Physical inspection</a></nav>{body}</main></body></html>''',encoding='utf8')
def package(output):
    names=sorted((ROOT/'release-files.txt').read_text(encoding='utf8').splitlines())
    if len(names)!=len(set(names)) or 'release-files.txt' not in names:raise ValueError('Invalid file list')
    for name in names:
        if not name or name.startswith('/') or '..' in Path(name).parts:raise ValueError(name)
    manifest=''.join(sha256((ROOT/n).read_bytes()).hexdigest()+'  '+n+'\n' for n in names)
    (ROOT/'SHA256SUMS.txt').write_bytes(manifest.encode())
    output.parent.mkdir(parents=True,exist_ok=True)
    with TemporaryDirectory(dir=output.parent) as tmp:
        archive=Path(tmp)/'release.zip'
        with ZipFile(archive,'w') as z:
            for n in sorted(names+['SHA256SUMS.txt']):
                info=ZipInfo(n,(2026,11,23,0,0,0))
                info.external_attr=0o100644<<16
                z.writestr(info,(ROOT/n).read_bytes(),compress_type=ZIP_DEFLATED,compresslevel=9)
        archive.replace(output)
    print(output.name, sha256(output.read_bytes()).hexdigest(), len(names)+1, 'files')
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--guides',action='store_true',help='regenerate printable HTML from Markdown using pandoc')
    parser.add_argument('--output',type=Path,help='write the offline ZIP and refresh SHA256SUMS.txt')
    args=parser.parse_args()
    if not (args.guides or args.output):parser.error('Use --guides, --output, or both.')
    if args.guides:guides()
    if args.output:package(args.output.resolve())
