"""Open seven Crest–Root applications on this computer using Python 3.11+."""
import argparse
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import os
import mimetypes
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
from urllib.parse import urlsplit
import webbrowser

ROOT = Path(__file__).resolve().parent
CATALOGUE = json.loads((ROOT/'workshop/catalogue.json').read_text(encoding='utf-8'))
CAPACITY = threading.BoundedSemaphore(2)
MAX_INPUT = 512 * 1024

def execute(app, payload, timeout=60, raw=None):
    """Run the same file-based CLI used from a terminal, outside the HTTP process."""
    with tempfile.TemporaryDirectory() as directory:
        case = Path(directory)/'case.json'
        if raw is None:
            case.write_text(json.dumps(payload, ensure_ascii=False, allow_nan=False), encoding='utf-8')
        else:
            case.write_bytes(raw)
        result = subprocess.run([sys.executable,'-B',str(ROOT/app/'app.py'),str(case)],
                                cwd=ROOT/app, capture_output=True, text=True, encoding='utf-8',
                                env={**os.environ,'PYTHONIOENCODING':'utf-8'},timeout=timeout)
    if result.returncode:
        lines = result.stderr.strip().splitlines() or result.stdout.strip().splitlines()
        raise ValueError(lines[-1][:500] if lines else 'The calculation could not complete.')
    return json.loads(result.stdout)

class Handler(BaseHTTPRequestHandler):
    def send_bytes(self, status, data, kind):
        self.send_response(status)
        self.send_header('Content-Type',kind)
        self.send_header('Content-Length',str(len(data)))
        self.send_header('Cache-Control','no-store')
        self.send_header('X-Content-Type-Options','nosniff')
        self.send_header('Referrer-Policy','no-referrer')
        self.end_headers()
        try:
            self.wfile.write(data)
        except (BrokenPipeError, ConnectionResetError):
            pass

    def json(self, status, data):
        self.send_bytes(status,json.dumps(data,allow_nan=False).encode(),'application/json; charset=utf-8')

    def local_host(self):
        port=self.server.server_address[1]
        return self.headers.get('Host') in {f'127.0.0.1:{port}',f'localhost:{port}'}

    def do_GET(self):
        if not self.local_host():
            return self.json(403,{'error':'Use the local address printed by the workshop.'})
        path=urlsplit(self.path).path
        static={'/':('workshop/index.html','text/html'),'/ui/style.css':('workshop/style.css','text/css'),
                '/ui/workshop.js':('workshop/workshop.js','text/javascript'),
                '/api/catalogue':('workshop/catalogue.json','application/json')}
        if path in static:
            filename,kind=static[path]
            return self.send_bytes(200,(ROOT/filename).read_bytes(),kind+'; charset=utf-8')
        parts=path.strip('/').split('/')
        if len(parts)==2 and parts[0]=='guide' and parts[1] in CATALOGUE:
            return self.send_bytes(200,(ROOT/'workshop/guides'/(parts[1]+'.html')).read_bytes(),'text/html; charset=utf-8')
        if len(parts)==4 and parts[:2]==['api','example'] and parts[2] in CATALOGUE:
            files={f['file'] for f in CATALOGUE[parts[2]]['cases']}
            if parts[3] in files:
                return self.send_bytes(200,(ROOT/parts[2]/parts[3]).read_bytes(),'application/json; charset=utf-8')
        if len(parts)>=3 and parts[0]=='file' and parts[1] in CATALOGUE:
            directory=(ROOT/parts[1]).resolve()
            asset=directory.joinpath(*parts[2:]).resolve()
            if asset.is_relative_to(directory) and not any(part.startswith('.') for part in parts[2:]):
                if asset.is_dir(): asset=asset/'README.md'
                if asset.is_file() and '__pycache__' not in asset.parts:
                    kind={'.py':'text/plain','.md':'text/plain','.json':'application/json'}.get(asset.suffix) or mimetypes.guess_type(asset.name)[0] or 'application/octet-stream'
                    return self.send_bytes(200,asset.read_bytes(),kind)
        if len(parts)==3 and parts[0]=='source' and parts[1] in CATALOGUE:
            resources=CATALOGUE[parts[1]].get('resources',{})
            if parts[2] in resources:
                asset=ROOT/parts[1]/resources[parts[2]]
                kind='application/json' if asset.suffix=='.json' else 'text/plain'
                return self.send_bytes(200,asset.read_bytes(),kind+'; charset=utf-8')
        self.json(404,{'error':'No such workshop page.'})

    def do_POST(self):
        parts=urlsplit(self.path).path.strip('/').split('/')
        origin=self.headers.get('Origin')
        if not self.local_host() or origin not in (None,'http://'+self.headers.get('Host','')):
            return self.json(403,{'error':'This workshop accepts requests from its own local page.'})
        if len(parts)!=3 or parts[:2]!=['api','run'] or parts[2] not in CATALOGUE:
            return self.json(404,{'error':'Unknown application.'})
        if self.headers.get_content_type()!='application/json':
            return self.json(415,{'error':'Supply a JSON case.'})
        try:
            size=int(self.headers.get('Content-Length','0'))
            if not 0<size<=MAX_INPUT:
                return self.json(413,{'error':'A case must be between 1 byte and 512 KiB.'})
            raw=self.rfile.read(size)
            payload=json.loads(raw)
            if not isinstance(payload,dict):
                raise ValueError('The case must be a JSON object.')
            if not CAPACITY.acquire(blocking=False):
                return self.json(429,{'error':'Two calculations are active. Retry when one finishes.'})
            try:
                answer=execute(parts[2],payload,raw=raw)
            finally:
                CAPACITY.release()
            self.json(200,answer)
        except subprocess.TimeoutExpired:
            self.json(422,{'error':'The 60-second interactive budget expired. Reduce the case or run the CLI.'})
        except (ValueError,TypeError,RecursionError) as error:
            self.json(400,{'error':str(error)[:500]})

    def log_message(self, fmt,*args):
        pass

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--port',type=int,default=0,help='local port, default chooses an available port')
    parser.add_argument('--no-browser',action='store_true')
    args=parser.parse_args()
    if sys.version_info<(3,11):
        parser.exit(1,'Install Python 3.11 or newer to use this workshop.\n')
    if not 0<=args.port<=65535:
        parser.error('Port must be 0–65535')
    with ThreadingHTTPServer(('127.0.0.1',args.port),Handler) as server:
        url=f'http://127.0.0.1:{server.server_address[1]}/'
        print(f'Open {url} in your browser. Press Ctrl+C here to close the workshop.',flush=True)
        if not args.no_browser: webbrowser.open(url)
        try: server.serve_forever()
        except KeyboardInterrupt: pass

if __name__=='__main__':main()
