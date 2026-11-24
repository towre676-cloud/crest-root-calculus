#!/usr/bin/env python3
"""Check the v3 viewer's bindings and interactions without claiming browser QA.

Python checks are standard-library only. When Node is present, a deliberately
small DOM/Canvas adapter executes the actual embedded application script. It
checks rendering paths and controls, not CSS, pixels, accessibility conformance
or actual mobile/touch behavior. It never downloads a browser.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

import build_viewer


ROOT = Path(__file__).resolve().parent


ADAPTER = r'''
const fs=require('fs'),vm=require('vm'),assert=require('assert');
const html=fs.readFileSync(process.argv[2],'utf8');
const data=html.match(/<script id="payload" type="application\/json">([\s\S]*?)<\/script>/)[1];
const code=html.match(/<script>\n([\s\S]*?)<\/script>/)[1];
const D=JSON.parse(data),results=[];
for(const viewport of ['desktop_adapter','narrow_adapter']){
 const ids=new Map(),groups={phaseIndex:[]},timers=new Map();let timerId=0,drawCalls=0;
 const context2d=new Proxy({measureText:s=>({width:String(s).length*6})},{get:(obj,key)=>{
  if(key in obj)return obj[key];return (...args)=>{if(['moveTo','lineTo','arc','fillRect','setTransform'].includes(key)){assert(args.every(v=>typeof v!=='number'||Number.isFinite(v)),'Nonfinite canvas coordinate: '+key);drawCalls++;}};
 }});
 class Element{
  constructor(tag='div',id=''){this.tagName=tag;this.id=id;this.dataset={};this.events={};this.style={};this.checked=false;this.value='';this.textContent='';this.width=0;this.height=0;this.attrs={};}
  get clientWidth(){return viewport==='narrow_adapter'?356:this.id==='map'?440:820;}
  get clientHeight(){return this.id==='scene'?(viewport==='narrow_adapter'?380:470):this.id==='map'?343:this.id==='side'?220:100;}
  getBoundingClientRect(){return {width:this.clientWidth,height:this.clientHeight,left:0,top:0};}
  getContext(){return context2d;}setAttribute(k,v){this.attrs[k]=v;}addEventListener(k,v){this.events[k]=v;}
  setPointerCapture(){}prepend(){}
  set innerHTML(v){this._innerHTML=v;if(this.id==='phaseButtons'){groups.phaseIndex=[];parse(v);}}
  get innerHTML(){return this._innerHTML||'';}
 }
 function parse(text){for(const match of text.matchAll(/<([a-z][a-z0-9]*)\b([^>]*)>/gi)){
  const a=match[2],id=a.match(/\bid="([^"]*)"/)?.[1]||'',e=new Element(match[1],id);if(id)ids.set(id,e);
  e.checked=/\bchecked\b/.test(a);for(const m of a.matchAll(/data-([a-z-]+)="([^"]*)"/g)){const key=m[1].replace(/-([a-z])/g,(_,c)=>c.toUpperCase());e.dataset[key]=m[2];if(groups[key])groups[key].push(e);}
 }}
 parse(html.split('<script id="payload"')[0]);const payload=new Element('script','payload');payload.textContent=data;ids.set('payload',payload);
 const document={getElementById:id=>ids.get(id),querySelectorAll:selector=>selector==='[data-phase-index]'?groups.phaseIndex:[],createElement:tag=>new Element(tag),querySelector:()=>new Element()};
 const win={devicePixelRatio:1.5,events:{},addEventListener(k,v){this.events[k]=v;}};
 const context=vm.createContext({console,document,window:win,setInterval:fn=>{timers.set(++timerId,fn);return timerId;},clearInterval:id=>timers.delete(id)});
 new vm.Script(code,{filename:'workbench.js'}).runInContext(context);
 assert(win.__workbenchReady);assert(!win.__workbenchError);assert.equal(win.__workbenchState.designSha256,D.designSha256);
 assert.equal(win.__workbenchState.frameCount,D.scene.frames.length);assert.equal(win.__workbenchState.memberCount,D.scene.members.length);
 assert.equal(ids.get('openingMetric').textContent,'Unresolved');assert(ids.get('openingBoundary').textContent.includes('not authorized'));
 const slider=ids.get('phase');let frames=0;
 for(let i=0;i<D.scene.frames.length;i++){slider.value=String(i);slider.events.input();assert.equal(win.__workbenchState.frameIndex,i);assert.equal(win.__workbenchState.stageId,D.scene.frames[i].stage_id);frames++;}
 for(const button of groups.phaseIndex){button.events.click();assert.equal(win.__workbenchState.frameIndex,Number(button.dataset.phaseIndex));}
 for(const id of ['parcels','canopy','widths','access','swept']){const e=ids.get(id);e.checked=true;e.events.change();e.checked=false;e.events.change();}
 for(const id of ['sceneZoomIn','sceneZoomOut','sceneReset','mapZoomIn','mapZoomOut','mapDistrict','mapFit'])ids.get(id).events.click();
 for(const id of ['scene','map']){
  const e=ids.get(id);for(const key of ['ArrowLeft','ArrowRight','ArrowUp','ArrowDown','+','-'])e.events.keydown({key,preventDefault(){}});
  e.events.pointerdown({clientX:100,clientY:100,pointerId:1});e.events.pointermove({clientX:140,clientY:120,pointerId:1});e.events.pointerup({});e.events.pointercancel({});
 }
 ids.get('map').events.wheel({deltaY:-1,preventDefault(){}});win.events.resize();
 slider.value='0';slider.events.input();ids.get('play').events.click();assert.equal(ids.get('play').attrs['aria-pressed'],'true');
 for(const fn of timers.values())fn();assert.equal(win.__workbenchState.frameIndex,1);ids.get('play').events.click();assert.equal(timers.size,0);
 assert(ids.get('evidenceRows').innerHTML.includes('ambiguous'));assert(ids.get('evidenceRows').innerHTML.includes('model_mismatch'));
 results.push({adapter:viewport,frames_executed:frames,phase_buttons:groups.phaseIndex.length,checkbox_layers:5,view_controls:7,keyboard_and_pointer_handlers:true,play_pause:true,finite_canvas_arguments:true,draw_calls:drawCalls,runtime_exceptions:0});
}
console.log(JSON.stringify({engine:'Node vm with a minimal explicit DOM/Canvas2D adapter',results,real_browser_tested:false,css_layout_tested:false,mobile_browser_tested:false}));
'''


def main():
    html_path = ROOT / "workbench.html"
    html = html_path.read_text(encoding="utf-8")
    payload = json.loads(re.search(r'<script id="payload" type="application/json">([\s\S]*?)</script>', html).group(1))
    report = build_viewer.read_json(ROOT / "integrated-report.json")
    parts = build_viewer.integrated_components(report)
    checks = {"all_six_design_bindings_match": True}
    for name, expected in payload["sourceHashes"].items():
        if build_viewer.sha256(ROOT / name) != expected:
            raise AssertionError("Embedded source is stale: " + name)
    checks["embedded_input_file_hashes_match"] = True
    checks["exactly_395_source_buildings"] = len(payload["atlas"]["buildings"]) == 395
    checks["same_motion_frame_count"] = len(payload["scene"]["frames"]) == len(parts["motion"]["frames"])
    checks["same_member_count"] = len(payload["scene"]["members"]) == parts["assembly"]["assembly"]["physical_member_count"]
    checks["canonical_nodes_preserved"] = all(view["world_nodes_m"] == source["world_nodes_m"]
        for view, source in zip(payload["scene"]["frames"], parts["motion"]["frames"]))
    checks["network_requests_blocked_by_csp"] = "connect-src 'none'" in html
    checks["no_external_script_or_style_resource"] = not re.search(r'<(?:script|link)\b[^>]*(?:src|href)\s*=', html, re.I)
    script = re.search(r'<script>\n([\s\S]*?)</script>', html).group(1)
    checks["no_network_api_in_application"] = not re.search(r'\b(?:fetch|XMLHttpRequest|WebSocket|EventSource)\b', script)
    local_links = [value for value in re.findall(r'href="([^"]+)"', html)
                   if not value.startswith(("#", "https:", "http:"))]
    checks["local_reading_links_exist"] = all((ROOT / value).exists() for value in local_links)
    mismatch = deepcopy(report)
    container = mismatch.get("reports", mismatch)
    motion_key = "motion" if "motion" in container else "motion_report"
    container[motion_key]["design_sha256"] = "0" * 64
    try:
        build_viewer.integrated_components(mismatch)
    except ValueError:
        checks["mismatched_design_rejected"] = True
    else:
        checks["mismatched_design_rejected"] = False
    if not all(checks.values()):
        raise AssertionError({key: value for key, value in checks.items() if not value})
    node = shutil.which("node")
    javascript = {"status": "not_run", "reason": "Node is unavailable"}
    if node:
        with tempfile.TemporaryDirectory(prefix="crest-v3-viewer-review-") as directory:
            path = Path(directory)
            (path / "workbench.js").write_text(script, encoding="utf-8")
            (path / "adapter.cjs").write_text(ADAPTER, encoding="utf-8")
            subprocess.run([node, "--check", str(path / "workbench.js")], check=True, capture_output=True, text=True)
            result = subprocess.run([node, str(path / "adapter.cjs"), str(html_path)], check=True,
                                    capture_output=True, text=True, timeout=45)
            javascript = {"status": "passed", "syntax": "passed", **json.loads(result.stdout)}
    browsers = {name: shutil.which(name) for name in ("chromium", "chromium-browser", "google-chrome", "firefox")}
    receipt = {
        "format": "crest-root-integrated-viewer-review/3", "design_sha256": report["design_sha256"],
        "checks": checks, "javascript": javascript, "browser_executables_on_path": browsers,
        "real_browser_rendering_verified": False, "mobile_layout_verified": False,
        "windows_launcher_verified": False,
        "limitation": "Source/binding checks and simulated DOM/Canvas interactions only; no actual browser pixels, CSS layout or mobile touch verification. No browser was installed or downloaded.",
        "file_sha256": {name: build_viewer.sha256(ROOT / name) for name in
            ("build_viewer.py", "workbench.html", "integration-review.md", "independent_review.py", "integrated-report.json")},
    }
    (ROOT / "verification-viewer.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"checks_passed": len(checks), "javascript": javascript["status"],
                      "physical_approval": False, "real_browser_verified": False}))


if __name__ == "__main__":
    main()
