'use strict';
const bench=window.LinkwrightBench, form=document.getElementById('bench'), key='linkwright-inspection:'+bench.revision;
const reads=new Set(['buildId','inspectionDate','instruments','limitingPair','contact','notes']);
const oldFormat='linkwright.physical-inspection/1';
let unsaved=false, legacyMeasurements=null;
function node(tag,text){const n=document.createElement(tag);if(text!==undefined)n.textContent=text;return n;}
function field(parent,name,title,hint){
  const l=node('label',title),n=node('input');n.type='number';n.step='any';n.inputMode='decimal';n.name=name;
  l.append(n);if(hint)l.append(node('small',hint));parent.append(l);reads.add(name);
}
const section=id=>document.getElementById(id);
for(const row of [
  ['printError','Print-rule error','± mm per 50 mm rule reading'],
  ['lengthError','Derived pin-spacing error','± mm on the calculated centre spacing'],
  ['gapError','Feeler/endplay error','± mm, largest applicable bound'],
  ['angleError','Angle reading error','± degrees per direction reading']
]) field(section('uncertainty'),...row);
for(const [key,sheet] of Object.entries(bench.sheets)) {
  for(const axis of ['X','Y']) field(section('lengths'),'scale'+key+axis,sheet+' · '+axis,'50.00 mm');
}
for(const row of [['ground','Base O–D','93.80 mm'],['crank','Crank O–A','46.60 mm'],['coupler','Coupler A–B','84.10 mm'],['rocker','Rocker D–B','111.30 mm']])field(section('lengths'),...row);
for(const j of ['O','D','A','B'])field(section('gaps'),'endplay'+j,'Joint '+j+' endplay','0.15–0.35 mm');
field(section('gaps'),'goBlade','Verified freely passing blade','Thickness in mm; confirmed at all listed checks');
field(section('gaps'),'limitingPhase','Crank at closest approach','0–359.99…°, counterclockwise from O→D');
for(let i=0;i<3;i++){
  const f=node('fieldset'),g=node('div');g.className='fields';f.append(node('legend','Pass '+(i+1)),g);
  for(const [k,t] of [['rockerLow','Rocker at low stop'],['crankAtLow','Crank at low stop'],['rockerHigh','Rocker at high stop'],['crankAtHigh','Crank at high stop']])field(g,k+i,t.replace('stop','reversal'),'degrees');
  section('sweeps').append(f);
}
for(const d of ['CW','CCW'])field(section('cycles'),'cycles'+d,d==='CW'?'Clockwise turns':'Counterclockwise turns','10 complete turns');
for(const [k,t] of Object.entries(bench.conditions)){
  const l=node('label',t),s=node('select');l.className='observed';s.name=k;
  for(const [v,m] of [['','Untested'],['yes','Pass'],['no','Fail']]){const opt=node('option',m);opt.value=v;s.append(opt);}
  l.append(s);section('conditions').append(l);reads.add(k);
}
function current(){
  const record={revision:bench.revision};
  for(const key of reads) record[key]=form.elements.namedItem(key).value;
  if (legacyMeasurements) record.legacyMeasurements=legacyMeasurements;
  return record;
}

function prepare(data, format) {
  if(!data || data.revision!==bench.revision)
    throw Error('The file must contain a '+bench.revision+' physical inspection.');
  const values={};
  for(const key of reads) {
    const value=data[key];
    if(value===undefined || value==='') continue;
    if(typeof value!=='string' && typeof value!=='number') throw Error('Invalid measurement: '+key);
    const input=form.elements.namedItem(key);
    if(input.type==='number' && (!Number.isFinite(Number(value)) || !String(value).trim()))
      throw Error('Expected a finite number for '+key+'.');
    if(input.type==='date' && !bench.validDay(String(value))) throw Error('Invalid calendar date in the saved record.');
    if(input.tagName==='SELECT' && !Array.from(input.options).some(o=>o.value===String(value)))
      throw Error('Invalid choice for '+key+'.');
    if(input.maxLength>0 && String(value).length>input.maxLength)
      throw Error('Value is too long for '+key+'.');
    values[key]=String(value);
  }
  const historical={};
  const prior=format===oldFormat?data:data.legacyMeasurements;
  if(prior) {
    for(const k of ['scaleX','scaleY','gap','lengthError']) {
      if(prior[k]!==undefined && prior[k]!=='' && ['string','number'].includes(typeof prior[k]))
        historical[k]=String(prior[k]).slice(0,100);
    }
  }
  if(format===oldFormat) delete values.lengthError;
  return {values,historical:Object.keys(historical).length?historical:null};
}

function fill(prepared) {
  form.reset();
  for(const [k,value] of Object.entries(prepared.values)) form.elements.namedItem(k).value=value;
  legacyMeasurements=prepared.historical;
}
function hasEntries(){return Array.from(reads).some(k=>form.elements.namedItem(k).value!=='') || !!legacyMeasurements;}
function confirmReplacement(){
  return !unsaved || !hasEntries() || window.confirm('Replace this unsaved physical inspection? Cancel to save its JSON first.');
}
const messages={pass:'Recorded acceptance checks passed',fail:'Build needs correction',review:'Remeasure before accepting',incomplete:'Inspection incomplete'};
let evaluated=false;
function render(){const report=bench.check(current()), panel=section('result');panel.dataset.state=report.state;
  const table=node('table'),head=table.createTHead().insertRow(),body=table.createTBody();for(const t of ['Check','Status','Reading or next action'])head.append(node('th',t));
  for(const row of report.rows){const tr=body.insertRow();for(const text of [row.label,row.state,row.detail])tr.append(node('td',text));}
  panel.replaceChildren(node('h2',messages[report.state]),node('p','Applies only to the identified, unloaded unit and the measurements recorded here.'),table);
  if(legacyMeasurements) {
    const title={scaleX:'aggregate X rule',scaleY:'aggregate Y rule',gap:'reported free gap',lengthError:'old combined spacing/print error'};
    const readings=Object.entries(legacyMeasurements).map(([k,v])=>title[k]+' '+v+' mm');
    panel.append(node('p','Earlier-format readings retained for reference: '+readings.join(', ')+'. The current assessment uses the new measurements above.'));
  }
}
function remember(){
  try {localStorage.setItem(key,JSON.stringify({format:bench.format,readings:current(),unsaved}));} catch {}
}
form.addEventListener('input',()=>{unsaved=true;remember();if(evaluated)render();});
form.addEventListener('submit',e=>{e.preventDefault();evaluated=true;render();section('result').focus();});
section('clear').onclick=()=>{
  if(!confirmReplacement())return;
  form.reset();legacyMeasurements=null;unsaved=false;remember();evaluated=false;
  section('io').textContent='Started a blank inspection.';
  delete section('result').dataset.state;
  section('result').replaceChildren(node('h2','Awaiting bench readings'));
};
section('save').onclick=()=>{
  const r=current(),data={format:bench.format,savedUtc:new Date().toISOString(),readings:r,assessment:bench.check(r)};
  const url=URL.createObjectURL(new Blob([JSON.stringify(data,null,2)+'\n'],{type:'application/json'}));const a=node('a');
  a.href=url;
  a.download='linkwright-'+(r.buildId.toLowerCase().replace(/[^a-z0-9]+/g,'-').replace(/^-|-$/g,'').slice(0,50)||'inspection')+'.json';
  document.body.append(a);a.click();a.remove();
  setTimeout(()=>URL.revokeObjectURL(url),1000);
  unsaved=false;remember();section('io').textContent='Requested a JSON download of the readings and current assessment.';
};
section('open').onchange=async e=>{
  const f=e.target.files[0];if(!f)return;
  try{
    if(f.size>100000)throw Error('The inspection file exceeds 100 kB.');
    const data=JSON.parse(await f.text());
    if(![bench.format,oldFormat].includes(data.format))throw Error('Choose a Linkwright physical inspection JSON.');
    const prepared=prepare(data.readings,data.format);
    if(!confirmReplacement()){section('io').textContent='Kept current inspection.';return;}
    fill(prepared);unsaved=data.format===oldFormat;remember();evaluated=true;render();
    section('io').textContent=data.format===oldFormat?
      'Loaded the older inspection. Its prior print, gap and spacing-error readings remain in the JSON history. Complete the sheet-by-sheet checks, passing-blade witness and observation details for this assessment.':
      'Loaded readings and recalculated the result.';
  }catch(err){section('io').textContent=err.message;}
  finally{e.target.value='';}
};
section('print').onclick=()=>{evaluated=true;render();window.print();};
window.addEventListener('beforeprint',()=>{evaluated=true;render();});
try{
  const stored=localStorage.getItem(key);
  if(stored) {
    const data=JSON.parse(stored);
    const modern=data.format===bench.format;
    fill(prepare(modern?data.readings:data,modern?bench.format:oldFormat));
    unsaved=modern?Boolean(data.unsaved):hasEntries();
    if(hasEntries())section('io').textContent=modern?'Restored this inspection from the browser.':'Restored older evidence. Complete the new checks before acceptance.';
  }
}catch{}
