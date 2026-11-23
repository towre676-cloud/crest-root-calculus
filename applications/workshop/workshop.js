'use strict';
const $=id=>document.getElementById(id);
const appNav=$('apps'), editor=$('editor'), out=$('out');
let catalogue, active, phase=0, revision=0, latest=null;
function el(tag,text='',cls=''){
  const e=document.createElement(tag); if(text!=='') e.textContent=text; if(cls)e.className=cls; return e;
}
function svg(tag,attrs={}){
  const e=document.createElementNS('http://www.w3.org/2000/svg',tag);
  for(const [k,v] of Object.entries(attrs))e.setAttribute(k,String(v)); return e;
}
function download(name,obj){
  const url=URL.createObjectURL(new Blob([JSON.stringify(obj,null,2)+'\n'],{type:'application/json'}));
  const a=el('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
}
function message(text){$('activity').textContent=text;}
function compact(n){return typeof n==='number'?new Intl.NumberFormat('en',{maximumFractionDigits:5}).format(n):String(n??'—');}
function human(s){return String(s??'').replaceAll('_',' ');}
function range(a){return Array.isArray(a)?a.map(compact).join(' – '):'—';}
function line(text){out.append(el('p',text,'result-lead'));}
function table(headers,rows){
  if(!rows.length)return;
  const wrap=el('div','','table-wrap'),t=el('table'),head=el('thead'),tr=el('tr'),body=el('tbody');
  headers.forEach(h=>{const th=el('th',h);th.scope='col';tr.append(th)});head.append(tr);t.append(head);
  rows.forEach(row=>{const r=el('tr');row.forEach(v=>r.append(el('td',v??'—')));body.append(r)});t.append(body);wrap.append(t);out.append(wrap);
}
function bars(rows,caption){
  if(!rows.length)return;const f=el('figure','','figure'),group=el('div','','bars');
  const max=Math.max(0,...rows.map(r=>Number(r[1])));
  for(const [name,value] of rows){
    const row=el('div','','barrow'),track=el('div','','track'),bar=el('div','','bar');
    bar.style.width=(max?100*Math.max(0,value)/max:0)+'%';track.append(bar);
    row.append(el('span',name),track,el('strong',compact(value)));group.append(row);
  }
  f.append(group,el('figcaption',caption));out.append(f);
}
function pointPatterns(candidates){
  if(!candidates?.length)return;
  const f=el('figure','','figure'), list=candidates.slice(0,6),width=Math.max(1,...list.map(p=>Math.max(...p)));
  const s=svg('svg',{viewBox:'0 0 620 '+(list.length*40+24),class:'result-chart',role:'img','aria-label':'Candidate site arrangements shown on a common position axis.'});
  for(const [row,points] of list.entries()){
    const y=20+row*40;
    s.append(svg('line',{x1:34,x2:590,y1:y,y2:y,stroke:'#a9bcb2'}));
    for(const p of points){const x=34+p/width*556; s.append(svg('circle',{cx:x,cy:y,r:4,fill:row%2?'#994929':'#00645d'}));const t=svg('text',{x,y:y+16,'font-size':11,'text-anchor':'middle'});t.textContent=p;s.append(t);}
  }
  f.append(s,el('figcaption','Positions in grid units. Each row is a distinct arrangement up to translation and reflection. Showing '+list.length+' of '+candidates.length+'.'));out.append(f);
}
function addWitnesses(ws){
  if(!ws?.length)return;
  table(['Example profile','Library row','Bottom / middle / top width (nm)'],
    ws.map(w=>[w.class,w.row,w.bottom_middle_top_nm.map(compact).join(' / ')]));
}
function subsequent(k,r,input){
  let title, fields, apply;
  if(k==='wafer' && r.next_reading && input.observations){
    const n=r.next_reading.channel; title='Record channel '+n;fields=['Normalized intensity'];
    apply=v=>{input.observations[n]=v[0]};
  } else if(k==='battery' && r.next_checkpoint_s){
    title='Record the '+r.next_checkpoint_s+' s checkpoint'; fields=['Measured voltage (V)'];
    apply=v=>{input.observations[r.next_checkpoint_s]=v[0]};
  } else if(k==='water' && r.next_logger){
    title='Record the reading at '+r.next_logger;fields=['Lower pressure change ('+(input.pressure_unit||'input units')+')','Upper pressure change'];
    apply=v=>{input.observations[r.next_logger]=v};
  } else if(k==='retrofit' && r.next_test){
    const test=input.tests.find(t=>t.id===r.next_test);
    title='Record '+r.next_test;fields=['Lower loss (W/K)','Upper loss (W/K)'];
    apply=v=>{(input.measured??=[]).push({quantity:test.quantity,bounds_W_K:v})};
  } else if(k==='materials' && r.next_probe?.sites){
    title='Record the local count';fields=['Groups observed'];
    apply=v=>{if(!Number.isSafeInteger(Number(v[0]))||Number(v[0])<0)throw Error('Enter a nonnegative whole count.');
      (input.observations??=[]).push({sites:r.next_probe.sites,span:r.next_probe.span,count:Number(v[0])})};
  } else if(k==='orbits' && r.recommended_window){
    title='Record an observed transit midpoint';fields=['Earliest time (input day scale)','Latest time'];
    apply=v=>{(input.events??=[]).push(v)};
  } else return;
  const form=el('form','','follow'),head=el('h3',title),inner=el('div','','controls'),inputs=[];
  form.append(head,inner);
  for(const name of fields){
    const l=el('label',name),i=el('input');i.required=true;i.type='number';i.step='any';l.append(i);inner.append(l);inputs.push(i);
  }
  const button=el('button','Apply measurement');button.type='submit';const error=el('span','','warn');
  form.append(button,error);out.append(form);
  form.onsubmit=e=>{e.preventDefault();try{
    const values=inputs.map(i=>i.value);
    if(values.length===2&&Number(values[0])>Number(values[1]))throw Error('Lower value exceeds upper value.');
    apply(values);editor.value=JSON.stringify(input,null,2);fillControls(input);invalidate();calculate();
  }catch(err){error.textContent=err.message}};
}

function render(k,r){
  out.replaceChildren();out.append(el('h2','Calculation result'));out.firstChild.id='answer';
  if(k==='wafer'){
    line(r.decision);bars([['Within width limit',r.within],['Over width limit',r.above]],'Number of compatible simulated profiles. Bars share a scale.');
    if(r.bottom_width_range_nm)line('Surviving bottom widths span '+range(r.bottom_width_range_nm)+' nm.');
    else line('Closest profile’s largest intensity discrepancy: '+compact(r.smallest_maximum_intensity_discrepancy)+'.');
    if(r.next_reading)line('Next reading: channel '+r.next_reading.channel+' · '+r.next_reading.setting+'. It separates '+compact(r.next_reading.opposite_decision_pairs_separated)+' of '+compact(r.opposite_decision_pairs)+' opposite-grade pairs.');
    addWitnesses(r.witnesses);
  } else if(k==='water'){
    line(r.zones?.length?'Search '+r.zones.join(' and ')+'. The surviving reaches cover '+compact(r.search_extent)+' '+r.extent_unit+'.':'No supplied leak scenario explains all the readings.');
    if(r.next_logger)line('Place the next temporary logger at '+r.next_logger+'.');
    bars((r.logger_comparison||[]).slice(0,12).map(x=>[x.sensor,x.worst_remaining_extent]),'Largest search extent that could remain after each logger visit ('+r.extent_unit+'). Shorter is better.');
    table(['Surviving reach','Leak strength range'],(r.surviving_scenarios||[]).map(x=>[x.zone,range(x.strength)]));
  } else if(k==='battery'){
    line(r.candidate_count?'The readings fit '+r.candidate_count+' reference traces. Their capacity interval is '+range(r.capacity_range_ah)+' Ah.':'The measured voltages fall outside every supplied trace allowance.');
    line('Capacity decision: '+human(r.decision)+'. '+(r.next_checkpoint_s?'Next checkpoint at '+r.next_checkpoint_s+' s.':''));
    table(['Elapsed time (s)','Opposite-grade pairs separated','Additional time (s)'],(r.checkpoint_options||[]).map(x=>[x.elapsed_s,x.separated_cross_grade_pairs,x.additional_s]));
    if(r.reason)out.append(el('p',r.reason,'fine'));
  } else if(k==='retrofit'){
    line(r.reason||(r.chosen_upgrade?(r.status==='one_upgrade_never_worse'?'This option never scores lower: ':'The measured ranges favor ')+r.chosen_upgrade+'.':human(r.status)+'.'));
    if(r.next_test)line('Next measurement: '+r.next_test+'.');
    table(['Proposal','Annual heat saved (kWh)','Annual kWh per cost unit'],(r.options||[]).map(x=>[x.id,range(x.heat_saved_kWh_per_year),range(x.annual_kWh_per_cost_unit)]));
    bars((r.test_comparison||[]).map(t=>[t.id,t.guaranteed_width_reduction]),'Guaranteed reduction in uncertainty about the advantage of one proposal over the other (annual kWh per cost unit).');
  } else if(k==='materials'){
    line(human(r.status)+'. '+(r.complete?'The search is complete for the supplied finite model.':'The computation stopped at its search budget. Further arrangements may exist.'));
    if(r.classes){
      const names={P:'Primitive',I:'Body-centered',F:'Face-centered'};
      table(['Compatible cubic indexing','Possible cell edge (Å)'],Object.entries(r.classes).filter(([,v])=>v.length).map(([c,v])=>[names[c]||c,v.map(w=>range(w.a_A_display_approx)).join(' / ')]));
    }else pointPatterns(r.candidates);if(r.next_probe?.question)line('Next local observation: '+r.next_probe.question);
    if(r.next_probe?.predicted_count_per_candidate)table(['Arrangement','Predicted count'],r.next_probe.predicted_count_per_candidate.slice(0,12).map((v,i)=>[i+1,v]));
  } else if(k==='orbits'){
    line(r.surviving_hypotheses?.length?'Surviving orbits: '+r.surviving_hypotheses.join(', ')+'.':'No listed orbit matches the supplied events.');
    if(r.recommended_window)line('Best available window: '+r.recommended_window+'.');
    table(['Observing window','Hours','Predicted outcome by orbit'],(r.windows||[]).map(w=>[w.name,compact(w.hours),Object.entries(w.outcomes).map(([h,v])=>h+' — '+v).join(' / ')]));
  } else if(k==='atlas'){
    const label=s=>({no_elliptic_map_over_Qbar:'All elliptic maps excluded',no_degree_three_map_over_Qbar:'Degree-three maps excluded',degree_three_map_count_certified:'Degree-three count certified',degree_three_elliptic_map_over_Q:'Degree-three elliptic map over Q',degree_three_genus_one_map:'Degree-three genus-one map',unresolved:'Open'}[s]||human(s));
    line(Object.entries(r.counts).map(([k,v])=>v+' · '+label(k)).join('   /   '));
    table(['Source curve','Conclusion','Proof detail'],r.curves.map(x=>[x.id,label(x.status),
      x.proof?'At prime '+x.proof.prime+', lower '+x.proof.lower+' and upper '+x.proof.upper:
        x.reason||((x.reductions||[]).map(p=>'p = '+p.p).join(', ')||'See the full result for the map')]));
  }
  subsequent(k,r,JSON.parse(editor.value));
  const save=el('button','Save calculation record');save.type='button';save.onclick=()=>download('crest-root-'+k+'-record.json',latest);
  const details=el('details'),summary=el('summary','Inspect the full result'),pre=el('pre',JSON.stringify(r,null,2));details.append(summary,pre);out.append(save,details);
  const links=Object.entries(catalogue[k].resources||{});if(links.length){
    const section=el('p','','fine');section.append('Recorded checks: ');links.forEach(([key],i)=>{if(i)section.append(' · ');const a=el('a',human(key.replaceAll('-',' ')));a.href='/source/'+k+'/'+key;a.target='_blank';section.append(a)});out.append(section);
  }
}
function invalidate(){revision++;latest=null;out.replaceChildren();message('Inputs changed. Run the calculation to update the result.');}
function locate(obj,path){const parts=path.split('/');let parent=obj;for(const key of parts.slice(0,-1)){if(parent==null)return [null,''];parent=parent[key]}return [parent,parts.at(-1)];}
function fillControls(input){
  const panel=$('controls');panel.replaceChildren();
  const fields=[...catalogue[active].fields];
  if(active==='water' && input.observations && !input.observations.Inlet){
    for(const name of Object.keys(input.observations).slice(0,2))fields.push({path:'observations/'+name,label:name+' pressure change ('+(input.pressure_unit||'input units')+')',type:'decimals'});
  }
  for(const f of fields){
    const [owner,key]=locate(input,f.path);if(!owner||!(key in owner))continue;
    const l=el('label',f.label),v=el('input');l.append(v);panel.append(l);
    if(f.type==='number'||f.type==='integer'){v.type='number';v.step=f.type==='integer'?'1':'any';v.value=owner[key];}
    else {v.type='text';v.value=Array.isArray(owner[key])?owner[key].join(', '):owner[key];}
    v.oninput=()=>{
      invalidate();
      try{
        const obj=JSON.parse(editor.value),[parent,k]=locate(obj,f.path);let x=v.value;
        if(f.type==='integer'||f.type==='number'){
          if(!x.trim()||!Number.isFinite(Number(x))||(f.type==='integer'&&!Number.isInteger(Number(x))))throw Error('Enter a valid number.');
          x=Number(x);
        }else if(f.type==='list'||f.type==='decimals'){
          const pieces=x.split(',').map(a=>a.trim());if(pieces.some(a=>!a||!Number.isFinite(Number(a))))throw Error('Separate finite numbers with commas.');
          x=f.type==='decimals'?pieces:pieces.map(Number);
        }
        parent[k]=x;editor.value=JSON.stringify(obj,null,2);v.setCustomValidity('');
      }catch(err){v.setCustomValidity(err.message);message(err.message);}
    };
  }
}
async function readJSON(url,opts){
  const resp=await fetch(url,opts);const body=await resp.json();if(!resp.ok)throw Error(body.error||'Calculation failed.');return body;
}
async function calculate(){
  const k=active, stamp=++phase, rev=revision;
  if(!$('controls').querySelectorAll('input:invalid').length){
    try{
      const raw=editor.value, input=JSON.parse(raw);$('run').disabled=true;out.replaceChildren();message('Calculating…');
      const result=await readJSON('/api/run/'+k,{method:'POST',headers:{'Content-Type':'application/json'},body:raw});
      if(stamp!==phase||rev!==revision)return;
      latest={format:'crest-root-application-record/v1',application:k,input,input_json:raw,result};render(k,result);message('Calculation complete.');
    }catch(error){if(stamp===phase){latest=null;out.replaceChildren();message(error.message);}}
    finally{if(stamp===phase)$('run').disabled=false;}
  }else message('Correct the input field before running.');
}
async function loadCase(name){
  const k=active,stamp=++phase;$('run').disabled=true;message('Opening case…');out.replaceChildren();latest=null;
  try{
    const obj=await readJSON('/api/example/'+k+'/'+encodeURIComponent(name));if(stamp!==phase)return;
    editor.value=JSON.stringify(obj,null,2);fillControls(obj);revision++;await calculate();
  }catch(error){if(stamp===phase)message(error.message);}
}
function choose(){
  const h=location.hash.slice(1);active=Object.hasOwn(catalogue,h)?h:Object.keys(catalogue)[0];const a=catalogue[active];
  document.title=a.name+' · Crest–Root';$('domain').textContent=a.name;$('question').textContent=a.question;$('context').textContent=a.context;
  for(const link of appNav.children){if(link.hash==='#'+active)link.setAttribute('aria-current','page');else link.removeAttribute('aria-current');}
  $('guide').href='/guide/'+active;$('guide').target='_blank';$('cases').replaceChildren();
  for(const c of a.cases){const o=el('option',c.label);o.value=c.file;$('cases').append(o);}
  $('advanced').open=false;loadCase(a.cases[0].file);
}
$('cases').onchange=e=>loadCase(e.target.value);$('run').onclick=calculate;
editor.oninput=()=>{invalidate();try{fillControls(JSON.parse(editor.value));}catch{$('controls').replaceChildren();}};
$('save-in').onclick=()=>{try{download(active+'-case.json',JSON.parse(editor.value))}catch(e){message(e.message)}};
$('import').onchange=async e=>{const file=e.target.files[0];if(!file)return;try{
  if(file.size>512*1024)throw Error('File exceeds 512 KiB.');
  const obj=JSON.parse(await file.text());if(!obj||typeof obj!=='object'||Array.isArray(obj))throw Error('Open a JSON case object.');
  editor.value=JSON.stringify(obj,null,2);fillControls(obj);invalidate();message('Case opened. Run it to calculate the result.');
}catch(err){message(err.message)}finally{e.target.value='';}};
(async()=>{try{catalogue=await readJSON('/api/catalogue');for(const [id,info] of Object.entries(catalogue)){
  const a=el('a',info.name);a.href='#'+id;appNav.append(a);}window.addEventListener('hashchange',choose);choose();
}catch(err){message(err.message)}})();
