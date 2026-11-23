/* Local inspection logic for drawing LW-159-R1. Lengths in mm, angles in degrees. */
(function (exports) {
  'use strict';
  const revision = 'LW-159-R1';
  const lengths = {ground:93.8, crank:46.6, coupler:84.1, rocker:111.3};
  const conditions = {
    fit: 'Shafts, bores, collars, flatness, spacing and alignment match the drawing.',
    layers: 'Layer stack and integral pointer match the illustrated assembly branch.',
    operation: 'Completed cycles by hand without a hard stop, rubbing or visible lifting.',
    retention: 'Grip, collars, spacers and shaft positions remained secure after cycling.'
  };
  const format = 'linkwright.physical-inspection/2';
  const sheets = {Base:'Base pattern', Links:'Moving-link pattern', Gauge:'Reversal-angle gauge'};
  function validDay(value) {
    if(typeof value!=='string' || !/^\d{4}-\d{2}-\d{2}$/.test(value) || value.startsWith('0000')) return false;
    const day = new Date(value+'T00:00:00.000Z');
    return Number.isFinite(day.getTime()) && day.toISOString().slice(0,10)===value;
  }
  function decoded(value) {
    if (typeof value !== 'number' && typeof value !== 'string') return null;
    if (typeof value === 'string' && !value.trim()) return null;
    const number = Number(value);
    return Number.isFinite(number) ? number : null;
  }
  function band(lo, hi, lower, upper) {
    // Floor the comparison precision far below any declared instrument resolution.
    const e = 1e-10;
    if (lo > upper + e || hi < lower - e) return 'fail';
    return lo >= lower - e && hi <= upper + e ? 'pass' : 'review';
  }
  function check(record) {
    const rows = [];
    function add(label, state, detail) { rows.push({label, state, detail}); }
    function measured(name, label, low, high, error) {
      const value = decoded(record[name]);
      if (value === null || error === null) { add(label, 'incomplete', 'Reading and uncertainty required.'); return; }
      if (value < 0) { add(label, 'fail', 'A length cannot be negative.'); return; }
      add(label, band(value-error, value+error, low, high),
        (value-error).toFixed(3) + ' to ' + (value+error).toFixed(3) + ' mm');
    }
    function uncertainty(key,label) {
      const n=decoded(record[key]);
      if(n===null || n<=0){ add(label,n===null?'incomplete':'fail','Enter a positive reading-error bound.');return null;}
      return n;
    }
    const printError = uncertainty('printError','Printed-rule reading uncertainty');
    const mm = uncertainty('lengthError','Derived pin-spacing uncertainty');
    const gapError = uncertainty('gapError','Feeler-gauge measurement uncertainty');
    const deg = uncertainty('angleError','Angular measurement uncertainty');
    function requiredText(key,label,missing) {
      const value=record[key];
      add(label,typeof value==='string' && value.trim()?'pass':'incomplete',
        typeof value==='string' && value.trim()?value:missing);
    }
    requiredText('buildId','Build identifier','Give this physical article an identifier.');
    requiredText('instruments','Instruments used','Identify your caliper, angle-reading method and feeler gauges.');
    if(!record.inspectionDate) add('Inspection date','incomplete','Enter the day the physical checks were made.');
    else add('Inspection date',validDay(record.inspectionDate)?'pass':'fail',
      validDay(record.inspectionDate)?record.inspectionDate:'Enter a valid calendar date.');
    if (record.revision!==revision) add('Drawing revision','fail','Record must use '+revision+'.');
    for (const [key,sheet] of Object.entries(sheets)) {
      for(const axis of ['X','Y']) measured('scale'+key+axis, sheet+' · '+axis+' rule',49.75,50.25,printError);
    }
    for (const [key, value] of Object.entries(lengths)) measured(key,key+' centres',value-0.25,value+0.25,mm);
    for (const joint of ['O','D','A','B']) measured('endplay'+joint,'Endplay '+joint,.15,.35,gapError);
    const goBlade=decoded(record.goBlade);
    if(goBlade===null || gapError===null) add('Passing feeler at closest approaches','incomplete','Record a freely passing blade thickness and its error.');
    else if(goBlade<=0) add('Passing feeler at closest approaches','fail','The blade thickness must be positive.');
    else {
      const lower=goBlade-gapError, passes=lower>=0.5-1e-10;
      add('Passing feeler at closest approaches',passes?'pass':'review',
        'Clearance at the checked positions is at least '+(Math.floor((lower+1e-12)*1000)/1000).toFixed(3)+' mm.'+
        (passes?'':' Try a sufficiently thick certified blade. A thinner passing blade leaves the actual gap unresolved.'));
    }
    requiredText('limitingPair','Closest inspected pair','Identify the two physical features at closest observed approach.');
    const phase=decoded(record.limitingPhase);
    add('Crank at closest observed approach',phase===null?'incomplete':phase>=0 && phase<360?'pass':'fail',
      phase===null?'Record the input angle from O→D.':phase>=0 && phase<360?phase+'° counterclockwise from O→D.':'Use degrees from zero up to, but excluding, 360.');
    const contact=record.contact;
    add('Unintended contact during the sweep',contact==='none'?'pass':contact==='observed'?'fail':'incomplete',
      contact==='none'?'Builder completed both sweeps without observing contact.':contact==='observed'?'Observed unintended contact. Correct the interference before accepting.':'Inspect both full sweeps and record any unintended contact.');
    for (let i=0; i<3; i++) {
      const label='Sweep '+(i+1);
      const [p1,p2,t1,t2]=['rockerLow','rockerHigh','crankAtLow','crankAtHigh'].map(k=>decoded(record[k+i]));
      if ([p1,p2,t1,t2,deg].includes(null)) { add(label,'incomplete','Four angles and their reading error are required.');continue; }
      if (p1 < 0 || p2 > 180 || p2 <= p1 || t1 < 0 || t1 >= 360 || t2 < 0 || t2 >= 360) {
        add(label,'fail','Use the illustrated upper branch, ascending rocker angles and crank readings in [0, 360).');continue;
      }
      const swing=p2-p1;
      add(label+' · swing',band(swing-2*deg,swing+2*deg,57,63),
        swing.toFixed(2)+'° with ±'+(2*deg).toFixed(2)+'° combined error');
      const forward=(t2-t1+360)%360, lo=forward-2*deg, hi=forward+2*deg;
      if (lo<=0 || hi>=360) {add(label+' · crank travel','review','The reading interval reaches the angle-wrap boundary. Remeasure both reversals.');continue;}
      const ratio=d=>d/(360-d);
      add(label+' · stroke ratio',band(ratio(lo),ratio(hi),1.65,1.95),ratio(lo).toFixed(3)+' to '+ratio(hi).toFixed(3)+' (increasing/decreasing stroke)');
    }
    for (const dir of ['CW','CCW']) {
      const cycles=decoded(record['cycles'+dir]);
      const valid=Number.isSafeInteger(cycles) && cycles>=0;
      add('Full turns '+dir,cycles===null?'incomplete':!valid?'fail':cycles<10?'incomplete':'pass',
        cycles===null?'Count completed turns.':!valid?'Use a whole, nonnegative count.':cycles+' of 10 required');
    }
    for (const [k,label] of Object.entries(conditions)) {
      const v=record[k];
      add(label,v==='yes'?'pass':v==='no'?'fail':'incomplete',v==='yes'?'Builder confirmed.':v==='no'?'Builder observed a failure.':'Inspection pending.');
    }
    const state=['fail','review','incomplete'].find(s=>rows.some(r=>r.state===s)) || 'pass';
    return {revision, scope:'Builder-entered inspection of one unloaded reference article',state,rows};
  }
  Object.assign(exports,{revision,format,lengths,conditions,sheets,validDay,check});
})(typeof module!=='undefined' && module.exports ? module.exports : (window.LinkwrightBench={}));
