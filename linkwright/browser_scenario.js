// Executed in the live page's JavaScript context by test_page.py.
(async () => {
  const get=id=>document.getElementById(id);
  const click=id=>get(id).click();
  const assert=(pred,msg)=>{if(!pred) throw Error(msg);};
  const change=(id,value,event)=>{
    get(id).value=value;
    get(id).dispatchEvent(new Event(event,{bubbles:true}));
  };
  const wait=ms=>new Promise(resolve=>setTimeout(resolve,ms));
  const baseline=window.__linkwrightSelfTest();
  assert(baseline.length===4 && baseline.every(v=>v.ok),'Original flap parameter cases');
  assert(get('closureText').textContent.includes('Every length combination closes'),'Default full rotation');
  assert(designRecord().workshop.roundedConfiguration.rockerPositions===2,'Exact assembly fibre');
  change('anglePick','1/0','change');
  assert(designRecord().workshop.roundedConfiguration.crankAnglePair[1]===0,'Half-turn projective coordinate');
  change('q',180,'input');click('synthFlap');change('tol',5,'input');
  assert(get('closureText').textContent.includes('cannot close'),'Tolerance obstruction');
  change('driveSense','-1','change');
  assert(designRecord().driveDirection==='clockwise' && get('checks').textContent.includes('Turning the input clockwise'),'Reversible drive direction');
  change('driveSense','1','change');
  const below=designRecord();
  assert(below.version===3 && below.cutResolutionMm===.1,'Manufacturing model version');
  const expectedG=below.workshop.roundedConfiguration.ground/10;
  assert(Math.abs(Math.hypot(below.rockerPivot.x-below.crankPivot.x,below.rockerPivot.y-below.crankPivot.y)-expectedG)<1e-9,'Printed ground matches export');
  assert(below.workshop.fitScales.min===160,'Rounded scale cutoff');
  change('scale',159,'input');
  assert(designRecord().workshop.closure.state==='at limit','Rounded legacy 159% hits toggle');
  change('tol',.25,'input'); click('loadCut');
  get('asBuilt').requestSubmit();
  assert(get('actualResult').textContent.includes('59.98°'),'As-built reference rocker swing');
  change('actualB','40','input');get('asBuilt').requestSubmit();
  assert(get('actualResult').textContent.includes('Full turn: fails'),'Measured failure must revoke pass');
  change('poseAngle','92.5','input');
  assert(get('play').textContent==='Play' && get('poseValue').textContent==='92.5°','Pose inspection pauses input');
  change('tol',5,'input'); change('scale',100,'input');
  assert(below.workshop.fitScales.min>100,'Inverse scale must increase design');
  click('applyScale');
  assert(get('closureText').textContent.includes('Every length combination closes'),'Suggested scale fixes folding');
  const repaired=designRecord();
  change('layer','shared','change');
  assert(get('sweepText').textContent.includes('one layer'),'Shared-layer rendering');
  change('mat','stick','change');
  assert(get('applyScale').disabled && !designRecord().workshop.fitScales.count,'Stock too short');
  change('mat','card','change');change('tol',.5,'input');change('scale',100,'input');
  click('tabPath');click('samplePath');click('synthPath');click('tabFlap');
  await wait(200);
  assert(designRecord().mode==='flap' && !get('synthPath').disabled,'Cancelled stale search');
  click('tabPath');click('samplePath');click('synthPath');
  for(let i=0;i<1800 && get('synthPath').disabled;i++) await wait(100);
  assert(get('status').textContent.startsWith('Best fit:'),'Path synthesis did not finish');
  const path=designRecord();
  assert(path.pathGap && path.pathGap.fraction<.12 && path.workshop.closure.state==='closes','Geometric fit');
  const gap=path.pathGap.mean;
  change('scale',145,'input');
  const scaled=designRecord();
  assert(scaled.pathGap.dimensionUnits==='mm at current cut scale','Scaled fit has current units');
  assert(Math.abs(scaled.pathGap.mean-1.45*gap)<.2,'Scaled fit recomputed from cut geometry');
  assert(get('partsWrap').querySelector('path[d*=" A "]'),'Path plate has a cuttable outside contour');
  assert(scaled.requested.scaledTargetMm.length===64,'Export includes target coordinates');
  change('scale',100,'input');
  assert(get('partsWrap').querySelector('svg').getAttribute('width').endsWith('mm'),'Scaled, printable SVG');
  assert(get('partsWrap').textContent.includes('Pin spacing'),'Parts sheet carries workshop settings');
  return {baseline,timingExample:{swingDeg:60,ratio:1.8,toleranceMm:5,
    before:below.workshop.closure.state,scalePercent:below.workshop.fitScales.min,after:repaired.workshop.closure.state},
    exampleTrace:{...path.pathGap,sweep:path.workshop.sweep.state,scaleRange:path.workshop.fitScales}};
})()
