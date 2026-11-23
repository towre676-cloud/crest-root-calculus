const assert=require('node:assert/strict');
const bench=require('./bench.js');
const spec=require('../kit/design.json');
const limits=spec.acceptance;
assert.equal(bench.revision,spec.revision);
for(const [key,value] of Object.entries(bench.lengths)) assert.equal(value,spec.spacings[key]);
assert.equal(spec.spacings.acceptance_plus_minus,.25);
assert.deepEqual([limits.calibration_rule_length,limits.calibration_rule_tol,limits.swing_deg,limits.stroke_angle_ratio,
  limits.angle_repetitions,limits.cycles_each_direction,limits.minimum_measured_unintended_gap,spec.joint.endplay],
  [50,.25,[57,63],[1.65,1.95],3,10,.5,[.15,.35]]);
// Simulated form values exercise the decision logic; these are not a prototype test.
const r={revision:bench.revision,buildId:'simulated test',inspectionDate:'2026-11-23',instruments:'Synthetic instrument fixture',
 printError:.03,lengthError:.02,gapError:.02,angleError:.5,
 scaleBaseX:50,scaleBaseY:50,scaleLinksX:50,scaleLinksY:50,scaleGaugeX:50,scaleGaugeY:50,
 ...bench.lengths,endplayO:.25,endplayD:.25,endplayA:.25,endplayB:.25,goBlade:.55,limitingPair:'crank and rocker',limitingPhase:94.5,contact:'none',cyclesCW:10,cyclesCCW:10,
 fit:'yes',layers:'yes',operation:'yes',retention:'yes'};
for(let i=0;i<3;i++) Object.assign(r,{['rockerLow'+i]:101.334,['rockerHigh'+i]:161.319,['crankAtLow'+i]:56.61,['crankAtHigh'+i]:288.08});
const result=changes=>bench.check({...r,...changes});
assert.equal(result({}).state,'pass');
assert.equal(result({ground:94.05}).state,'review');
assert.equal(result({ground:95}).state,'fail');
assert.equal(result({goBlade:.52}).state,'pass');
assert.equal(result({goBlade:.5}).state,'review');
assert.equal(result({cyclesCW:9}).state,'incomplete');
assert.equal(result({endplayD:' '}).state,'incomplete');
assert.equal(result({crankAtHigh2:200}).state,'fail');
assert.equal(result({crankAtLow2:320,crankAtHigh2:191.47}).state,'pass');
assert.equal(result({crankAtLow0:10,crankAtHigh0:10.4}).state,'review');
assert.equal(result({angleError:-1}).state,'fail');
assert.equal(result({retention:'no'}).state,'fail');
assert.equal(result({revision:'LW-1'}).state,'fail');
assert.equal(result({cyclesCCW:10.5}).state,'fail');
assert.equal(result({scaleBaseX:undefined}).state,'incomplete');
assert.equal(bench.check({revision:bench.revision}).state,'incomplete');
assert.equal(result({goBlade:.48}).state,'review'); // A freely passing blade is a one-sided bound.
assert.equal(result({contact:'observed'}).state,'fail');
assert.equal(result({goBlade:0}).state,'fail');
assert.equal(result({limitingPhase:360}).state,'fail');
assert.equal(result({limitingPhase:0}).state,'pass');
assert.equal(result({scaleGaugeY:50.26}).state,'review');
assert.equal(result({scaleLinksY:49.5}).state,'fail');
assert.equal(result({inspectionDate:'2025-02-29'}).state,'fail');
assert.equal(result({inspectionDate:'2024-02-29'}).state,'pass');
assert.equal(result({lengthError:.3}).state,'review');
for(const key of ['inspectionDate','instruments','printError','lengthError','limitingPair','limitingPhase','contact','goBlade','scaleBaseX','scaleBaseY','scaleLinksX','scaleLinksY','scaleGaugeX','scaleGaugeY']) {
  assert.equal(result({[key]:''}).state,'incomplete',key+' must be required');
}
console.log('Passed physical-record classification and reference-value checks (synthetic test inputs).');
