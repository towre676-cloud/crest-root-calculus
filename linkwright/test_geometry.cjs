/* Run with Node.js. The tests exercise the geometry used by the actual page. */
'use strict';
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const html = fs.readFileSync(path.join(__dirname,'index.html'),'utf8');
const script = html.split('<script>')[1].split('// ---------------------------------------------------------------- state & UI')[0];
const sandbox = {};
vm.runInNewContext(script + '\nglobalThis.api={V,rot,dist,solve,simulate,closureReport,barDistance,sweepReport,synthFlap,traceGap,scaleChoices,roundedCurve,projectiveFibre,evenQuarticJ,goodPrimeCounts,fabricationModel,strokeReport,minimumTransmission,plateOutline};',sandbox);
const {V,dist,simulate,closureReport,barDistance,sweepReport,synthFlap,traceGap,scaleChoices,roundedCurve,projectiveFibre,evenQuarticJ,goodPrimeCounts,fabricationModel,strokeReport,minimumTransmission,plateOutline}=sandbox.api;
const close=(x,y,eps=1e-8)=>assert.ok(Math.abs(x-y)<=eps, x+' vs '+y);
const model=(g,a,b,c,angle=0,branch=1)=>({A0:V(0,0),B0:V(g*Math.cos(angle),g*Math.sin(angle)),a,b,c,s:branch});
const deg=Math.PI/180;
const receipts={};
const hidden=model(100,30,75,54.99994,.25*deg);
assert.equal(simulate(hidden,720).fails,0);
const obstruction=closureReport(hidden);
assert.equal(obstruction.state,'fails');
close(obstruction.witnessDeg,180.25);
close(obstruction.margin,-.00006);
receipts.missedBySampling={positions:720,sampledFailures:0,witnessDeg:obstruction.witnessDeg,triangleMarginMm:obstruction.margin};
const near=model(100,30,75,55.4);
assert.equal(closureReport(near).state,'closes');
close(closureReport(near).limit,.1);
assert.equal(closureReport(near,.09).state,'closes');
assert.equal(closureReport(near,.1).state,'at limit');
assert.equal(closureReport(near,.11).state,'fails');
assert.equal(closureReport(near,30).state,'invalid');
assert.equal(closureReport(near,-1).state,'invalid');
assert.equal(closureReport(model(2,2,4,4)).state,'at limit');
receipts.toleranceExample={lengthsMm:[100,30,75,55.4],limitMm:closureReport(near).limit,adverseLengths:closureReport(near,.11).lengths};
let seed=19070;
const rand=()=>{seed=(Math.imul(1664525,seed)+1013904223)>>>0;return seed/4294967296;};
for(let i=0;i<160;i++){
  const nominal=[0,1,2,3].map(()=>30+rand()*100);
  const err=rand()*3, report=closureReport(model(...nominal),err);
  const z=report.lengths;
  const check=closureReport(model(z.ground,z.crank,z.coupler,z.rocker));
  close(report.margin,check.margin);
  for(let k=0;k<40;k++){
    const sample=nominal.map(x=>x+(2*rand()-1)*err);
    assert.ok(closureReport(model(...sample)).margin>=report.margin-1e-8);
  }
}
close(barDistance(V(0,0),V(3,0),V(1,-2),V(1,2)).gap,0);
close(barDistance(V(0,0),V(1,0),V(2,0),V(4,0)).gap,1);
close(barDistance(V(0,0),V(1,0),V(0,3),V(1,3)).gap,3);
const cases=[[40,1],[70,1.3],[90,1.5],[60,1.8],[60,1.2]];
receipts.flaps=cases.map(([psi,q])=>{
  const m=synthFlap(psi,q,70);assert.ok(m);
  const sim=simulate(m,720);
  assert.equal(sim.fails,0);
  const ph=sim.pts.map(p=>Math.atan2(p.B.y-m.B0.y,p.B.x-m.B0.x));
  const swing=(Math.max(...ph)-Math.min(...ph))/deg;
  close(swing,psi,.02);
  close(m.timeRatio,q,.011);
  const cr=closureReport(m,.5), sr=sweepReport(m,14);
  const fine=simulate(m,8192).pts;
  let minimum=Infinity;
  for(const p of fine)minimum=Math.min(minimum,barDistance(m.A0,p.A,m.B0,p.B).gap-14);
  assert.ok(sr.lower<=minimum+1e-8 && minimum<=sr.upper+.02);
  assert.ok(sr.upper-sr.lower<.17);
  return {swingDeg:psi,requestedRatio:q,sampledSwingDeg:swing,sampledRatio:m.timeRatio,toleranceLimitMm:cr.limit,marginAtHalfMm:cr.margin,clearanceMm:[sr.lower,sr.upper],state:sr.state};
});
const crossing=model(60,32,71,72,0,-1);
const sweep=sweepReport(crossing,43);
assert.equal(sweep.state,'overlap');
assert.ok(sweep.witness.gap<0);
const wide=sweepReport(near,100), bad=sweepReport(hidden,10);
assert.equal(wide.state,'overlap');
assert.equal(bad.state,'unresolved');
receipts.interference={widthMm:43,linkLengthsMm:[60,32,71,72],intervalMm:[sweep.lower,sweep.upper],witnessDeg:sweep.witness.angle/deg};
const stock={w:14,maxL:400};
const pins={hole:3};
const specimen=synthFlap(60,1.8,70);
const scales=scaleChoices(specimen,5,stock,pins);
assert.equal(scales[0],160);
const reference=fabricationModel(specimen,1.59);
assert.equal(closureReport(reference,5).state,'at limit');
assert.equal(closureReport(fabricationModel(specimen,1.60),5).state,'closes');
const stroke=strokeReport(reference);
close(stroke.swingDeg,59.9847637275,1e-8);
close(stroke.longToShortTimeRatio,1.80081466188,1e-8);
close(minimumTransmission(reference),22.9984207532,1e-8);
close(stroke.ccwRockerStrokeCrankDeg+stroke.cwRockerStrokeCrankDeg,360);
assert.equal(strokeReport(hidden),null);
assert.equal(minimumTransmission(hidden),null);
assert.equal(strokeReport(model(10,55,50,45)),null); // output rotates through 360 degrees
const reversed=strokeReport({...reference,s:-reference.s});
close(reversed.swingDeg,stroke.swingDeg);
close(reversed.cwRockerStrokeCrankDeg,stroke.cwRockerStrokeCrankDeg);
for(const [swing,ratio] of cases){
  const cut=fabricationModel(synthFlap(swing,ratio,70),1.23),out=strokeReport(cut);
  // Second route: resolve rocker direction continuously through 8,192 crank steps.
  const samples=simulate(cut,8192).pts.map(p=>Math.atan2(p.B.y-cut.B0.y,p.B.x-cut.B0.x));
  const unwrapped=[samples[0]];
  for(let i=1;i<samples.length;i++){
    const change=samples[i]-samples[i-1];
    unwrapped.push(unwrapped[i-1]+Math.atan2(Math.sin(change),Math.cos(change)));
  }
  const sampledSwing=(Math.max(...unwrapped)-Math.min(...unwrapped))/deg;
  close(sampledSwing,out.swingDeg,.001);
  const iMin=unwrapped.indexOf(Math.min(...unwrapped)),iMax=unwrapped.indexOf(Math.max(...unwrapped));
  const n=(iMax-iMin+samples.length)%samples.length;
  close(n*360/samples.length,out.ccwRockerStrokeCrankDeg,.065);
}
const capsule=plateOutline([V(5,0),V(20,0),V(5,0)],3);
assert(!capsule.includes('NaN') && !capsule.includes('Infinity'));
receipts.roundedFabrication={at159:closureReport(reference,5).state,at160:closureReport(fabricationModel(specimen,1.60),5).state,
  referenceSwingDeg:stroke.swingDeg,referenceTiming:stroke.longToShortTimeRatio,
  referenceMinTransmissionDeg:minimumTransmission(reference),relativeStopAngles:stroke.stops};
assert.equal(scaleChoices(specimen,5,{w:10,maxL:114,isStick:true},pins).length,0);
receipts.inverseBuild={swingDeg:60,ratio:1.8,toleranceMm:5,cardboardMinPercent:scales[0],
  cardboardMaxPercent:scales[scales.length-1],craftStickCount:0};
const m=synthFlap(70,1.3,50);
const gap=traceGap(m,simulate(m,720).pts.map(p=>p.P));
close(gap.mean,0,1e-9);
const curve=roundedCurve(model(100,30,75,55.4));
for(const [n,d] of [[0,1],[1,2],[1,1],[1,0],[-1,1],[3,7],[-5,13]]){
  const fibre=projectiveFibre(curve,n,d);
  assert.equal(fibre.count,2);
  const N=BigInt(n),D=BigInt(d);
  assert.equal(fibre.delta,curve.D0*D**4n+curve.D2*N*N*D*D+curve.D4*N**4n);
}
const tangent=roundedCurve(model(100,30,75,55));
assert.equal(projectiveFibre(tangent,1,0).count,1);
assert.equal(tangent.smoothQuartic,false);
assert.equal(projectiveFibre(roundedCurve(model(100,30,75,54.9)),1,0).count,0);
assert.equal(projectiveFibre(roundedCurve(model(80,80,90,90)),0,1).count,'all');
assert.equal(curve.smoothQuartic,true);
const reconstructed=roundedCurve(hidden);
assert.equal(projectiveFibre(reconstructed,1,0).count,1);
receipts.exactRoundedCurves={lengthUnit:'0.1 mm',halfTurn:[
  {couplerMm:75,rockerMm:55.4,positions:2},
  {couplerMm:75,rockerMm:55,positions:1},
  {couplerMm:75,rockerMm:54.9,positions:0}
],groundMm:100,crankMm:30};
// Direct count on P1 x P1, independently of square-root/character counting.
function countClosure(curve,p) {
  const P=BigInt(p), pts=Array.from({length:p},(_,i)=>[BigInt(i),1n]).concat([[1n,0n]]);
  let count=0, div=0n;
  function gcd(x,y){return y===0n?(x<0n?-x:x):gcd(y,x%y);}
  for (const c of [curve.C00,curve.C02,curve.C20,curve.C22,curve.C11]) div=gcd(div,c);
  for(const [n,d] of pts)for(const [x,y] of pts){
    const f=(curve.C02*d*d+curve.C22*n*n)*x*x+curve.C11*n*d*x*y+(curve.C00*d*d+curve.C20*n*n)*y*y;
    if(f/div%P===0n)count++;
  }
  return count;
}
const j=evenQuarticJ(curve.D0,curve.D2,curve.D4);
const twice=roundedCurve(model(200,60,150,110.8));
assert.deepEqual(evenQuarticJ(twice.D0,twice.D2,twice.D4),j);
assert.deepEqual(goodPrimeCounts(curve),goodPrimeCounts(twice));
assert.equal(evenQuarticJ(1n,0n,1n).numerator,'1728');
assert.equal(evenQuarticJ(1n,0n,1n).denominator,'1');
// t = -2,-1,1,2 have Legendre cross-ratio 9. The Legendre formula yields 4*73^3/81.
const cross=evenQuarticJ(4n,-5n,1n);
assert.equal(BigInt(cross.numerator)*81n,4n*73n**3n*BigInt(cross.denominator));
assert.equal(evenQuarticJ(tangent.D0,tangent.D2,tangent.D4),null);
for(const test of [curve,roundedCurve(model(75,25,60,50)),roundedCurve(model(90,30,65,56))]){
  const counts=goodPrimeCounts(test);assert.equal(counts.length,3);
  for(const row of counts)assert.equal(row.total,countClosure(test,row.prime));
}
receipts.exactRoundedCurves.j=j;
receipts.exactRoundedCurves.primeCounts=goodPrimeCounts(curve);
const output=process.argv[2];
if(output)fs.writeFileSync(output, JSON.stringify(receipts,null,2)+'\n');
console.log('All four original flap cases, the default, interval witnesses, clearance bounds, and the missed-angle regression passed.');
