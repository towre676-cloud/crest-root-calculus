const vm=require('node:vm'),fs=require('node:fs'),assert=require('node:assert/strict');
const path=__dirname+require('node:path').sep;
const ctx={},items={},uni={};let calls=[];
const webgl=new Proxy({}, {get(o,k){
    if(k==='getShaderParameter'||k==='getProgramParameter')return ()=>true;
    if(k==='getUniformLocation')return (_,n)=>n;
    if(k==='getAttribLocation')return ()=>1;
    if(k==='uniform2fv'||k==='uniform3fv'||k==='uniform1f')return (n,v)=>{uni[n]=v};
    if(k==='uniform2f')return (n,a,b)=>{uni[n]=[a,b]};
    if(k==='drawArrays')return ()=>calls.push({...uni});
    return ()=>({});
}});
function elem(id){
 if(!items[id]) items[id]={value:0,clientWidth:640,clientHeight:420,getContext:()=>webgl,style:{},children:Array.from({length:4},()=>({style:{}})),addEventListener(){}};
 return items[id];
}
Object.assign(ctx,{document:{getElementById:elem,querySelectorAll:()=>[]},window:{devicePixelRatio:1},atob,ResizeObserver:class{observe(){}},Math,Float32Array,DataView,Uint8Array,cancelAnimationFrame(){}});
vm.createContext(ctx);
vm.runInContext(fs.readFileSync(path+'mesh.js','utf8'),ctx);
vm.runInContext(fs.readFileSync(path+'view.js','utf8'),ctx);
const model=ctx.window.LinkwrightAssembly;
const near=(actual,expected,tol)=>assert.ok(Math.abs(actual-expected)<tol,actual+' vs '+expected);
const g=93.8,a=46.6,b=84.1,c=111.3;
const t0B=model.B0;
for(let deg=0;deg<360;deg+=1){
  elem('theta').value=deg;calls=[];elem('theta').oninput();
  const t=deg*Math.PI/180,A=[a*Math.cos(t),a*Math.sin(t)];
  const d=Math.hypot(A[0]-g,A[1]), ray=Math.atan2(A[1],A[0]-g)-Math.acos((d*d+c*c-b*b)/(2*d*c));
  const B=[g+c*Math.cos(ray),c*Math.sin(ray)];
  assert.equal(calls.length,39);
  for(let i=0;i<39;i++){
    const name=model.parts[i].name,draw=calls[i];
    const move=p=>{const [px,py]=p.map((v,j)=>v-draw.oldPin[j]),[co,si]=draw.rotation;
      return [draw.newPin[0]+px*co-py*si,draw.newPin[1]+px*si+py*co]};
    const test=(p,q)=>move(p).forEach((v,j)=>near(v,q[j],1e-8));
    if(name==='crank'){test([0,0],[0,0]);test([a,0],A);}
    else if(name==='coupler'){test([a,0],A);test(t0B,B);}
    else if(name==='rocker'){test([g,0],[g,0]);test(t0B,B);}
    else if(name.startsWith('A_'))test([a,0],A);
    else if(name.startsWith('B_'))test(t0B,B);
    else test([0,0],[0,0]);
    near(draw.lift,0,1e-10);
  }
}
console.log('Actual view.js WebGL draw uniforms: 14,040 physical transforms across 360 crank settings passed independent rocker-circle oracle.');
