'use strict';
(() => {
  const $=id=>document.getElementById(id), canvas=$('model'), controls=['theta','spread','play','top','home'];
  const focusButtons=Array.from(document.querySelectorAll('.focusJoint'));
  let focusJoint='all';
  const info={all:'Choose a joint to see its washers and keepers up close.',
    O:'O anchors the crank. Two collars clamp the board. The 3 mm spacer and lower washer support the crank. Set endplay under the top keeper.',
    D:'D anchors the rocker. The rocker and pointer swing on this shaft. Set the lower support stack to the same height as the crank at O.',
    A:'A joins the crank to the coupler. The 10 mm spacer separates the levels. The finger grip is bonded to the extended shaft above the top keeper.',
    B:'B joins the rocker to the coupler. A 10 mm spacer separates the levels. The lower collar and shaft tail pass above the base during the turn.'};
  function setFocus(joint){focusJoint=joint;for(const button of focusButtons)button.setAttribute('aria-pressed',String(button.dataset.joint===joint));$('jointInfo').textContent=info[joint];}
  function fallback(){for(const button of focusButtons)button.disabled=true;canvas.hidden=true;canvas.style.display='none';$('pins').hidden=true;$('fallback').hidden=false;for(const k of controls)$(k).disabled=true;}
  const gl=canvas.getContext('webgl',{alpha:false,antialias:true});
  if(!gl){fallback();return;}
  const data=window.LinkwrightAssembly, {ground:g,crank:a,coupler:b,rocker:c}=data.lengths;
  const D=[g,0], O=[0,0], A0=[a,0], B0=data.B0, rad=Math.PI/180;
  const angle=(first,last)=>Math.atan2(last[1]-first[1],last[0]-first[0]);
  const phi0=angle(D,B0),psi0=angle(A0,B0);
  const vSource=`
    attribute vec3 pos,normal;
    uniform vec2 oldPin,newPin,rotation;uniform float lift,halfHeight,aspect;
    uniform vec3 right,up,eye,look;
    varying float diffuse;
    vec2 turn(vec2 v){return vec2(rotation.x*v.x-rotation.y*v.y,rotation.y*v.x+rotation.x*v.y);}
    void main(){
      vec3 p=vec3(newPin+turn(pos.xy-oldPin),pos.z+lift)-look;
      gl_Position=vec4(dot(p,right)/(halfHeight*aspect),dot(p,up)/halfHeight,-dot(p,eye)/430.,1.);
      vec3 n=normalize(vec3(turn(normal.xy),normal.z));
      diffuse=.62+.38*max(dot(n,normalize(vec3(-.45,-.45,1.))),0.);
    }`;
  const fSource='precision mediump float;uniform vec3 color;varying float diffuse;void main(){gl_FragColor=vec4(color*diffuse,1.);}';
  function shader(type,source){const sh=gl.createShader(type);gl.shaderSource(sh,source);gl.compileShader(sh);if(!gl.getShaderParameter(sh,gl.COMPILE_STATUS))throw Error(gl.getShaderInfoLog(sh));return sh;}
  const p=gl.createProgram();
  try{
    gl.attachShader(p,shader(gl.VERTEX_SHADER,vSource));gl.attachShader(p,shader(gl.FRAGMENT_SHADER,fSource));gl.linkProgram(p);
    if(!gl.getProgramParameter(p,gl.LINK_STATUS))throw Error('Graphics program failed to link.');
  }catch{fallback();return;}
  gl.useProgram(p);
  const uv={};for(const n of ['oldPin','newPin','rotation','lift','halfHeight','aspect','right','up','eye','look','color'])uv[n]=gl.getUniformLocation(p,n);
  const locPos=gl.getAttribLocation(p,'pos'),locNorm=gl.getAttribLocation(p,'normal');
  function unpack(src,signed){const raw=atob(src),bytes=Uint8Array.from(raw,ch=>ch.charCodeAt(0)),v=new DataView(bytes.buffer),out=[];
    for(let i=0;i<bytes.length;i+=2)out.push(signed?v.getInt16(i,true)/100:v.getUint16(i,true));return out;}
  function cross(u,v){return [u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]];}
  const minus=(u,v)=>u.map((n,i)=>n-v[i]), dot=(u,v)=>u.reduce((sum,n,i)=>sum+n*v[i],0);
  const buffers=data.parts.map(part=>{
    const pos=unpack(part.xyz,true),idx=unpack(part.triangles,false),flat=[];
    for(let i=0;i<idx.length;i+=3){const face=idx.slice(i,i+3).map(j=>pos.slice(3*j,3*j+3)),normal=cross(minus(face[1],face[0]),minus(face[2],face[0]));
      const length=Math.hypot(...normal);if(length<1e-8)continue;const n=normal.map(k=>k/length);
      for(const v of face)flat.push(...v,...n);}
    const buffer=gl.createBuffer();gl.bindBuffer(gl.ARRAY_BUFFER,buffer);gl.bufferData(gl.ARRAY_BUFFER,new Float32Array(flat),gl.STATIC_DRAW);
    return {...part,buffer,count:flat.length/6};
  });
  gl.clearColor(.93,.95,.93,1);gl.enable(gl.DEPTH_TEST);gl.enableVertexAttribArray(locPos);gl.enableVertexAttribArray(locNorm);
  let az=-1.23,elev=.77,drag=null,animation=null,startTime=0,initialTheta=0;
  function pose(theta){
    const A=[a*Math.cos(theta),a*Math.sin(theta)],delta=minus(D,A),d=Math.hypot(...delta),v=delta.map(z=>z/d);
    const along=(b*b-c*c+d*d)/(2*d),h=Math.sqrt(Math.max(0,b*b-along*along));
    return {A,B:[A[0]+v[0]*along-v[1]*h,A[1]+v[1]*along+v[0]*h]};
  }
  function draw(){
    const w=canvas.clientWidth,h=canvas.clientHeight,ratio=Math.min(window.devicePixelRatio||1,2);
    if(!w||!h)return;const width=Math.round(w*ratio),height=Math.round(h*ratio);
    if(width!==canvas.width||height!==canvas.height){canvas.width=width;canvas.height=height;}
    gl.viewport(0,0,width,height);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT);
    const theta=+$('theta').value*rad,spread=+$('spread').value/100,{A,B}=pose(theta);
    const right=[-Math.sin(az),Math.cos(az),0],eye=[Math.cos(az)*Math.cos(elev),Math.sin(az)*Math.cos(elev),Math.sin(elev)],up=cross(eye,right);
    const joints={O,D,A,B},fixedFocus=focusJoint==='O'||focusJoint==='D';
    const look=focusJoint==='all'?[25,41,8+22*spread]:[...joints[focusJoint],(fixedFocus?5:22)+(fixedFocus?20:43)*spread];
    const frame=focusJoint==='all'?132+12*spread:(38+34*spread);
    const half=Math.max(frame,frame/(w/h));
    for(const [n,val] of Object.entries({right,eye,up,look}))gl.uniform3fv(uv[n],val);
    gl.uniform1f(uv.aspect,w/h);gl.uniform1f(uv.halfHeight,half);
    const move={fixed:[O,O,0],A:[A0,A,0],B:[B0,B,0],crank:[O,O,theta],coupler:[A0,A,angle(A,B)-psi0],rocker:[D,D,angle(D,B)-phi0]};
    for(const part of buffers){const [oldPin,newPin,rot]=move[part.motion];
      gl.uniform2fv(uv.oldPin,oldPin);gl.uniform2fv(uv.newPin,newPin);gl.uniform2f(uv.rotation,Math.cos(rot),Math.sin(rot));
      gl.uniform1f(uv.lift,part.spread*spread);gl.uniform3fv(uv.color,data.colors[part.category]);
      gl.bindBuffer(gl.ARRAY_BUFFER,part.buffer);gl.vertexAttribPointer(locPos,3,gl.FLOAT,false,24,0);gl.vertexAttribPointer(locNorm,3,gl.FLOAT,false,24,12);
      gl.drawArrays(gl.TRIANGLES,0,part.count);}
    [[O,18+48*spread],[D,18+48*spread],[A,47+76*spread],[B,32+68*spread]].forEach(([pt,z],i)=>{
      const p=minus([...pt,z],look),tag=$('pins').children[i];
      tag.hidden=focusJoint!=='all'&&'ODAB'[i]!==focusJoint;tag.style.left=(w/2+dot(p,right)*h/(2*half))+'px';tag.style.top=(h/2-dot(p,up)*h/(2*half))+'px';
    });
    $('thetaValue').value=Math.round(theta/rad)+'°';$('spreadValue').value=spread?Math.round(100*spread)+'% separated':'Assembled';}
  function stop(){if(animation!==null)cancelAnimationFrame(animation);animation=null;$('play').textContent='Play one turn';}
  function frame(time){if(!startTime)startTime=time;const delta=(time-startTime)*.045;
    $('theta').value=delta>=360?initialTheta:(initialTheta+delta)%360;draw();
    if(delta>=360)stop();else animation=requestAnimationFrame(frame);}
  $('play').onclick=()=>{if(animation!==null){stop();return;}initialTheta=+$('theta').value;startTime=0;$('play').textContent='Pause turn';animation=requestAnimationFrame(frame);};
  $('theta').oninput=()=>{stop();draw();};$('spread').oninput=draw;
  for(const button of focusButtons)button.onclick=()=>{setFocus(button.dataset.joint);draw();};
  $('top').onclick=()=>{az=-Math.PI/2;elev=Math.PI/2;draw();};
  $('home').onclick=()=>{stop();setFocus('all');az=-1.23;elev=.77;$('theta').value=0;$('spread').value=0;draw();};
  function rotate(dx,dy){az-=dx*.008;elev=Math.max(-1.42,Math.min(Math.PI/2,elev+dy*.008));draw();}
  canvas.onpointerdown=e=>{drag=[e.clientX,e.clientY];canvas.setPointerCapture(e.pointerId);};
  canvas.onpointermove=e=>{if(!drag)return;rotate(e.clientX-drag[0],e.clientY-drag[1]);drag=[e.clientX,e.clientY];};
  canvas.onpointerup=canvas.onpointercancel=()=>{drag=null;};
  canvas.onkeydown=e=>{const d={ArrowLeft:[-12,0],ArrowRight:[12,0],ArrowUp:[0,12],ArrowDown:[0,-12]}[e.key];if(d){e.preventDefault();rotate(...d);}};
  canvas.addEventListener('webglcontextlost',e=>{e.preventDefault();stop();fallback();});
  new ResizeObserver(draw).observe(canvas);draw();
})();
