#!/usr/bin/env python
"""Build viewer.html: interactive pan/zoom mask viewer with um measurements."""
import json

with open("out/masks.json") as f:
    masks = json.load(f)

html = """<!DOCTYPE html><html><head><meta charset="utf-8">
<title>N=50 silica converter — mask viewer</title>
<style>
body{margin:0;font-family:-apple-system,Helvetica,sans-serif;background:#111;color:#ddd}
#bar{padding:8px 14px;background:#1b1b1b;display:flex;gap:18px;align-items:center;
flex-wrap:wrap;border-bottom:1px solid #333}
#bar label{cursor:pointer;font-size:13px}
#bar .hint{color:#888;font-size:12px}
#readout{font-family:ui-monospace,monospace;font-size:13px;color:#8fd}
#svgwrap{width:100vw;height:calc(100vh - 46px);overflow:hidden;cursor:grab}
svg{display:block}
#scalebar{position:fixed;bottom:14px;left:16px;background:#000a;padding:6px 10px;
border-radius:4px;font-family:ui-monospace,monospace;font-size:12px}
#sbline{height:4px;background:#fff;margin-top:3px}
#measure{position:fixed;bottom:14px;right:16px;background:#000a;padding:6px 10px;
border-radius:4px;font-family:ui-monospace,monospace;font-size:13px;color:#ff6}
</style></head><body>
<div id="bar">
<b>N=50 silica EO&rarr;photonic converter — photomasks</b>
<span id="layerboxes"></span>
<span class="hint">wheel = zoom &middot; drag = pan &middot; click 2 points = measure &middot; M = clear measurement</span>
<span id="readout"></span>
</div>
<div id="svgwrap"><svg id="svg" xmlns="http://www.w3.org/2000/svg"></svg></div>
<div id="scalebar"><span id="sbtext"></span><div id="sbline"></div></div>
<div id="measure">click two points to measure</div>
<script>
const DATA = __DATA__;
const svg = document.getElementById('svg');
const wrap = document.getElementById('svgwrap');
const chip = DATA.chip;
let vb = {x: chip.x0, y: -chip.y1, w: chip.x1-chip.x0, h: chip.y1-chip.y0};
function setVB(){ svg.setAttribute('viewBox', `${vb.x} ${vb.y} ${vb.w} ${vb.h}`);
  svg.setAttribute('width', wrap.clientWidth); svg.setAttribute('height', wrap.clientHeight);
  updateScalebar(); }
function el(n,a){const e=document.createElementNS('http://www.w3.org/2000/svg',n);
  for(const k in a)e.setAttribute(k,a[k]);return e;}
// note: y flipped (SVG y down) -> plot at -y
const gRoot = el('g',{});svg.appendChild(gRoot);
DATA.layers.forEach((lay,i)=>{
  const g = el('g',{id:'lay'+i});
  lay.features.forEach(ft=>{
    if(ft.t==='line'){
      const pts = ft.p.map(p=>`${p[0]},${-p[1]}`).join(' ');
      g.appendChild(el('polyline',{points:pts,fill:'none',stroke:lay.color,
        'stroke-width':ft.w,'stroke-linecap':'butt'}));
    } else if(ft.t==='arc'){
      if(ft.a1-ft.a0>=360){
        g.appendChild(el('circle',{cx:ft.c[0],cy:-ft.c[1],r:ft.r,fill:'none',
          stroke:lay.color,'stroke-width':ft.w}));
      } else {
        const a0=ft.a0*Math.PI/180, a1=ft.a1*Math.PI/180;
        const x0=ft.c[0]+ft.r*Math.cos(a0), y0=-(ft.c[1]+ft.r*Math.sin(a0));
        const x1=ft.c[0]+ft.r*Math.cos(a1), y1=-(ft.c[1]+ft.r*Math.sin(a1));
        const large=(ft.a1-ft.a0)>180?1:0;
        g.appendChild(el('path',{d:`M ${x0} ${y0} A ${ft.r} ${ft.r} 0 ${large} 0 ${x1} ${y1}`,
          fill:'none',stroke:lay.color,'stroke-width':ft.w}));
      }
    } else {
      g.appendChild(el('rect',{x:ft.x,y:-(ft.y+ft.dy),width:ft.dx,height:ft.dy,
        fill:lay.color,opacity:0.85}));
    }
  });
  gRoot.appendChild(g);
  const cb=document.createElement('label');
  cb.innerHTML=`<input type="checkbox" checked onchange="document.getElementById('lay${i}').style.display=this.checked?'':'none'"> <span style="color:${lay.color}">${lay.name}</span>`;
  document.getElementById('layerboxes').appendChild(cb);
});
const gM = el('g',{});svg.appendChild(gM);
function toUser(evt){const pt=svg.createSVGPoint();pt.x=evt.clientX;pt.y=evt.clientY;
  const p=pt.matrixTransform(svg.getScreenCTM().inverse());return {x:p.x,y:-p.y};}
// zoom
wrap.addEventListener('wheel',e=>{e.preventDefault();
  const p=toUser(e);const f=e.deltaY>0?1.15:1/1.15;
  vb.x=p.x-(p.x-vb.x)*f; vb.y=(-p.y)-((-p.y)-vb.y)*f; vb.w*=f; vb.h*=f; setVB();},{passive:false});
// pan
let drag=null;
wrap.addEventListener('mousedown',e=>{drag={x:e.clientX,y:e.clientY,vx:vb.x,vy:vb.y};
  wrap.style.cursor='grabbing';moved=false;});
let moved=false;
window.addEventListener('mousemove',e=>{
  const u=toUser(e);
  document.getElementById('readout').textContent=`x=${u.x.toFixed(2)} um  y=${u.y.toFixed(2)} um`;
  if(!drag)return;
  const sc=vb.w/wrap.clientWidth;
  const dx=(e.clientX-drag.x)*sc, dy=(e.clientY-drag.y)*sc;
  if(Math.abs(dx)+Math.abs(dy)>0.5*sc*3)moved=true;
  vb.x=drag.vx-dx; vb.y=drag.vy-dy; setVB();});
window.addEventListener('mouseup',e=>{wrap.style.cursor='grab';
  if(drag&&!moved)clickMeasure(e); drag=null;});
// measure
let m1=null;
function clickMeasure(e){
  const u=toUser(e);
  if(!m1){m1=u; gM.innerHTML='';
    gM.appendChild(el('circle',{cx:u.x,cy:-u.y,r:vb.w/300,fill:'#ff6'}));
    document.getElementById('measure').textContent=`P1 = (${u.x.toFixed(2)}, ${u.y.toFixed(2)}) — click P2`;
  } else {
    const d=Math.hypot(u.x-m1.x,u.y-m1.y);
    gM.appendChild(el('line',{x1:m1.x,y1:-m1.y,x2:u.x,y2:-u.y,stroke:'#ff6',
      'stroke-width':vb.w/600,'stroke-dasharray':`${vb.w/200} ${vb.w/400}`}));
    gM.appendChild(el('circle',{cx:u.x,cy:-u.y,r:vb.w/300,fill:'#ff6'}));
    const t=el('text',{x:(m1.x+u.x)/2,y:-(m1.y+u.y)/2-vb.w/150,fill:'#ff6',
      'font-size':vb.w/60,'text-anchor':'middle'});
    t.textContent=`${d.toFixed(3)} um`;
    gM.appendChild(t);
    document.getElementById('measure').textContent=
      `d = ${d.toFixed(3)} um   (dx=${(u.x-m1.x).toFixed(3)}, dy=${(u.y-m1.y).toFixed(3)})`;
    m1=null;
  }
}
window.addEventListener('keydown',e=>{if(e.key==='m'||e.key==='M'){gM.innerHTML='';m1=null;
  document.getElementById('measure').textContent='click two points to measure';}});
function updateScalebar(){
  const target=vb.w/6; const pow=Math.pow(10,Math.floor(Math.log10(target)));
  let nice=pow; for(const m of [1,2,5,10]) if(m*pow<=target) nice=m*pow;
  const px=nice/vb.w*wrap.clientWidth;
  document.getElementById('sbline').style.width=px+'px';
  document.getElementById('sbtext').textContent=
    nice>=1000?`${nice/1000} mm`:`${nice} um`;
}
window.addEventListener('resize',setVB);
setVB();
</script></body></html>"""
html = html.replace("__DATA__", json.dumps(masks))
with open("viewer.html", "w") as f:
    f.write(html)
print("viewer.html written,", len(html) // 1024, "KB")
