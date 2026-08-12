#!/usr/bin/env python
"""Build viewer.html: interactive materials viewer for the N=50 converter.

Pan/zoom plan view, per-MATERIAL select/deselect, click-2-points distance
measurement in um, live coordinates, scale bar, and a process cross-section.
"""
import json

with open("out/masks.json") as f:
    masks = json.load(f)

chip = masks["chip"]
lay = {l["name"].split()[0]: l["features"] for l in masks["layers"]}

# derived material regions
X0, X1 = chip["x0"], chip["x1"]
si_field = [dict(t="rect", x=X0, y=chip["y0"], dx=X1 - X0, dy=chip["y1"] - chip["y0"])]
# released (air-clad) zone: band around the photonic circuit + MRR disc
air = [dict(t="rect", x=X0 + 100, y=-38, dx=X1 - X0 - 130, dy=50),
       dict(t="arc", c=[60.0, 17.7], r=28.0, a0=0, a1=360, w=56.0)]

layers = [
    dict(key="si",  name="Si substrate (carrier wafer)", color="#3a3f46",
         mode="fill", op=0.55, features=si_field,
         desc="⟨100⟩ silicon carrier; unetched field supports rails & pads"),
    dict(key="air", name="Released cavity — air cladding (XeF2 undercut)",
         color="#1c5d63", mode="fill", op=0.5, features=air,
         desc="XeF2-undercut zone: silica guides suspended, air above AND below"),
    dict(key="wg",  name="SiO2 waveguide core (thermal oxide, n=1.444, t=1.0 um)",
         color="#4da3ff", mode="stroke", op=1.0, features=lay["WG"],
         desc="the only optical material; mask SEO-WG-01, e-beam/DUV, full etch"),
    dict(key="met", name="Cu electrodes (Ti/TiN/Cu 200 nm/TiN, liftoff)",
         color="#c46b3a", mode="stroke", op=0.95, features=lay["MET"],
         desc="MRR fast EO pair + RF/DC pads; 50 weight pairs + rails; TiN barriers "
              "stop Cu drift into SiO2 under DC bias; mask SEO-MET-02"),
    dict(key="rel", name="Release windows (mask SEO-REL-03)",
         color="#e377c2", mode="fill", op=0.85, features=lay["REL"],
         desc="openings in resist for the XeF2 isotropic Si etch; 24 um pitch"),
]
DATA = dict(units="um", chip=chip, layers=layers)

html = """<!DOCTYPE html><html><head><meta charset="utf-8">
<title>N=50 silica converter — materials viewer</title>
<style>
body{margin:0;font-family:-apple-system,Helvetica,sans-serif;background:#0d0d0f;color:#ddd;overflow:hidden}
#side{position:fixed;top:0;left:0;bottom:0;width:295px;background:#17181c;
border-right:1px solid #2a2c33;padding:12px 14px;overflow-y:auto;box-sizing:border-box;z-index:5}
#side h1{font-size:14px;margin:0 0 4px}
#side .sub{font-size:11px;color:#888;margin-bottom:12px}
.mat{margin:7px 0;padding:7px 8px;background:#1e2026;border-radius:6px;cursor:pointer;
border:1px solid transparent}
.mat:hover{border-color:#3a3d46}
.mat .row{display:flex;align-items:center;gap:8px;font-size:12.5px}
.mat .sw{width:14px;height:14px;border-radius:3px;flex:none}
.mat .desc{font-size:10.5px;color:#9aa;margin:4px 0 0 22px;line-height:1.35}
.mat.off{opacity:0.38}
#xsec{margin-top:14px;border-top:1px solid #2a2c33;padding-top:10px}
#xsec h2{font-size:12px;margin:0 0 6px;color:#bbb}
#svgwrap{position:fixed;top:0;left:295px;right:0;bottom:0;overflow:hidden;cursor:grab}
svg{display:block}
#hud{position:fixed;top:10px;right:14px;background:#000b;padding:6px 10px;border-radius:5px;
font-family:ui-monospace,monospace;font-size:12px;color:#8fd;z-index:6}
#scalebar{position:fixed;bottom:14px;left:315px;background:#000b;padding:6px 10px;
border-radius:4px;font-family:ui-monospace,monospace;font-size:12px;z-index:6}
#sbline{height:4px;background:#fff;margin-top:3px}
#measure{position:fixed;bottom:14px;right:16px;background:#000b;padding:6px 10px;
border-radius:4px;font-family:ui-monospace,monospace;font-size:13px;color:#ff6;z-index:6}
#hint{position:fixed;top:10px;left:315px;font-size:11px;color:#777;z-index:6}
</style></head><body>
<div id="side">
<h1>N=50 silica EO&rarr;photonic converter</h1>
<div class="sub">materials view &mdash; click a card to show/hide that material</div>
<div id="mats"></div>
<div id="xsec"><h2>process cross-section (A&ndash;A&prime; through a suspended guide)</h2>
<svg viewBox="0 0 260 120" style="width:100%;background:#101014;border-radius:6px">
<rect x="0" y="78" width="260" height="42" fill="#3a3f46"/>
<rect x="30" y="52" width="200" height="26" fill="#101014"/>
<rect x="112" y="52" width="36" height="14" fill="#4da3ff"/>
<rect x="60" y="30" width="34" height="14" fill="#3a3f46"/>
<rect x="60" y="16" width="34" height="14" fill="#c46b3a"/>
<rect x="170" y="30" width="34" height="14" fill="#3a3f46"/>
<rect x="170" y="16" width="34" height="14" fill="#c46b3a"/>
<text x="130" y="46" fill="#4da3ff" font-size="8" text-anchor="middle">SiO2 core 0.991 x 1.0 um</text>
<text x="130" y="74" fill="#1c9da8" font-size="8" text-anchor="middle">air (XeF2 cavity)</text>
<text x="77" y="12" fill="#c46b3a" font-size="8" text-anchor="middle">Cu electrode</text>
<text x="187" y="12" fill="#c46b3a" font-size="8" text-anchor="middle">Cu electrode</text>
<text x="24" y="112" fill="#889" font-size="8">Si substrate</text>
<text x="77" y="40" fill="#889" font-size="7" text-anchor="middle">Si mesa</text>
</svg>
<div class="desc" style="font-size:10.5px;color:#9aa;margin-top:6px">
Electrodes sit on unreleased Si mesas flanking the cavity; the field E crosses
the suspended silica core (DC Kerr). Membrane held by 24 um-pitch tether ribs.</div>
</div>
</div>
<div id="hint">wheel = zoom &middot; drag = pan &middot; click 2 pts = measure &middot; M = clear</div>
<div id="hud"></div>
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
const gRoot = el('g',{});svg.appendChild(gRoot);
DATA.layers.forEach((lay,i)=>{
  const g = el('g',{id:'lay'+i,opacity:lay.op});
  lay.features.forEach(ft=>{
    if(ft.t==='line'){
      const pts = ft.p.map(p=>`${p[0]},${-p[1]}`).join(' ');
      g.appendChild(el('polyline',{points:pts,fill:'none',stroke:lay.color,
        'stroke-width':ft.w,'stroke-linecap':'butt'}));
    } else if(ft.t==='arc'){
      if(ft.a1-ft.a0>=360){
        g.appendChild(el('circle',{cx:ft.c[0],cy:-ft.c[1],r:ft.r,
          fill: lay.mode==='fill'?lay.color:'none',
          stroke:lay.mode==='fill'?'none':lay.color,'stroke-width':ft.w||0}));
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
        fill:lay.color}));
    }
  });
  gRoot.appendChild(g);
  const card=document.createElement('div');card.className='mat';card.id='card'+i;
  card.innerHTML=`<div class="row"><div class="sw" style="background:${lay.color}"></div>
    <b>${lay.name}</b></div><div class="desc">${lay.desc}</div>`;
  card.onclick=()=>{const g=document.getElementById('lay'+i);
    const off=g.style.display==='none';
    g.style.display=off?'':'none';card.classList.toggle('off',!off);};
  document.getElementById('mats').appendChild(card);
});
const gM = el('g',{});svg.appendChild(gM);
function toUser(evt){const pt=svg.createSVGPoint();pt.x=evt.clientX;pt.y=evt.clientY;
  const p=pt.matrixTransform(svg.getScreenCTM().inverse());return {x:p.x,y:-p.y};}
wrap.addEventListener('wheel',e=>{e.preventDefault();
  const p=toUser(e);const f=e.deltaY>0?1.15:1/1.15;
  vb.x=p.x-(p.x-vb.x)*f; vb.y=(-p.y)-((-p.y)-vb.y)*f; vb.w*=f; vb.h*=f; setVB();},{passive:false});
let drag=null,moved=false;
wrap.addEventListener('mousedown',e=>{drag={x:e.clientX,y:e.clientY,vx:vb.x,vy:vb.y};
  wrap.style.cursor='grabbing';moved=false;});
window.addEventListener('mousemove',e=>{
  const u=toUser(e);
  document.getElementById('hud').textContent=`x=${u.x.toFixed(2)}  y=${u.y.toFixed(2)} um`;
  if(!drag)return;
  const sc=vb.w/wrap.clientWidth;
  const dx=(e.clientX-drag.x)*sc, dy=(e.clientY-drag.y)*sc;
  if(Math.abs(e.clientX-drag.x)+Math.abs(e.clientY-drag.y)>3)moved=true;
  vb.x=drag.vx-dx; vb.y=drag.vy-dy; setVB();});
window.addEventListener('mouseup',e=>{wrap.style.cursor='grab';
  if(drag&&!moved&&e.target.closest('#svgwrap'))clickMeasure(e); drag=null;});
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
  document.getElementById('sbtext').textContent=nice>=1000?`${nice/1000} mm`:`${nice} um`;
}
window.addEventListener('resize',setVB);
setVB();
</script></body></html>"""
html = html.replace("__DATA__", json.dumps(DATA))
with open("viewer.html", "w") as f:
    f.write(html)
print("viewer.html written,", len(html) // 1024, "KB")
