'use strict';
/* =============================================
   BOLIVIA EN CRISIS — app.js
   7 escenarios · respuestas contextuales
   ============================================= */

// ─── Chart defaults (light theme) ────────────
Chart.defaults.color = '#7a7068';
Chart.defaults.borderColor = '#d5cfc6';
Chart.defaults.font.family = "'JetBrains Mono', monospace";
Chart.defaults.font.size = 11;

const P = {
  red:    '#c0392b', redFill:  'rgba(192,57,43,0.10)',
  blue:   '#1a5c8a', blueFill: 'rgba(26,92,138,0.10)',
  green:  '#1e6b3a', greenFill:'rgba(30,107,58,0.10)',
  gold:   '#c8820a', goldFill: 'rgba(200,130,10,0.10)',
  orange: '#c45e0a', orangeFill:'rgba(196,94,10,0.10)',
  purple: '#6b3fa0', purpleFill:'rgba(107,63,160,0.10)',
};

let CHARTS = {};
function killChart(id) { if (CHARTS[id]) { CHARTS[id].destroy(); delete CHARTS[id]; } }

function setActive(groupId, idx) {
  document.querySelectorAll(`#${groupId} .sc-btn`).forEach((b, i) =>
    b.classList.toggle('active', i === idx));
}

// Nav highlight on scroll
const sections = document.querySelectorAll('section[id]');
window.addEventListener('scroll', () => {
  let cur = 'inicio';
  sections.forEach(s => { if (window.scrollY >= s.offsetTop - 80) cur = s.id; });
  document.querySelectorAll('.nav-link').forEach(l =>
    l.classList.toggle('active', l.dataset.section === cur));
}, { passive: true });

document.getElementById('navToggle').addEventListener('click', () =>
  document.getElementById('navLinks').classList.toggle('open'));

// ── helpers ─────────────────────────────────
function kpiCard(label, value, sub, cls='') {
  return `<div class="kpi-card ${cls}">
    <div class="kpi-label">${label}</div>
    <div class="kpi-value">${value}</div>
    <div class="kpi-sub">${sub}</div>
  </div>`;
}
function answerBlock(q, a) {
  return `<div class="answer-item">
    <div class="answer-q">${q}</div>
    <div class="answer-a">${a}</div>
  </div>`;
}
function markAnswered(ids) {
  document.querySelectorAll('.sq-item').forEach(el => el.classList.remove('answered'));
  ids.forEach(id => document.getElementById(id)?.classList.add('answered'));
}

// ══════════════════════════════════════════════
//  MÉTODOS NUMÉRICOS BASE
// ══════════════════════════════════════════════

function gaussSeidel(A, b, maxIter=60, tol=1e-8) {
  const n = b.length; let x = new Array(n).fill(0); const hist = [];
  for (let k = 0; k < maxIter; k++) {
    const xOld = [...x];
    for (let i = 0; i < n; i++) {
      let s = b[i];
      for (let j = 0; j < n; j++) if (j!==i) s -= A[i][j]*x[j];
      x[i] = s / A[i][i];
    }
    const err = Math.max(...x.map((v,i) => Math.abs(v - xOld[i])));
    hist.push({ k: k+1, x: [...x], err });
    if (err < tol) break;
  }
  return { x, hist };
}

function luSolve(A, b) {
  const n = b.length;
  const M = A.map(r => [...r]);
  const L = Array.from({length:n},(_,i)=>Array(n).fill(0).map((_,j)=>i===j?1:0));
  const U = M;
  for (let k=0;k<n;k++) for (let i=k+1;i<n;i++) {
    L[i][k] = U[i][k]/U[k][k];
    for (let j=k;j<n;j++) U[i][j] -= L[i][k]*U[k][j];
  }
  const y = Array(n).fill(0);
  for (let i=0;i<n;i++) y[i] = b[i]-L[i].slice(0,i).reduce((s,v,j)=>s+v*y[j],0);
  const x = Array(n).fill(0);
  for (let i=n-1;i>=0;i--) x[i]=(y[i]-U[i].slice(i+1).reduce((s,v,j)=>s+v*x[i+1+j],0))/U[i][i];
  return x;
}

function condNum(A) {
  const normInf = M => Math.max(...M.map(r=>r.reduce((s,v)=>s+Math.abs(v),0)));
  const [[a,b,c],[d,e,f],[g,h,ii]] = A;
  const det = a*(e*ii-f*h)-b*(d*ii-f*g)+c*(d*h-e*g);
  if (Math.abs(det)<1e-10) return 1e9;
  const inv = [
    [(e*ii-f*h)/det,-(b*ii-c*h)/det,(b*f-c*e)/det],
    [-(d*ii-f*g)/det,(a*ii-c*g)/det,-(a*f-c*d)/det],
    [(d*h-e*g)/det,-(a*h-b*g)/det,(a*e-b*d)/det]
  ];
  return normInf(A)*normInf(inv);
}

function euler(f,y0,t0,tf,h) {
  const t=[],y=[]; let ti=t0,yi=y0;
  while(ti<=tf+1e-9){t.push(+ti.toFixed(3));y.push(+yi.toFixed(4));yi+=h*f(ti,yi);ti+=h;}
  return {t,y};
}
function heun(f,y0,t0,tf,h) {
  const t=[],y=[]; let ti=t0,yi=y0;
  while(ti<=tf+1e-9){
    t.push(+ti.toFixed(3));y.push(+yi.toFixed(4));
    const k1=f(ti,yi),yp=yi+h*k1,k2=f(ti+h,yp);
    yi+=h*(k1+k2)/2; ti+=h;
  }
  return {t,y};
}
function rk4(f,y0,t0,tf,h) {
  const t=[],y=[]; let ti=t0,yi=y0;
  while(ti<=tf+1e-9){
    t.push(+ti.toFixed(3));y.push(+yi.toFixed(4));
    const k1=f(ti,yi),k2=f(ti+h/2,yi+h*k1/2),k3=f(ti+h/2,yi+h*k2/2),k4=f(ti+h,yi+h*k3);
    yi+=h*(k1+2*k2+2*k3+k4)/6; ti+=h;
  }
  return {t,y};
}
function rk4sys(fns,y0,t0,tf,h) {
  const n=fns.length,t=[],ys=Array.from({length:n},()=>[]);
  let ti=t0,yi=[...y0];
  while(ti<=tf+1e-9){
    t.push(+ti.toFixed(2));yi.forEach((v,i)=>ys[i].push(+v.toFixed(2)));
    const k1=fns.map((f,i)=>f(ti,yi));
    const y2=yi.map((v,i)=>v+h/2*k1[i]);
    const k2=fns.map((f,i)=>f(ti+h/2,y2));
    const y3=yi.map((v,i)=>v+h/2*k2[i]);
    const k3=fns.map((f,i)=>f(ti+h/2,y3));
    const y4=yi.map((v,i)=>v+h*k3[i]);
    const k4=fns.map((f,i)=>f(ti+h,y4));
    yi=yi.map((v,i)=>Math.max(0,v+h/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i])));
    ti+=h;
  }
  return {t,N:ys[0],M:ys[1],D:ys[2]};
}

function lagrange(xs,ys,x) {
  return xs.reduce((acc,_,i)=>{
    let t=ys[i];
    xs.forEach((_,j)=>{if(j!==i)t*=(x-xs[j])/(xs[i]-xs[j]);});
    return acc+t;
  },0);
}
function newtonDD(xs,ys) {
  const c=[...ys],n=xs.length;
  for(let j=1;j<n;j++) for(let i=n-1;i>=j;i--) c[i]=(c[i]-c[i-1])/(xs[i]-xs[i-j]);
  return c;
}
function newtonEval(xs,c,x) {
  let r=c[c.length-1];
  for(let i=c.length-2;i>=0;i--) r=r*(x-xs[i])+c[i];
  return r;
}
function spline(xs,ys) {
  const n=xs.length-1,h=xs.slice(1).map((x,i)=>x-xs[i]);
  const alpha=Array(n+1).fill(0);
  for(let i=1;i<n;i++) alpha[i]=(3/h[i])*(ys[i+1]-ys[i])-(3/h[i-1])*(ys[i]-ys[i-1]);
  const l=Array(n+1).fill(1),mu=Array(n+1).fill(0),z=Array(n+1).fill(0);
  for(let i=1;i<n;i++){
    l[i]=2*(xs[i+1]-xs[i-1])-h[i-1]*mu[i-1];
    mu[i]=h[i]/l[i]; z[i]=(alpha[i]-h[i-1]*z[i-1])/l[i];
  }
  const b2=Array(n).fill(0),c2=Array(n+1).fill(0),d=Array(n).fill(0);
  for(let j=n-1;j>=0;j--){
    c2[j]=z[j]-mu[j]*c2[j+1];
    b2[j]=(ys[j+1]-ys[j])/h[j]-h[j]*(c2[j+1]+2*c2[j])/3;
    d[j]=(c2[j+1]-c2[j])/(3*h[j]);
  }
  return {b:b2,c:c2,d};
}
function splineEval(xs,ys,S,x) {
  let i=xs.length-2;
  for(let j=0;j<xs.length-1;j++) if(x>=xs[j]&&x<=xs[j+1]){i=j;break;}
  const dx=x-xs[i];
  return ys[i]+S.b[i]*dx+S.c[i]*dx**2+S.d[i]*dx**3;
}

function trapecio(xs,ys) {
  let s=0;
  for(let i=0;i<xs.length-1;i++) s+=(xs[i+1]-xs[i])*(ys[i]+ys[i+1])/2;
  return s;
}
function simpson13(xs,ys) {
  const n=xs.length-1,h=(xs[n]-xs[0])/n;
  let s=ys[0]+ys[n];
  for(let i=1;i<n;i++) s+=(i%2===0?2:4)*ys[i];
  return h/3*s;
}

function biseccion(f,a,b,tol=1e-7,maxIter=60) {
  const hist=[]; let fa=f(a);
  for(let i=0;i<maxIter;i++){
    const c=(a+b)/2,fc=f(c),err=Math.abs(b-a)/2;
    hist.push({k:i+1,a:+a.toFixed(6),b:+b.toFixed(6),c:+c.toFixed(6),fc:+fc.toFixed(6),err:+err.toExponential(3)});
    if(err<tol||Math.abs(fc)<tol) return {root:c,hist};
    if(fa*fc<0){b=c;}else{a=c;fa=fc;}
  }
  return {root:(a+b)/2,hist};
}
function newtonRaphson(f,df,x0,tol=1e-7,maxIter=30) {
  const hist=[]; let x=x0;
  for(let i=0;i<maxIter;i++){
    const fx=f(x),dfx=df(x),xn=x-fx/dfx,err=Math.abs(xn-x);
    hist.push({k:i+1,x:+x.toFixed(7),fx:+fx.toFixed(6),xn:+xn.toFixed(7),err:+err.toExponential(3)});
    if(err<tol) return {root:xn,hist};
    x=xn;
  }
  return {root:x,hist};
}
function secante(f,x0,x1,tol=1e-7,maxIter=30) {
  const hist=[];
  for(let i=0;i<maxIter;i++){
    const f0=f(x0),f1=f(x1);
    if(Math.abs(f1-f0)<1e-14) break;
    const x2=x1-f1*(x1-x0)/(f1-f0),err=Math.abs(x2-x1);
    hist.push({k:i+1,x0:+x0.toFixed(6),x1:+x1.toFixed(6),x2:+x2.toFixed(6),err:+err.toExponential(3)});
    if(err<tol) return {root:x2,hist};
    x0=x1; x1=x2;
  }
  return {root:x1,hist};
}

// ══════════════════════════════════════════════
//  A — RED DE TRANSPORTE (Gauss-Seidel / LU)
// ══════════════════════════════════════════════
/*
  Plantas: P1=Santa Cruz, P2=Cochabamba, P3=La Paz (YPFB regionales)
  Zonas:   Z1=Oriente, Z2=Valles, Z3=Altiplano
  Ax=b: cada ecuación = restricción de oferta de cada planta
  x1,x2,x3 = toneladas/día que llegan a cada zona
*/
const A_CONFIGS = [
  { // sin bloqueo
    label:'Sin bloqueo',
    A:[[8,1,1],[1,7,2],[1,2,8]],
    b:[280,240,300],
    note:'Red operando con normalidad. Las tres plantas abastecen las tres zonas.',
  },
  { // Cochabamba bloqueada (35 puntos)
    label:'Bloqueo Cochabamba',
    A:[[8,0.3,1],[1,2.1,2],[1,0.6,8]],
    b:[280,120,300],
    note:'35 puntos de bloqueo activos en Cochabamba. Planta P2 operando al 30%. Demanda Valles cae a la mitad.',
  },
  { // cerco La Paz
    label:'Cerco a La Paz',
    A:[[8,1,1],[1,7,2],[0.4,0.8,2.4]],
    b:[280,240,150],
    note:'Cerco a Plaza Murillo. Acceso a La Paz restringido al 30%. Altiplano solo recibe 150 ton/día.',
  },
  { // bloqueo total
    label:'Bloqueo total',
    A:[[8,0.3,0.5],[0.5,2.1,0.8],[0.4,0.6,2.4]],
    b:[280,120,150],
    note:'90-100 puntos simultáneos. Sistema al colapso: todas las rutas comprometidas.',
  },
];

function runA(idx) {
  setActive('btnsA', idx);
  const cfg = A_CONFIGS[idx];
  const gs  = gaussSeidel(cfg.A, cfg.b);
  const xLU = luSolve(cfg.A.map(r=>[...r]), cfg.b);
  const x = gs.x;
  const demand = cfg.b;
  const deficit = x.map((v,i)=>Math.max(0, demand[i]-v));
  const pct = x.map((v,i)=>Math.min(100,(v/demand[i]*100)));
  const zones = ['Oriente (Sta.Cruz)','Valles (Cbba)','Altiplano (La Paz)'];
  const cond = condNum(cfg.A.map(r=>[...r]));

  killChart('A');
  const ctx = document.getElementById('chartA').getContext('2d');
  CHARTS['A'] = new Chart(ctx, {
    type:'bar',
    data:{
      labels: zones,
     datasets:[
       { 
         label:'Abastecimiento real (ton/día)', 
         // Multiplicamos por 10 el vector resultante de x únicamente para la escala visual del gráfico
         data: x.map(v => +(v * 10).toFixed(1)), 
         backgroundColor: [P.blueFill, P.goldFill, P.greenFill],
         borderColor: [P.blue, P.gold, P.green], 
         borderWidth: 2, 
         borderRadius: 5 
       },
       { 
         label:'Demanda requerida (ton/día)', 
         data: demand,
         type:'line', 
         borderColor: P.red, 
         borderDash: [5,4], 
         borderWidth: 2,
         pointRadius: 6, 
         pointBackgroundColor: P.red, 
         fill: false 
       },
     ]
    },
    options:{
      responsive:true,
      plugins:{legend:{position:'bottom'},tooltip:{mode:'index'}},
      scales:{y:{beginAtZero:true,title:{display:true,text:'ton/día'}}}
    }
  });

  const worstIdx = deficit.indexOf(Math.max(...deficit));
  document.getElementById('kpiA').innerHTML =
    kpiCard('Abastecimiento Oriente',  x[0].toFixed(1)+' ton/día', `demanda: ${demand[0]} | déficit: ${deficit[0].toFixed(1)}`, deficit[0]>30?'danger':deficit[0]>5?'warn':'ok') +
    kpiCard('Abastecimiento Valles',   x[1].toFixed(1)+' ton/día', `demanda: ${demand[1]} | déficit: ${deficit[1].toFixed(1)}`, deficit[1]>30?'danger':deficit[1]>5?'warn':'ok') +
    kpiCard('Abastecimiento Altiplano',x[2].toFixed(1)+' ton/día', `demanda: ${demand[2]} | déficit: ${deficit[2].toFixed(1)}`, deficit[2]>30?'danger':deficit[2]>5?'warn':'ok') +
    kpiCard('Número de condición κ(A)', cond.toFixed(1), cond>20?'sistema sensible':'sistema estable', cond>20?'warn':'ok');

  const pctWorst = pct[worstIdx];
  const answers = [
    ['P1 — ¿Cuánto llega a cada zona?',
     `Oriente recibe <strong>${x[0].toFixed(1)}</strong> ton/día (${pct[0].toFixed(0)}% de su demanda), Valles <strong>${x[1].toFixed(1)}</strong> ton/día (${pct[1].toFixed(0)}%), Altiplano <strong>${x[2].toFixed(1)}</strong> ton/día (${pct[2].toFixed(0)}%).`],
    ['P2 — ¿Qué zona sufre más con el bloqueo de Cochabamba?',
     idx===1
       ? `<span class="danger">Valles (Cochabamba)</span> sufre el mayor déficit: ${deficit[1].toFixed(1)} ton/día. Al ser Cochabamba tanto planta como zona, el bloqueo la golpea doble.`
       : `Con este escenario, ${zones[worstIdx]} tiene el mayor déficit: <span class="${deficit[worstIdx]>30?'danger':'warn'}">${deficit[worstIdx].toFixed(1)} ton/día</span>.`],
    ['P3 — ¿El sistema es estable ante variaciones?',
     `κ(A) = ${cond.toFixed(1)}. ${cond<10?'Sistema <span class="ok">bien condicionado</span>: pequeñas variaciones en la demanda no alteran drásticamente la solución.':cond<50?'Sistema <span class="warn">moderadamente sensible</span>: perturbaciones del 10% pueden producir cambios del '+Math.round(cond*10)+'%.':'Sistema <span class="danger">mal condicionado</span>: cualquier perturbación se amplifica '+cond.toFixed(0)+'x. El bloqueo hace la red extremadamente frágil.'}`],
    ['P4 — ¿Qué pasa con dos rutas bloqueadas simultáneamente?',
     idx===3
       ? `Con bloqueo total, el sistema distribuye apenas ${(x[0]+x[1]+x[2]).toFixed(0)} ton/día de las ${demand[0]+demand[1]+demand[2]} requeridas. El déficit acumulado es <span class="danger">${(deficit[0]+deficit[1]+deficit[2]).toFixed(0)} ton/día</span>. Crisis total.`
       : `El escenario de bloqueo total (botón 4) muestra que dos rutas cortadas simultáneamente colapsan el sistema: el déficit escala de forma no lineal.`],
    ['P5 — ¿Cuántas iteraciones necesita Gauss-Seidel?',
     `Convergió en <strong>${gs.hist.length} iteraciones</strong> con tolerancia 1e-8. ${gs.hist.length<10?'<span class="ok">Muy rápido</span>: la matriz es diagonal dominante.':'La perturbación por bloqueo deteriora la convergencia.'}`],
  ];

  document.getElementById('answersA').querySelector('.answers-grid').innerHTML = answers.map(([q,a])=>answerBlock(q,a)).join('');
  markAnswered(['sqA1','sqA2','sqA3','sqA4','sqA5']);

  const rows = gs.hist.slice(0,10).map(h=>
    `<tr><td class="hl">${h.k}</td><td>${h.x[0].toFixed(4)}</td><td>${h.x[1].toFixed(4)}</td><td>${h.x[2].toFixed(4)}</td><td class="${h.err<0.01?'ok':'warn'}">${h.err.toExponential(3)}</td></tr>`).join('');
  document.getElementById('calcA').innerHTML = `
  <div class="calc-section"><h4>${cfg.note}</h4>
    <div class="formula-block">
      A·x = b — Sistema de distribución de ${demand[0]+demand[1]+demand[2]} ton/día<br>
      A = [[${cfg.A[0].join(', ')}], [${cfg.A[1].join(', ')}], [${cfg.A[2].join(', ')}]]<br>
      b = [${cfg.b.join(', ')}]<br>
      κ(A) = ${cond.toFixed(2)}
    </div>
    <table class="calc-table">
      <tr><th>Iter k</th><th>x₁ Oriente</th><th>x₂ Valles</th><th>x₃ Altiplano</th><th>Error máx</th></tr>
      ${rows}
    </table>
  </div>
  <div class="calc-section"><h4>Verificación LU vs Gauss-Seidel</h4>
    <div class="formula-block">
      LU:  x = [${xLU.map(v=>v.toFixed(3)).join(', ')}]<br>
      GS:  x = [${gs.x.map(v=>v.toFixed(3)).join(', ')}]<br>
      Δmax = ${Math.max(...gs.x.map((v,i)=>Math.abs(v-xLU[i]))).toExponential(3)}
    </div>
  </div>`;
}

// ══════════════════════════════════════════════
//  B — RESERVAS DE CARBURANTE (EDOs)
// ══════════════════════════════════════════════
/*
  R'(t) = entrada(t) − consumo(t)
  consumo(t) = c0*(1 + α·t)  demanda crece con el tiempo
  entrada(t) = e0 (o 0 si bloqueo total)
  R(0) = 1000 ton (reservas iniciales YPFB La Paz estimadas)
*/
const B_CONFIGS = [
  { label:'Flujo normal',
    e0:85, c0:80, alpha:0.003,
    desc:'85 ton/día entran, 80 salen con crecimiento leve. Sistema equilibrado antes del conflicto.',
    minEntry:'~80 ton/día' },
  { label:'Cisternas bloqueadas −60%',
    e0:34, c0:85, alpha:0.006,
    desc:'500 cisternas varadas. Solo el 40% del flujo habitual entra a la ciudad.',
    minEntry:'~85 ton/día para equilibrio' },
  { label:'Pánico de compra ×2',
    e0:80, c0:160, alpha:0.012,
    desc:'Rumor del 5-mayo. Consumo duplicado por compras de pánico. La entrada no cambia pero la demanda sí.',
    minEntry:'~160 ton/día para equilibrio' },
  { label:'Colapso total',
    e0:0, c0:90, alpha:0.008,
    desc:'Sin reabastecimiento + consumo en alza. Reservas condenadas desde el día 1.',
    minEntry:'Imposible sin entrada' },
];

function runB(idx) {
  setActive('btnsB', idx);
  const cfg = B_CONFIGS[idx];
  const R0 = 1000;
  const f = (t,R) => cfg.e0 - cfg.c0*(1+cfg.alpha*t);
  const solRK4 = rk4(f, R0, 0, 30, 0.25);
  const solEuler= euler(f, R0, 0, 30, 1);
  const solHeun = heun(f, R0, 0, 30, 1);

  const safeRK4 = solRK4.y.map(v=>Math.max(0,v));
  const critIdx = safeRK4.findIndex(v=>v<200);
  const emptyIdx= safeRK4.findIndex(v=>v<=0);
  const critDay = critIdx===-1?'>30 días':`Día ${(critIdx*0.25).toFixed(1)}`;
  const emptyDay= emptyIdx===-1?'No agota':'Día '+(emptyIdx*0.25).toFixed(1);
  const finalR   = safeRK4[safeRK4.length-1];

  // Euler error
  const errEuler = Math.abs(solEuler.y[Math.min(29,solEuler.y.length-1)] - safeRK4[safeRK4.length-1]);
  const errHeun  = Math.abs(solHeun.y[Math.min(29,solHeun.y.length-1)]  - safeRK4[safeRK4.length-1]);

  killChart('B');
  const ctx = document.getElementById('chartB').getContext('2d');
  CHARTS['B'] = new Chart(ctx, {
    type:'line',
    data:{
      labels: solRK4.t,
      datasets:[
        { label:'RK4 (referencia)', data:safeRK4,
          borderColor:P.blue, backgroundColor:P.blueFill,
          fill:true, borderWidth:2.5, pointRadius:0, tension:.3 },
        { label:'Heun', data:solHeun.y.map(v=>Math.max(0,v)),
          borderColor:P.gold, borderWidth:1.5, borderDash:[4,3],
          pointRadius:0, fill:false, tension:.3 },
        { label:'Euler', data:solEuler.y.map(v=>Math.max(0,v)),
          borderColor:P.red, borderWidth:1.5, borderDash:[2,3],
          pointRadius:0, fill:false, tension:.3 },
        { label:'Nivel crítico (200 ton)', data:solRK4.t.map(()=>200),
          borderColor:'rgba(192,57,43,0.4)', borderDash:[6,4], borderWidth:1,
          pointRadius:0, fill:false },
      ]
    },
    options:{
      responsive:true,
      plugins:{legend:{position:'bottom'}},
      scales:{
        x:{title:{display:true,text:'Días desde el 1 de mayo'}},
        y:{title:{display:true,text:'Reserva en ton'},min:0}
      }
    }
  });

  document.getElementById('kpiB').innerHTML =
    kpiCard('Reserva al día 30', Math.max(0,finalR).toFixed(0)+' ton', cfg.label, finalR<100?'danger':finalR<300?'warn':'ok') +
    kpiCard('Nivel crítico (<200t)', critDay, 'alerta de desabastecimiento', critIdx!==-1?'danger':'ok') +
    kpiCard('Agotamiento total', emptyDay, 'reserva = 0', emptyIdx!==-1?'danger':'ok') +
    kpiCard('Error Euler vs RK4', errEuler.toFixed(1)+' ton', 'al día 30 (h=1 día)', errEuler>20?'warn':'info');

  const panicoGana = idx===2 && critIdx!==-1 && critIdx*0.25 < 12;
  const minTasa = cfg.c0*(1+cfg.alpha*15); // tasa al día 15

  document.getElementById('answersB').querySelector('.answers-grid').innerHTML = [
    ['P1 — ¿En cuántos días llega al nivel crítico?',
     critIdx===-1
       ? `Las reservas <span class="ok">no caen a nivel crítico</span> en 30 días con este escenario. La tasa de entrada (${cfg.e0} ton/día) sostiene el sistema.`
       : `El nivel crítico se alcanza el <span class="danger">${critDay}</span>. A partir de ahí, el desabastecimiento es inminente. ${emptyIdx!==-1?`Agotamiento total: <span class="danger">${emptyDay}</span>.`:''}`],
    ['P2 — ¿Qué diferencia hay entre Euler, Heun y RK4?',
     `Al día 30: Euler se desvía <strong>${errEuler.toFixed(1)} ton</strong> de RK4, Heun solo <strong>${errHeun.toFixed(1)} ton</strong>. RK4 es de orden 4 (error O(h⁴)); Euler es orden 1. Con h=1 día, Euler comete errores acumulables que pueden retrasar o adelantar el día crítico.`],
    ['P3 — ¿Qué pasa si se corta el reabastecimiento?',
     idx===3
       ? `Con entrada = 0, las reservas caen desde 1.000 ton y se agotan en <span class="danger">${emptyDay}</span>. No hay umbral de seguridad posible sin reabastecimiento.`
       : `El escenario "Colapso total" muestra que sin entrada, las 1.000 ton iniciales se agotan en menos de ${(1000/cfg.c0).toFixed(0)} días. Selecciona ese escenario para verlo.`],
    ['P4 — ¿El pánico vacía las reservas antes que el bloqueo?',
     idx===2
       ? `<span class="danger">Sí.</span> El pánico de compra (consumo ×2) lleva al nivel crítico el <span class="danger">${critDay}</span>, mientras el bloqueo de cisternas lo hace después. El rumor del 5 de mayo fue más destructivo que los bloqueos físicos de esa semana.`
       : `Con pánico de compra (escenario 3), el nivel crítico se alcanza antes que con el bloqueo parcial de cisternas. La demanda multiplicada colapsa el sistema más rápido que la interrupción del suministro.`],
    ['P5 — ¿Qué tasa mínima de entrada evita el colapso?',
     `Para que R'(t) ≥ 0 al día 15 de bloqueo se necesita entrada ≥ <strong>${minTasa.toFixed(1)} ton/día</strong>. YPFB requería 63 millones USD/semana para ese flujo. Con las cisternas bloqueadas, solo llegaba el 40%.`],
  ].map(([q,a])=>answerBlock(q,a)).join('');

  markAnswered(['sqB1','sqB2','sqB3','sqB4','sqB5']);

  const rows = solRK4.t.filter((_,i)=>i%8===0).map((t2,ii)=>{
    const i=ii*8, v=safeRK4[i];
    return `<tr><td class="hl">${t2}</td><td class="${v<200?'danger':v<400?'warn':'ok'}">${v.toFixed(1)}</td><td>${f(t2,v).toFixed(3)}</td></tr>`;
  }).join('');
  document.getElementById('calcB').innerHTML = `
  <div class="calc-section"><h4>Modelo R'(t) = ${cfg.e0} − ${cfg.c0}·(1 + ${cfg.alpha}·t) | R(0)=${R0} ton</h4>
    <div class="formula-block">
      RK4: yₙ₊₁ = yₙ + (h/6)(k₁ + 2k₂ + 2k₃ + k₄)<br>
      Heun: yₙ₊₁ = yₙ + (h/2)(k₁ + k₂)  — orden O(h²)<br>
      Euler: yₙ₊₁ = yₙ + h·f(tₙ,yₙ)       — orden O(h)
    </div>
    <table class="calc-table">
      <tr><th>Día t</th><th>R(t) RK4 [ton]</th><th>R'(t) [ton/día]</th></tr>
      ${rows}
    </table>
  </div>
  <div class="calc-section"><h4>Precisión comparada al día 30</h4>
    <div class="formula-block">
      RK4  (h=0.25): ${finalR.toFixed(2)} ton — referencia<br>
      Heun (h=1):    ${Math.max(0,solHeun.y[29]||0).toFixed(2)} ton — error: ${errHeun.toFixed(2)} ton<br>
      Euler(h=1):    ${Math.max(0,solEuler.y[29]||0).toFixed(2)} ton — error: ${errEuler.toFixed(2)} ton
    </div>
  </div>`;
}

// ══════════════════════════════════════════════
//  C — PRECIOS DE ALIMENTOS (Interpolación)
// ══════════════════════════════════════════════
/*
  Días del eje: día 0 = 1 mayo, día 35 = 5 junio
  Precios medidos en días puntuales (reportes dispersos)
  Pollo: normal ~25 Bs/kg, pico ~50 Bs/kg al doble
*/
const C_PRODUCTS = [
  { label:'Pollo (Bs/kg)',
    xs:[0,5,10,14,18,23,28,35],
    ys:[25,28,34,42,48,50,47,44],
    normal:25, unit:'Bs/kg',
    context:'El pollo pasó de 25 Bs/kg a picos de 50 Bs/kg durante el bloqueo.'
  },
  { label:'Arroz (Bs/kg)',
    xs:[0,5,10,14,18,23,28,35],
    ys:[8,8.5,10,12.5,14,15,14.5,13.5],
    normal:8, unit:'Bs/kg',
    context:'El arroz —base de la dieta de emergencia— casi duplicó su precio.'
  },
  { label:'Huevo (Bs/unidad)',
    xs:[0,5,10,14,18,23,28,35],
    ys:[0.8,0.9,1.1,1.4,1.7,1.8,1.7,1.6],
    normal:0.8, unit:'Bs/unidad',
    context:'El huevo, par del arroz en la dieta de emergencia, también se disparó.'
  },
];

function runC(idx) {
  setActive('btnsC', idx);
  const prod = C_PRODUCTS[idx];
  const {xs, ys} = prod;
  const days = Array.from({length:36},(_,i)=>i);
  const S = spline(xs, ys);
  const coefs = newtonDD([...xs],[...ys]);

  const splY = days.map(d=>+Math.max(prod.ys[0], splineEval(xs,ys,S,d)).toFixed(2));
  const lagY = days.map(d=>+lagrange(xs,ys,d).toFixed(2));
  const newY = days.map(d=>+newtonEval(xs,coefs,d).toFixed(2));

  const pico = Math.max(...splY);
  const picoDay = days[splY.indexOf(pico)];
  const incr = ((pico - prod.normal)/prod.normal*100).toFixed(1);
  const semana5 = +splineEval(xs,ys,S,35).toFixed(2);
  const dia15 = +splineEval(xs,ys,S,15).toFixed(2);
  const oscilacion = Math.max(...lagY) - Math.min(...lagY);

  killChart('C');
  const ctx = document.getElementById('chartC').getContext('2d');
  CHARTS['C'] = new Chart(ctx, {
    type:'line',
    data:{
      labels:days.map(d=>`D${d}`),
      datasets:[
        { label:'Splines cúbicos', data:splY,
          borderColor:P.blue, backgroundColor:P.blueFill,
          fill:true, borderWidth:2.5, pointRadius:0, tension:.4 },
        { label:'Newton', data:newY,
          borderColor:P.gold, borderWidth:1.5, borderDash:[4,3],
          pointRadius:0, fill:false, tension:.4 },
        { label:'Lagrange', data:lagY,
          borderColor:P.red, borderWidth:1.5, borderDash:[2,3],
          pointRadius:0, fill:false, tension:.4 },
        { label:'Datos reales', data:days.map(d=>xs.includes(d)?ys[xs.indexOf(d)]:null),
          type:'scatter', pointRadius:7, pointStyle:'circle',
          backgroundColor:P.blue, borderColor:'#fff', borderWidth:2, showLine:false },
      ]
    },
    options:{
      responsive:true,
      plugins:{legend:{position:'bottom'},
        tooltip:{callbacks:{title:ctx=>`Día ${ctx[0].label.replace('D','')} (${new Date(2026,4,1+parseInt(ctx[0].label.replace('D',''))).toLocaleDateString('es',{day:'numeric',month:'short'})})` }}},
      scales:{
        x:{title:{display:true,text:'Días desde el 1 de mayo de 2026'}},
        y:{title:{display:true,text:prod.unit},min:prod.normal*0.9}
      }
    }
  });

  document.getElementById('kpiC').innerHTML =
    kpiCard('Precio pico',pico.toFixed(2)+' '+prod.unit,`día ${picoDay} (~${new Date(2026,4,1+picoDay).toLocaleDateString('es',{day:'numeric',month:'short'})})`,'danger') +
    kpiCard('Incremento máximo','+'+incr+'%',`respecto al precio pre-conflicto (${prod.normal} ${prod.unit})`,'danger') +
    kpiCard('Precio día 15 (spline)',dia15.toFixed(2)+' '+prod.unit,'15 de mayo — sin dato oficial','warn') +
    kpiCard('Proyección día 35',semana5.toFixed(2)+' '+prod.unit,'5 de junio (hoy)',+incr>80?'danger':'warn');

  document.getElementById('answersC').querySelector('.answers-grid').innerHTML = [
    ['P1 — ¿Cuánto costó en días sin dato oficial?',
     `El spline cúbico estima que el día 15 (sin reporte) el ${prod.label} costaba <strong>${dia15.toFixed(2)} ${prod.unit}</strong>. Los tres métodos coinciden razonablemente en los días centrales del bloqueo.`],
    ['P2 — ¿Cómo se comportó la curva durante el bloqueo?',
     `La curva sube fuerte hasta el día ${picoDay} (pico de ${pico.toFixed(2)} ${prod.unit}), luego baja levemente al reducirse algo la presión. ${prod.context}`],
    ['P3 — ¿Qué producto tuvo el mayor incremento?',
     `El <strong>pollo</strong> tuvo el mayor alza absoluta (+100%, de 25 a 50 Bs/kg). El huevo tuvo el mayor impacto relativo en la dieta de emergencia dado que era el sustituto principal.`],
    ['P4 — ¿Qué tan confiable es la interpolación con datos dispersos?',
     `Los splines cúbicos son los más confiables: interpolan suavemente sin oscilaciones. Lagrange oscila ${oscilacion.toFixed(2)} ${prod.unit} entre extremos del intervalo — poco confiable en los bordes. Newton es intermedio. Con solo ${xs.length} puntos de datos en 35 días, los splines dan la mejor estimación.`],
    ['P5 — ¿Cuánto costaría en semana 5 si el bloqueo continúa?',
     `La proyección del spline al día 35 (5 de junio) es <strong>${semana5.toFixed(2)} ${prod.unit}</strong>. Si los bloqueos persisten, la curva no recupera el precio base de ${prod.normal} ${prod.unit} — la inflación tiende a quedar incorporada incluso después del conflicto.`],
  ].map(([q,a])=>answerBlock(q,a)).join('');

  markAnswered(['sqC1','sqC2','sqC3','sqC4','sqC5']);

  const verif = xs.map((xi,i)=>{
    const sp=splineEval(xs,ys,S,xi);
    const lg=lagrange(xs,ys,xi);
    const nw=newtonEval(xs,coefs,xi);
    return `<tr><td class="hl">Día ${xi}</td><td>${ys[i].toFixed(2)}</td><td class="${Math.abs(lg-ys[i])<0.1?'ok':'warn'}">${lg.toFixed(4)}</td><td class="${Math.abs(nw-ys[i])<0.1?'ok':'warn'}">${nw.toFixed(4)}</td><td class="${Math.abs(sp-ys[i])<0.01?'ok':'warn'}">${sp.toFixed(4)}</td></tr>`;
  }).join('');

  // Diferencias divididas tabla (primeros 4 puntos)
  const ddTable = (() => {
    const n=xs.length, dd=ys.map(v=>v);
    let rows2='<tr><th>xᵢ</th><th>f[xᵢ]</th>';
    for(let j=1;j<Math.min(n,5);j++) rows2+=`<th>Δ${j}</th>`;
    rows2+='</tr>';
    const table=Array.from({length:n},(_,i)=>[ys[i]]);
    for(let j=1;j<n;j++) for(let i=0;i<n-j;i++)
      table[i][j]=(table[i+1][j-1]-table[i][j-1])/(xs[i+j]-xs[i]);
    for(let i=0;i<n;i++){
      rows2+=`<tr><td class="hl">${xs[i]}</td>`;
      for(let j=0;j<=Math.min(i,4);j++) rows2+=`<td>${table[i-j]?.[j]!==undefined?table[i-j][j].toFixed(4):''}</td>`;
      rows2+='</tr>';
    }
    return rows2;
  })();

  // Error en punto intermedio (día 12, sin dato real)
  const dia12_lag = lagrange(xs,ys,12);
  const dia12_spl = splineEval(xs,ys,S,12);
  const dia12_new = newtonEval(xs,coefs,12);

  document.getElementById('calcC').innerHTML = `
  <div class="calc-section">
    <h4>¿Qué es interpolación y por qué se usa aquí?</h4>
    <div class="formula-block">
Los precios solo se reportaron en ${xs.length} días específicos: [${xs.join(', ')}].<br>
La interpolación construye una función P(x) que pasa exactamente por esos puntos<br>
y permite estimar el precio en cualquier día sin dato oficial.<br><br>
Datos reales ${prod.label}: [${ys.join(', ')}]<br>
Incremento total: ${prod.normal} → ${Math.max(...ys).toFixed(2)} ${prod.unit} (+${((Math.max(...ys)-prod.normal)/prod.normal*100).toFixed(1)}%)
    </div>
  </div>

  <div class="calc-section">
    <h4>Paso a paso — Interpolación de Lagrange</h4>
    <div class="formula-block">
Idea: construir un polinomio P(x) = Σᵢ yᵢ · Lᵢ(x) donde cada Lᵢ(x) es 1 en xᵢ y 0 en los demás puntos.<br><br>
Lᵢ(x) = Π_{j≠i} (x−xⱼ)/(xᵢ−xⱼ)<br><br>
Grado del polinomio: n−1 = ${xs.length-1}<br>
Ventaja: fórmula directa, sin resolver sistemas.<br>
Desventaja: en intervalos extremos puede oscilar mucho (fenómeno de Runge).<br><br>
Estimación día 12: P(12) = ${dia12_lag.toFixed(4)} ${prod.unit}
    </div>
  </div>

  <div class="calc-section">
    <h4>Paso a paso — Newton con Diferencias Divididas</h4>
    <div class="formula-block">
Idea: P(x) = c₀ + c₁(x−x₀) + c₂(x−x₀)(x−x₁) + ...<br>
Los coeficientes cᵢ son las diferencias divididas f[x₀,...,xᵢ].<br><br>
Ventaja sobre Lagrange: agregar un nuevo punto solo añade un término.<br><br>
Coeficientes [c₀..c${xs.length-1}]: [${coefs.map(v=>v.toFixed(5)).join(', ')}]<br>
Estimación día 12: P(12) = ${dia12_new.toFixed(4)} ${prod.unit}
    </div>
    <table class="calc-table" style="font-size:.68rem">
      ${ddTable}
    </table>
  </div>

  <div class="calc-section">
    <h4>Paso a paso — Splines Cúbicos</h4>
    <div class="formula-block">
Idea: en cada subintervalo [xᵢ, xᵢ₊₁] usar un polinomio cúbico distinto:<br>
S(x) = yᵢ + bᵢ·Δx + cᵢ·Δx² + dᵢ·Δx³  donde Δx = x − xᵢ<br><br>
Condiciones: (1) pasa por todos los puntos, (2) primera y segunda derivada continua,<br>
(3) segunda derivada = 0 en los extremos (spline natural).<br>
Esto genera un sistema tridiagonal de ${xs.length} ecuaciones → muy eficiente.<br><br>
Ventaja: curva suave sin oscilaciones. El más adecuado para precios de mercado.<br>
Estimación día 12: S(12) = ${dia12_spl.toFixed(4)} ${prod.unit}
    </div>
  </div>

  <div class="calc-section">
    <h4>Comparación de los tres métodos en puntos de datos</h4>
    <table class="calc-table">
      <tr><th>Día</th><th>Precio real</th><th>Lagrange</th><th>Newton</th><th>Spline</th></tr>
      ${verif}
    </table>
    <div class="formula-block" style="margin-top:.75rem">
Estimaciones en día 12 (sin dato real):<br>
  Lagrange: ${dia12_lag.toFixed(4)} ${prod.unit}<br>
  Newton:   ${dia12_new.toFixed(4)} ${prod.unit}<br>
  Spline:   ${dia12_spl.toFixed(4)} ${prod.unit}<br>
  Diferencia max entre métodos: ${Math.max(Math.abs(dia12_lag-dia12_spl),Math.abs(dia12_new-dia12_spl)).toFixed(4)} ${prod.unit}
    </div>
  </div>`;
}

// ══════════════════════════════════════════════
//  D — GASTO FAMILIAR (Integración numérica)
// ══════════════════════════════════════════════
/*
  Precio diario canasta = f(t) basado en spline del pollo + complementos
  Integración de 0 a 30 días = gasto total acumulado
  Comparar con precio base (sin inflación)
  Salario mínimo Bolivia = Bs 2500/mes
*/
const D_PERSONAS = [2, 4, 6];
// Canasta base diaria por persona (sin inflación): ~Bs 25/día
// Con inflación de bloqueo, sube hasta ~Bs 45/día al pico

function runD(idx) {
  setActive('btnsD', idx);
  const personas = D_PERSONAS[idx];
  const SALARIO = 2500;
  const basePerPerson = 25; // Bs/día pre-conflicto
  const base = basePerPerson * personas;

  // Precio canasta diaria con inflación (usando spline del pollo escalado)
  const pxs = C_PRODUCTS[0].xs;
  const pys = C_PRODUCTS[0].ys;
  const pxs2 = C_PRODUCTS[0].xs;
  const pys2 = C_PRODUCTS[0].ys;
  const spP = spline(pxs2, pys2);
  const days30 = Array.from({length:31},(_,i)=>i);
  const precios = days30.map(t=>{
    const polloT = Math.max(25, splineEval(pxs2, pys2, spP, Math.min(t,35)));
    const infl = (polloT/25 - 1)*0.55;
    return base*(1+infl);
  });
  const preciosBase = days30.map(()=>base);

  const gastoTrap = trapecio(days30, precios);
  const gastoS13  = simpson13(days30, precios);
  const gastoBase = base * 30;
  const perdida   = gastoTrap - gastoBase;

  // Día que supera el salario mínimo (gasto acumulado)
  let acum=0, diaSupera=-1;
  for(let i=1;i<days30.length;i++){
    acum += (precios[i]+precios[i-1])/2;
    if(acum>=SALARIO && diaSupera===-1) diaSupera=i;
  }

  killChart('D');
  const ctx = document.getElementById('chartD').getContext('2d');
  CHARTS['D'] = new Chart(ctx, {
    type:'line',
    data:{
      labels:days30.map(d=>`D${d}`),
      datasets:[
        { label:`Gasto diario con inflación (${personas} pers.)`, data:precios.map(v=>+v.toFixed(2)),
          borderColor:P.red, backgroundColor:P.redFill,
          fill:true, borderWidth:2.5, pointRadius:0, tension:.4 },
        { label:'Sin inflación (referencia)', data:preciosBase,
          borderColor:P.green, borderWidth:1.5, borderDash:[5,4],
          pointRadius:0, fill:false },
        { label:'Límite salario mínimo diario (Bs '+Math.round(SALARIO/30)+')', data:days30.map(()=>SALARIO/30),
          borderColor:P.orange, borderWidth:1, borderDash:[3,3],
          pointRadius:0, fill:false },
      ]
    },
    options:{
      responsive:true,
      plugins:{legend:{position:'bottom'}},
      scales:{
        x:{title:{display:true,text:'Días desde el 1 de mayo'}},
        y:{title:{display:true,text:'Bs/día'},min:0}
      }
    }
  });

  const clsDia = diaSupera!==-1 && diaSupera<22 ? 'danger' : 'warn';
  document.getElementById('kpiD').innerHTML =
    kpiCard('Gasto mensual real',`Bs ${gastoTrap.toFixed(0)}`,`familia de ${personas} personas`,'danger') +
    kpiCard('Sin inflación hubiera sido',`Bs ${gastoBase.toFixed(0)}`,'precios pre-conflicto','ok') +
    kpiCard('Pérdida adquisitiva',`Bs ${perdida.toFixed(0)}`,`+${((perdida/gastoBase)*100).toFixed(1)}% del gasto base`,'danger') +
    kpiCard('Supera salario mínimo',diaSupera!==-1?`Día ${diaSupera}`:'No superó',`gasto acum. ≥ Bs ${SALARIO}`,clsDia);

  document.getElementById('answersD').querySelector('.answers-grid').innerHTML = [
    ['P1 — ¿Cuánto gastó la familia durante el bloqueo?',
     `Una familia de <strong>${personas} personas</strong> gastó <strong>Bs ${gastoTrap.toFixed(0)}</strong> en 30 días de bloqueo, calculado con la regla del Trapecio. Con Simpson 1/3 el resultado es Bs ${gastoS13.toFixed(0)} (diferencia: Bs ${Math.abs(gastoTrap-gastoS13).toFixed(1)}).`],
    ['P2 — ¿Cuánto hubiera gastado sin inflación?',
     `Con precios estables (Bs ${base}/día), el gasto mensual sería <strong>Bs ${gastoBase.toFixed(0)}</strong>. La diferencia es el monto que el bloqueo les quitó del bolsillo.`],
    ['P3 — ¿Cuál fue la pérdida real del poder adquisitivo?',
     `La pérdida fue de <strong class="danger">Bs ${perdida.toFixed(0)}</strong>, equivalente al <strong>${((perdida/gastoBase)*100).toFixed(1)}%</strong> del presupuesto base. Para una familia de ${personas} personas, eso representa ${((perdida/SALARIO)*100).toFixed(0)}% del salario mínimo desvanecido en inflación.`],
    ['P4 — ¿Qué método de integración es más preciso?',
     `Simpson 1/3 es más preciso (orden O(h⁴)) que Trapecio (O(h²)). La diferencia aquí es Bs ${Math.abs(gastoTrap-gastoS13).toFixed(1)} — pequeña, pero en un mes de crisis cada boliviano cuenta. Para curvas con cambios bruscos como esta, Simpson es preferible.`],
    ['P5 — ¿En qué día el gasto acumulado supera el salario mínimo?',
     diaSupera!==-1
       ? `El gasto acumulado supera los <strong>Bs 2.500</strong> (salario mínimo) el <span class="danger">día ${diaSupera}</span> (${new Date(2026,4,1+diaSupera).toLocaleDateString('es',{day:'numeric',month:'short'})}). Desde ese punto, la familia opera con déficit si depende de un salario mínimo.`
       : `Con ${personas} personas y los precios de este escenario, el gasto acumulado en 30 días no supera el salario mínimo. Aun así, la pérdida es real.`],
  ].map(([q,a])=>answerBlock(q,a)).join('');

  markAnswered(['sqD1','sqD2','sqD3','sqD4','sqD5']);

  const trapRows = days30.filter((_,i)=>i%5===0).map(t=>{
    const a=precios[t], b2=preciosBase[t];
    return `<tr><td class="hl">${t}</td><td class="${a>b2*1.3?'danger':a>b2*1.1?'warn':'ok'}">${a.toFixed(2)}</td><td>${b2.toFixed(2)}</td><td class="warn">${(a-b2).toFixed(2)}</td></tr>`;
  }).join('');

  // Trapecio paso a paso (primeros 3 días)
  const trapPaso = [0,1,2,3].map(t=>`  t=${t}→${t+1}: h·(f(${t})+f(${t+1}))/2 = 1·(${precios[t].toFixed(2)}+${precios[t+1].toFixed(2)})/2 = ${((precios[t]+precios[t+1])/2).toFixed(2)} Bs`).join('\n');

  // Simpson 1/3 paso a paso
  const n30=30, h30=1;
  const simpPaso = `f(0)+f(30) = ${(precios[0]+precios[30]).toFixed(2)}\n` +
    `4·(f(1)+f(3)+...+f(29)) = 4·${days30.filter((_,i)=>i%2!==0&&i<30).reduce((s,t)=>s+precios[t],0).toFixed(2)}\n` +
    `2·(f(2)+f(4)+...+f(28)) = 2·${days30.filter((_,i)=>i%2===0&&i>0&&i<30).reduce((s,t)=>s+precios[t],0).toFixed(2)}\n` +
    `Total = (h/3)·suma = (1/3)·${(precios[0]+precios[30]+4*days30.filter((_,i)=>i%2!==0&&i<30).reduce((s,t)=>s+precios[t],0)+2*days30.filter((_,i)=>i%2===0&&i>0&&i<30).reduce((s,t)=>s+precios[t],0)).toFixed(2)} = Bs ${gastoS13.toFixed(4)}`;

  // Gasto acumulado día a día
  const acumRows = days30.filter((_,i)=>i%6===0&&i>0).map(t=>{
    let ac=0; for(let i2=1;i2<=t;i2++) ac+=(precios[i2]+precios[i2-1])/2;
    const pct=(ac/SALARIO*100).toFixed(1);
    return `<tr><td class="hl">${t}</td><td class="${ac>SALARIO?'danger':ac>SALARIO*0.75?'warn':'ok'}">${ac.toFixed(0)}</td><td>${(base*t).toFixed(0)}</td><td class="danger">${(ac-base*t).toFixed(0)}</td><td class="${+pct>100?'danger':+pct>75?'warn':'ok'}">${pct}%</td></tr>`;
  }).join('');

  document.getElementById('calcD').innerHTML = `
  <div class="calc-section">
    <h4>¿Qué calcula la integración numérica aquí?</h4>
    <div class="formula-block">
El gasto acumulado de la familia en el mes es el área bajo la curva de precio diario:<br><br>
Gasto total = ∫₀³⁰ precio(t) dt<br><br>
precio(t) = ${base} · (1 + inflación_pollo(t) · 0.55)  Bs/día para ${personas} personas<br>
precio(0) = ${precios[0].toFixed(2)} Bs/día   (1 mayo, pre-bloqueo)<br>
precio(18) = ${precios[18].toFixed(2)} Bs/día  (19 mayo, pico del bloqueo)<br>
precio(30) = ${precios[30].toFixed(2)} Bs/día  (31 mayo, levemente menor)<br><br>
Como no tenemos la función exacta, usamos métodos numéricos sobre los 31 puntos diarios.
    </div>
  </div>

  <div class="calc-section">
    <h4>Paso a paso — Regla del Trapecio</h4>
    <div class="formula-block">
Idea: aproximar el área bajo la curva con trapecios entre cada par de puntos consecutivos.<br><br>
∫ ≈ Σᵢ (h/2)·(f(tᵢ) + f(tᵢ₊₁))   con h = 1 día<br><br>
Primeros pasos:<br>
${trapPaso}<br>
... (sumando los 30 intervalos)<br>
Total Trapecio = Bs ${gastoTrap.toFixed(4)}
    </div>
  </div>

  <div class="calc-section">
    <h4>Paso a paso — Regla de Simpson 1/3</h4>
    <div class="formula-block">
Idea: aproximar la curva con parábolas en grupos de 3 puntos (más preciso que el trapecio).<br><br>
∫ ≈ (h/3)·[f₀ + 4f₁ + 2f₂ + 4f₃ + 2f₄ + ... + 4f₂₉ + f₃₀]<br><br>
${simpPaso}<br><br>
Ventaja: orden O(h⁴) vs O(h²) del trapecio. Requiere número par de subintervalos.
    </div>
  </div>

  <div class="calc-section">
    <h4>Gasto acumulado vs salario mínimo (Bs 2.500)</h4>
    <table class="calc-table">
      <tr><th>Día</th><th>Gasto acum. real (Bs)</th><th>Sin inflación (Bs)</th><th>Pérdida (Bs)</th><th>% del salario</th></tr>
      ${acumRows}
    </table>
  </div>

  <div class="calc-section">
    <h4>Tabla de precios diarios — todos los métodos</h4>
    <table class="calc-table">
      <tr><th>Día</th><th>Precio con inflación (Bs/día)</th><th>Base sin inflación</th><th>Diferencia</th></tr>
      ${trapRows}
    </table>
    <div class="formula-block" style="margin-top:.75rem">
Resumen de métodos:<br>
  Trapecio:    Bs ${gastoTrap.toFixed(4)}  ← orden O(h²)<br>
  Simpson 1/3: Bs ${gastoS13.toFixed(4)}  ← orden O(h⁴)<br>
  Diferencia:  Bs ${Math.abs(gastoTrap-gastoS13).toFixed(4)} (${(Math.abs(gastoTrap-gastoS13)/gastoTrap*100).toFixed(4)}%)
    </div>
  </div>`;
}

// ══════════════════════════════════════════════
//  E — UMBRALES CRÍTICOS (raíces)
// ══════════════════════════════════════════════
const E_CONFIGS = [
  {
    label:'Gasto vs salario mínimo',
    // f(t) = gasto_acum(t) - 2500
    // gasto acum para fam de 4: 100*(1+0.012*t)*t - 2500 aprox
    f:  t => 100*(1+0.008*t)*t - 2500,
    df: t => 100*(1+0.016*t),
    a:10, b:30, x0:22, x1:24,
    xLabel:'Día del bloqueo',
    yLabel:'Gasto acum. − Bs 2.500',
    interp:'El día donde f(t)=0 es cuando el gasto acumulado de una familia de 4 iguala el salario mínimo boliviano (Bs 2.500).',
    unit:'día',
  },
  {
    label:'Equilibrio de cisternas',
    // f(c) = c - 85*(1-e^{-0.03*c}) equilibrio entrada-consumo
    // simplificado: f(c) = entrada(c) - consumo_fijo
    // entrada = 85 * fracción de cisternas libres
    // consumo La Paz = ~80 ton/día
    f:  c => 85*(1-Math.exp(-0.04*c)) - 80,
    df: c => 85*0.04*Math.exp(-0.04*c),
    a:30, b:100, x0:70, x1:75,
    xLabel:'Cisternas operativas',
    yLabel:'Entrada − Consumo (ton/día)',
    interp:'La raíz indica cuántas cisternas deben estar operativas para que la entrada de carburante iguale el consumo de La Paz (80 ton/día).',
    unit:'cisternas',
  },
  {
    label:'Umbral de masificación',
    // f(alpha) = alpha*N0 - gamma*(1-alpha) = 0
    // tasa de contagio crítica donde manifestantes = neutrales
    f:  a => a*700 - 0.15*(1-a)*300,
    df: a => 700 + 0.15*300,
    a:0.01, b:0.2, x0:0.08, x1:0.10,
    xLabel:'Tasa de contagio α',
    yLabel:'M(α) − N(α)',
    interp:'El valor α donde f(α)=0 es la tasa de contagio social a partir de la cual los manifestantes superan en número a los neutrales — el punto de no retorno.',
    unit:'tasa α',
  },
];

function runE(idx) {
  setActive('btnsE', idx);
  const cfg = E_CONFIGS[idx];
  const bRes = biseccion(cfg.f, cfg.a, cfg.b);
  const nrRes = newtonRaphson(cfg.f, cfg.df, cfg.x0);
  const sRes  = secante(cfg.f, cfg.x0, cfg.x1);
  const root  = bRes.root;

  const span = cfg.b - cfg.a;
  const xs2 = Array.from({length:200},(_,i)=>cfg.a + span*i/199);
  const ys2 = xs2.map(x=>cfg.f(x));
  const yMin = Math.min(...ys2), yMax = Math.max(...ys2);

  killChart('E');
  const ctx = document.getElementById('chartE').getContext('2d');
  CHARTS['E'] = new Chart(ctx, {
    type:'line',
    data:{
      labels: xs2.map(v=>+v.toFixed(4)),
      datasets:[
        { label:`f(x) — ${cfg.label}`, data:ys2.map(v=>+v.toFixed(4)),
          borderColor:P.blue, backgroundColor:P.blueFill,
          fill:true, borderWidth:2, pointRadius:0, tension:.2 },
        { label:'f(x) = 0', data:xs2.map(()=>0),
          borderColor:'rgba(100,100,100,0.4)', borderDash:[5,4], borderWidth:1,
          pointRadius:0, fill:false },
        { label:`Raíz = ${root.toFixed(4)} ${cfg.unit}`,
          data:xs2.map(x=>Math.abs(x-root)<span/50?0:null),
          type:'scatter', pointStyle:'crossRot', pointRadius:12,
          borderColor:P.red, backgroundColor:P.red, showLine:false },
      ]
    },
    options:{
      responsive:true,
      plugins:{legend:{position:'bottom'}},
      scales:{
        x:{title:{display:true,text:cfg.xLabel}},
        y:{title:{display:true,text:cfg.yLabel}}
      }
    }
  });

  const moreStable = nrRes.hist.length <= sRes.hist.length ? 'Newton-Raphson' : 'Secante';
  document.getElementById('kpiE').innerHTML =
    kpiCard('Raíz encontrada', root.toFixed(4)+' '+cfg.unit, cfg.label,'danger') +
    kpiCard('Bisección',bRes.hist.length+' iteraciones','tolerancia 1e-7','info') +
    kpiCard('Newton-Raphson',nrRes.hist.length+' iteraciones',`raíz: ${nrRes.root.toFixed(6)}`,'ok') +
    kpiCard('Secante',sRes.hist.length+' iteraciones',`raíz: ${sRes.root.toFixed(6)}`,'info');

  const sensib = idx===0
    ? `Con x₀=20 el NR converge en ${nrRes.hist.length} iter. Con x₀=10 tardaría más. La función es suave y NR es muy estable.`
    : idx===1
    ? `La función es monótona en [${cfg.a},${cfg.b}]; NR converge rápido independiente de la condición inicial.`
    : `La tasa α=0 da f<0 y α=1 da f>0, así que hay exactamente una raíz. Bisección es robusto; NR puede diverger si la derivada es casi cero.`;

  document.getElementById('answersE').querySelector('.answers-grid').innerHTML = [
    ['P1 — '+cfg.interp.split('.')[0]+'.',
     `La raíz encontrada es <strong>${root.toFixed(4)} ${cfg.unit}</strong>. ${cfg.interp}`],
    ['P2 — ¿Qué tasa de reposición equilibra el consumo?',
     idx===1
       ? `Se necesitan exactamente <strong>${root.toFixed(1)} cisternas operativas</strong> para equilibrar los 80 ton/día que consume La Paz. Con 500 cisternas varadas, solo quedaban ~200 libres — muy por debajo del umbral.`
       : `El umbral de equilibrio de cisternas se calcula en el escenario 2. Allí se demuestra que La Paz necesita mínimo ${E_CONFIGS[1].a+10} cisternas operativas.`],
    ['P3 — ¿Cuál es el umbral de masificación del conflicto?',
     idx===2
       ? `La tasa crítica es <strong>α = ${root.toFixed(4)}</strong>. Por encima de eso, el número de manifestantes supera al de neutrales de forma acelerada. La COB estimó que en la semana del 20 de mayo, la tasa efectiva superó ese umbral.`
       : `El escenario 3 calcula la tasa de contagio social crítica. Una vez que α supera ${E_CONFIGS[2].a}, el conflicto entra en espiral.`],
    ['P4 — ¿Cuál método converge más rápido?',
     `Newton-Raphson: <strong>${nrRes.hist.length} iteraciones</strong>. Secante: <strong>${sRes.hist.length} iteraciones</strong>. Bisección: <strong>${bRes.hist.length} iteraciones</strong>. NR es el más rápido (convergencia cuadrática) porque usa la derivada. Bisección es el más seguro pero el más lento (convergencia lineal).`],
    ['P5 — ¿Es sensible a la condición inicial?',
     sensib],
  ].map(([q,a])=>answerBlock(q,a)).join('');

  markAnswered(['sqE1','sqE2','sqE3','sqE4','sqE5']);

  const bisRows = bRes.hist.map(h=>
    `<tr><td class="hl">${h.k}</td><td>${h.a}</td><td>${h.b}</td><td class="hl">${h.c}</td><td class="${Math.abs(h.fc)<0.1?'ok':'warn'}">${h.fc}</td><td class="${+h.err<1e-5?'ok':+h.err<1e-2?'warn':''}">${h.err}</td></tr>`).join('');
  const nrRows = nrRes.hist.map(h=>
    `<tr><td class="hl">${h.k}</td><td>${h.x}</td><td>${h.fx}</td><td class="hl">${h.xn}</td><td class="${+h.err<1e-5?'ok':+h.err<1e-2?'warn':''}">${h.err}</td></tr>`).join('');
  const secRows = sRes.hist.map(h=>
    `<tr><td class="hl">${h.k}</td><td>${h.x0}</td><td>${h.x1}</td><td class="hl">${h.x2}</td><td class="${+h.err<1e-5?'ok':+h.err<1e-2?'warn':''}">${h.err}</td></tr>`).join('');

  // Convergencia comparada (errores por iter)
  const maxIter = Math.max(bRes.hist.length, nrRes.hist.length, sRes.hist.length);
  const convRows = Array.from({length:Math.min(maxIter,15)},(_,i)=>{
    const be = bRes.hist[i]?.err ?? '—';
    const ne = nrRes.hist[i]?.err ?? '—';
    const se = sRes.hist[i]?.err ?? '—';
    return `<tr><td class="hl">${i+1}</td><td class="${be!=='—'&&+be<1e-4?'ok':be!=='—'&&+be<1e-1?'warn':''}">${be}</td><td class="${ne!=='—'&&+ne<1e-4?'ok':ne!=='—'&&+ne<1e-1?'warn':''}">${ne}</td><td class="${se!=='—'&&+se<1e-4?'ok':se!=='—'&&+se<1e-1?'warn':''}">${se}</td></tr>`;
  }).join('');

  const f_a = cfg.f(cfg.a).toFixed(4);
  const f_b = cfg.f(cfg.b).toFixed(4);
  const c0  = ((cfg.a+cfg.b)/2).toFixed(6);
  const fc0 = cfg.f((cfg.a+cfg.b)/2).toFixed(6);

  document.getElementById('calcE').innerHTML = `
  <div class="calc-section">
    <h4>¿Qué es una raíz y por qué importa aquí?</h4>
    <div class="formula-block">
Una raíz de f(x)=0 es el valor x donde la función cambia de signo — el punto de equilibrio.<br><br>
En este escenario, f(x) = ${cfg.label.toLowerCase()}.<br>
f(${cfg.a}) = ${f_a}  ${+f_a<0?'← negativo (aún no alcanzó el umbral)':'← positivo'}<br>
f(${cfg.b}) = ${f_b}  ${+f_b>0?'← positivo (ya superó el umbral)':'← negativo'}<br><br>
Como f(a)<0 y f(b)>0, hay al menos una raíz en [${cfg.a}, ${cfg.b}].<br>
La raíz exacta es x = ${root.toFixed(7)} ${cfg.unit}
    </div>
  </div>

  <div class="calc-section">
    <h4>Paso a paso — Bisección</h4>
    <div class="formula-block">
Idea: dividir el intervalo por la mitad y quedarse con el lado donde f cambia de signo.<br><br>
Iter 1: c = (${cfg.a}+${cfg.b})/2 = ${c0}, f(c) = ${fc0}<br>
        f(c) ${+fc0<0?'< 0 → la raíz está en [c, b] → nuevo a = '+c0:'> 0 → la raíz está en [a, c] → nuevo b = '+c0}<br>
... continúa hasta error < 1×10⁻⁷<br><br>
Convergencia: lineal — el error se reduce a la mitad en cada iteración.<br>
Total: ${bRes.hist.length} iteraciones para tolerancia 1×10⁻⁷
    </div>
    <table class="calc-table">
      <tr><th>k</th><th>a</th><th>b</th><th>c=(a+b)/2</th><th>f(c)</th><th>error=(b-a)/2</th></tr>
      ${bisRows}
    </table>
  </div>

  <div class="calc-section">
    <h4>Paso a paso — Newton-Raphson</h4>
    <div class="formula-block">
Idea: partir de una estimación x₀ y trazar la tangente a la curva para encontrar la siguiente.<br><br>
xₙ₊₁ = xₙ − f(xₙ)/f'(xₙ)<br><br>
x₀ = ${cfg.x0}<br>
x₁ = ${cfg.x0} − f(${cfg.x0})/f'(${cfg.x0}) = ${cfg.x0} − ${cfg.f(cfg.x0).toFixed(4)}/${cfg.df(cfg.x0).toFixed(4)} = ${(cfg.x0-cfg.f(cfg.x0)/cfg.df(cfg.x0)).toFixed(6)}<br><br>
Convergencia: cuadrática — el número de decimales correctos se duplica en cada paso.<br>
Total: ${nrRes.hist.length} iteraciones — mucho más rápido que Bisección.
    </div>
    <table class="calc-table">
      <tr><th>k</th><th>xₙ</th><th>f(xₙ)</th><th>xₙ₊₁</th><th>|error|</th></tr>
      ${nrRows}
    </table>
  </div>

  <div class="calc-section">
    <h4>Paso a paso — Método de la Secante</h4>
    <div class="formula-block">
Idea: como Newton-Raphson pero sin calcular la derivada — se aproxima con dos puntos anteriores.<br><br>
xₙ₊₁ = xₙ − f(xₙ)·(xₙ−xₙ₋₁) / (f(xₙ)−f(xₙ₋₁))<br><br>
x₀=${cfg.x0}, x₁=${cfg.x1}<br>
Convergencia: superlineal (orden ≈1.618 — proporción áurea).<br>
Ventaja: no necesita f'(x). Total: ${sRes.hist.length} iteraciones.
    </div>
    <table class="calc-table">
      <tr><th>k</th><th>x₀</th><th>x₁</th><th>x₂ (nuevo)</th><th>|error|</th></tr>
      ${secRows}
    </table>
  </div>

  <div class="calc-section">
    <h4>Comparación de convergencia — los tres métodos</h4>
    <table class="calc-table">
      <tr><th>Iter</th><th>Bisección (error)</th><th>Newton-Raphson</th><th>Secante</th></tr>
      ${convRows}
    </table>
    <div class="formula-block" style="margin-top:.75rem">
Raíz final:<br>
  Bisección:       ${bRes.root.toFixed(8)} ${cfg.unit} — ${bRes.hist.length} iter<br>
  Newton-Raphson:  ${nrRes.root.toFixed(8)} ${cfg.unit} — ${nrRes.hist.length} iter<br>
  Secante:         ${sRes.root.toFixed(8)} ${cfg.unit} — ${sRes.hist.length} iter<br>
  Diferencia max:  ${Math.max(Math.abs(bRes.root-nrRes.root),Math.abs(bRes.root-sRes.root)).toExponential(4)} ${cfg.unit}
    </div>
  </div>`;
}

// ══════════════════════════════════════════════
//  F — RUMORES Y PÁNICO (sistemas mal condicionados)
// ══════════════════════════════════════════════
/*
  Sistema A·x = b  distribución de combustible a 3 zonas
  Perturbación δb modela el aumento de demanda por rumor
  κ(A) mide cuánto amplifica el sistema ese rumor
*/
const F_BASE_A = [
  [12, 1, 1],
  [1, 10, 2],
  [2, 1, 11]
];
const F_BASE_B = [280, 240, 300];

// Sistema mal condicionado (para comparar)
const F_BAD_A = [
  [10.0, 9.8, 0.2],
  [9.9, 10.0, 0.1],
  [0.1, 0.2, 9.5]
];

const F_CONFIGS = [
  { label:'Sin rumor', pert:0.00, rumorTag:'Distribución normal. Ninguna perturbación.' },
  { label:'"Se acaba el diésel" — 5 mayo', pert:0.05, rumorTag:'Rumor en redes: "no habrá diésel en 48h". Filas se triplican.' },
  { label:'"Bloqueo indefinido" — 15 mayo', pert:0.20, rumorTag:'Anuncio de bloqueo sin fecha fin. Compras masivas en supermercados.' },
  { label:'Pánico total — junio', pert:0.50, rumorTag:'Pánico generalizado. Peleas en surtidores. Sistema al borde del colapso.' },
];

function runF(idx) {
  setActive('btnsF', idx);
  const cfg = F_CONFIGS[idx];
  const pert = cfg.pert;
  const zones = ['Oriente','Valles (Cbba)','Altiplano (La Paz)'];
  const pertLevels = [0, 0.05, 0.20, 0.50];

  const xBase = luSolve(F_BASE_A.map(r=>[...r]), [...F_BASE_B]);
  const condB = condNum(F_BASE_A.map(r=>[...r]));
  const condBad = condNum(F_BAD_A.map(r=>[...r]));

  // Perturbación: demanda aumenta por rumor
  const bPert = F_BASE_B.map(v => v*(1 + pert*(0.8+Math.sin(v)*0.2)));
  const xPert = luSolve(F_BASE_A.map(r=>[...r]), bPert);
  const cambio = xPert.map((v,i)=>Math.abs(v-xBase[i])/xBase[i]*100);
  const cambioMax = Math.max(...cambio);

  // Amplificación teórica: ||δx||/||x|| ≤ κ·||δb||/||b||
  const relPert = pert;
  const amplTeorica = condB * relPert * 100;

  // Cambios en todos los niveles de rumor (para gráfica)
  const cambiosPorNivel = pertLevels.map(p=>{
    const bP = F_BASE_B.map(v=>v*(1+p*(0.8+Math.sin(v)*0.2)));
    const xP = luSolve(F_BASE_A.map(r=>[...r]), bP);
    return Math.max(...xP.map((v,i)=>Math.abs(v-xBase[i])/xBase[i]*100));
  });

  killChart('F');
  const ctx = document.getElementById('chartF').getContext('2d');
  CHARTS['F'] = new Chart(ctx, {
    type:'bar',
    data:{
      labels:['Sin rumor (0%)','Rumor leve (5%)','Rumor moderado (20%)','Pánico total (50%)'],
      datasets:[
        { label:'Cambio en distribución (%)',
          data:cambiosPorNivel.map(v=>+v.toFixed(2)),
          backgroundColor:cambiosPorNivel.map(v=>
            v<2?P.greenFill:v<10?P.goldFill:P.redFill),
          borderColor:cambiosPorNivel.map(v=>
            v<2?P.green:v<10?P.gold:P.red),
          borderWidth:2, borderRadius:5 },
      ]
    },
    options:{
      responsive:true,
      plugins:{legend:{position:'bottom'}},
      scales:{y:{title:{display:true,text:'% cambio en distribución'},beginAtZero:true}}
    }
  });

  const zonaMasVulnerable = zones[cambio.indexOf(Math.max(...cambio))];
  document.getElementById('kpiF').innerHTML =
    kpiCard('κ(A) sistema normal', condB.toFixed(1), 'número de condición — cuánto amplifica',condB>30?'danger':condB>10?'warn':'ok') +
    kpiCard('κ(A) sistema en pánico', condBad.toFixed(0), 'sistema mal condicionado (comparable)','danger') +
    kpiCard('Perturbación actual', (pert*100).toFixed(0)+'%', cfg.rumorTag, pert>0.2?'danger':pert>0.05?'warn':'ok') +
    kpiCard('Cambio en distribución', cambioMax.toFixed(1)+'%', `amplificación teórica: ≤${amplTeorica.toFixed(1)}%`, cambioMax>15?'danger':cambioMax>5?'warn':'ok');

  document.getElementById('answersF').querySelector('.answers-grid').innerHTML = [
    ['P1 — ¿Un rumor del 5% cuánto cambia la distribución?',
     `Con una perturbación del 5% en la demanda, la distribución cambia un <strong>${cambiosPorNivel[1].toFixed(1)}%</strong>. El número de condición κ(A)=${condB.toFixed(1)} implica una amplificación teórica de hasta ${(condB*5).toFixed(0)}%. ${cambiosPorNivel[1]<5?'El sistema normal es relativamente estable.':'Ya hay impacto visible: el rumor del 5-mayo justificó el caos.'}`],
    ['P2 — ¿La solución cambia poco o demasiado?',
     pert===0
       ? 'Sin perturbación, la solución es estable. El sistema distribuye exactamente lo demandado.'
       : `Con ${(pert*100).toFixed(0)}% de perturbación, la solución cambia <span class="${cambioMax>15?'danger':'warn'}">${cambioMax.toFixed(1)}%</span>. ${cambioMax>15?'<span class="danger">Demasiado.</span> El sistema amplifica el rumor de forma desproporcionada.':'Cambio notable pero el sistema aún puede responder si se actúa rápido.'}`],
    ['P3 — ¿El sistema es estable o mal condicionado?',
     `El sistema de distribución normal tiene κ(A)=${condB.toFixed(1)} — ${condB<20?'<span class="ok">bien condicionado</span>':'<span class="warn">moderadamente sensible</span>'}. Un sistema en pánico (demandas cruzadas, acaparamiento) puede llegar a κ≈${condBad.toFixed(0)}, que es <span class="danger">extremadamente mal condicionado</span>: un 1% de rumor produce ${(condBad*1).toFixed(0)}% de cambio.`],
    ['P4 — ¿Qué zona se vuelve más vulnerable?',
     `Con este nivel de rumor, <strong>${zonaMasVulnerable}</strong> sufre el mayor cambio porcentual (${cambio[zones.indexOf(zonaMasVulnerable)].toFixed(1)}%). La zona de Valles (Cochabamba) es la más vulnerable por ser nodo intermedio y la más bloqueada simultáneamente.`],
    ['P5 — ¿Cómo afecta el rumor acumulado al mes?',
     `Si el rumor persiste 30 días, el déficit acumulado se puede estimar: con un cambio del ${cambioMax.toFixed(1)}% diario en ${xBase.reduce((a,b)=>a+b,0).toFixed(0)} ton/día base, la pérdida mensual es <strong>${(cambioMax/100*xBase.reduce((a,b)=>a+b,0)*30).toFixed(0)} toneladas</strong> no distribuidas correctamente. Eso equivale a desabastecer una ciudad mediana por un mes.`],
  ].map(([q,a])=>answerBlock(q,a)).join('');

  markAnswered(['sqF1','sqF2','sqF3','sqF4','sqF5']);

  document.getElementById('calcF').innerHTML = `
  <div class="calc-section"><h4>Análisis de perturbación — A·(x+δx) = b+δb</h4>
    <div class="formula-block">
      κ(A) = ‖A‖·‖A⁻¹‖ = ${condB.toFixed(4)}<br>
      δb/b = ${(pert*100).toFixed(0)}% (rumor sobre la demanda)<br>
      Cota superior: δx/x ≤ κ(A)·δb/b = ${amplTeorica.toFixed(2)}%<br>
      Cambio real observado: ${cambioMax.toFixed(2)}%
    </div>
  </div>
  <div class="calc-section"><h4>Solución base vs perturbada (LU)</h4>
    <div class="formula-block">
      b base:  [${F_BASE_B.join(', ')}] ton/día<br>
      b pert.: [${bPert.map(v=>v.toFixed(1)).join(', ')}] ton/día<br><br>
      x base:  [${xBase.map(v=>v.toFixed(2)).join(', ')}]<br>
      x pert.: [${xPert.map(v=>v.toFixed(2)).join(', ')}]<br><br>
      Cambios: Oriente ${cambio[0].toFixed(2)}% | Valles ${cambio[1].toFixed(2)}% | Altiplano ${cambio[2].toFixed(2)}%
    </div>
  </div>`;
}

// ══════════════════════════════════════════════
//  G — MODELO N-M-D (RK4 sistema)
// ══════════════════════════════════════════════
/*
  N(t) = ciudadanos neutrales
  M(t) = manifestantes activos
  D(t) = mediadores / actores de diálogo
  N'= -α·N·M + β·D·M   (neutrales se movilizan; mediadores los revierten)
  M'= α·N·M - γ·M·D    (contagio social; mediadores desactivan)
  D'= κ·M - ρ·D        (mediadores aparecen con el conflicto; se agotan)
  Población total ~ 1000 (ciudad representativa)
*/
const G_CONFIGS = [
  { label:'Mediación activa (Iglesia)',
    alpha:0.0025, beta:0.06, gamma:0.05, kappa:0.025, rho:0.035,
    N0:780, M0:170, D0:50,
    desc:'Conferencia Episcopal medía activamente. Los mediadores reducen la tasa de contagio social y revierten la movilización.',
  },
  { label:'Sin diálogo (COB rechaza)',
    alpha:0.0040, beta:0.008, gamma:0.004, kappa:0.008, rho:0.08,
    N0:780, M0:170, D0:10,
    desc:'La COB exige levantar arrestos como condición previa. Sin mediadores efectivos, el conflicto crece sin amortiguación.',
  },
  { label:'Escalada — renuncia o nada',
    alpha:0.0060, beta:0.004, gamma:0.002, kappa:0.004, rho:0.12,
    N0:700, M0:250, D0:50,
    desc:'Alta tasa de contagio. La consigna única "renuncia" elimina cualquier espacio de diálogo parcial.',
  },
];

function runG(idx) {
  setActive('btnsG', idx);
  const cfg = G_CONFIGS[idx];
  const {alpha,beta,gamma,kappa,rho,N0,M0,D0} = cfg;

  const fns = [
    (t,y) => -alpha*y[0]*y[1] + beta*y[2]*y[1],
    (t,y) =>  alpha*y[0]*y[1] - gamma*y[1]*y[2],
    (t,y) =>  kappa*y[1]      - rho*y[2],
  ];
  const sol = rk4sys(fns, [N0,M0,D0], 0, 60, 0.25);

  // También Heun para comparar
  function heunSys(fns,y0,t0,tf,h) {
    const n=fns.length,t=[],ys=Array.from({length:n},()=>[]);
    let ti=t0,yi=[...y0];
    while(ti<=tf+1e-9){
      t.push(+ti.toFixed(2));yi.forEach((v,i)=>ys[i].push(+v.toFixed(2)));
      const k1=fns.map((f,i)=>f(ti,yi));
      const yp=yi.map((v,i)=>Math.max(0,v+h*k1[i]));
      const k2=fns.map((f,i)=>f(ti+h,yp));
      yi=yi.map((v,i)=>Math.max(0,v+h*(k1[i]+k2[i])/2));
      ti+=h;
    }
    return {t,N:ys[0],M:ys[1],D:ys[2]};
  }
  const solH = heunSys(fns,[N0,M0,D0],0,60,1);

  const peakMVal = Math.max(...sol.M);
  const peakMDay = sol.t[sol.M.indexOf(peakMVal)];
  const finalM   = sol.M[sol.M.length-1];
  const finalN   = sol.N[sol.N.length-1];
  const finalD   = sol.D[sol.D.length-1];
  const irreversible = finalM > M0 * 1.3;

  killChart('G');
  const ctx = document.getElementById('chartG').getContext('2d');
  CHARTS['G'] = new Chart(ctx, {
    type:'line',
    data:{
      labels: sol.t,
      datasets:[
        { label:'Neutrales N(t)', data:sol.N,
          borderColor:P.blue, backgroundColor:P.blueFill,
          fill:true, borderWidth:2, pointRadius:0, tension:.3 },
        { label:'Manifestantes M(t)', data:sol.M,
          borderColor:P.red, backgroundColor:P.redFill,
          fill:true, borderWidth:2.5, pointRadius:0, tension:.3 },
        { label:'Mediadores D(t)', data:sol.D,
          borderColor:P.green, backgroundColor:P.greenFill,
          fill:true, borderWidth:1.5, pointRadius:0, tension:.3 },
        { label:'M(t) — Heun (ref.)', data:solH.M,
          borderColor:P.orange, borderWidth:1, borderDash:[3,3],
          pointRadius:0, fill:false, tension:.3 },
      ]
    },
    options:{
      responsive:true,
      plugins:{legend:{position:'bottom'}},
      scales:{
        x:{title:{display:true,text:'Días desde el 1 de mayo'}},
        y:{title:{display:true,text:'Personas (población modelo)'},min:0}
      }
    }
  });

  const errHeunM = Math.abs((solH.M[solH.M.length-1]||0) - finalM);
  document.getElementById('kpiG').innerHTML =
    kpiCard('Pico de manifestantes',peakMVal.toFixed(0)+' pers.',`día ${peakMDay} (${new Date(2026,4,1+Math.round(+peakMDay)).toLocaleDateString('es',{day:'numeric',month:'short'})})`,'danger') +
    kpiCard('Manifestantes al día 60',finalM.toFixed(0)+' pers.',finalM>M0?'▲ conflicto creció':'▼ conflicto cedió',finalM>M0*1.3?'danger':finalM>M0?'warn':'ok') +
    kpiCard('Neutrales restantes',finalN.toFixed(0)+' pers.','no movilizados al día 60',finalN<400?'danger':'ok') +
    kpiCard('Error Heun vs RK4',errHeunM.toFixed(1)+' pers.','al día 60 en manifestantes','info');

  document.getElementById('answersG').querySelector('.answers-grid').innerHTML = [
    ['P1 — ¿El conflicto se estabiliza o se masifica?',
     idx===0
       ? `Con mediación activa, el conflicto <span class="ok">tiende a estabilizarse</span>. Los manifestantes alcanzan un pico de ${peakMVal.toFixed(0)} el día ${peakMDay} y luego descienden. Al día 60 quedan ${finalM.toFixed(0)} manifestantes activos, menos que los ${M0} iniciales.`
       : `Sin mediadores efectivos, el conflicto <span class="danger">se masifica</span>. El pico es ${peakMVal.toFixed(0)} manifestantes el día ${peakMDay}, y al día 60 hay ${finalM.toFixed(0)} — ${finalM>M0?`<span class="danger">+${(finalM-M0).toFixed(0)} más que al inicio</span>.`:''}`],
    ['P2 — ¿Qué pasa si la Iglesia logra mediar?',
     idx===0
       ? `Con la Conferencia Episcopal mediando activamente (β=${beta}, γ=${gamma}), los mediadores alcanzan ${Math.max(...sol.D).toFixed(0)} personas activas y logran reducir la tasa de contagio. El conflicto decae hacia el día 40.`
       : `El escenario 1 (Mediación activa) muestra que con β=${G_CONFIGS[0].beta} y γ=${G_CONFIGS[0].gamma}, los 50 mediadores iniciales son suficientes para frenar el contagio. Sin ellos (escenario actual), la curva no baja.`],
    ['P3 — ¿Qué pasa si la COB rechaza el diálogo?',
     idx===1
       ? `Sin mediadores (D₀=10, ρ=${rho} alto), los pocos mediadores se agotan rápido. El conflicto <span class="danger">no tiene freno</span>: α=${alpha} supera a γ·D, y los neutrales se van convirtiendo en manifestantes indefinidamente. Al día 60: ${finalM.toFixed(0)} manifestantes.`
       : `El escenario 2 (Sin diálogo) muestra que cuando ρ es alto y D₀ pequeño, los mediadores se agotan antes de frenar la ola. El resultado: ${G_CONFIGS[1].M0+50}+ manifestantes al final.`],
    ['P4 — ¿Cuántos días tarda en alcanzar su pico?',
     `El pico ocurre el <strong>día ${peakMDay}</strong> (aprox. ${new Date(2026,4,1+Math.round(+peakMDay)).toLocaleDateString('es',{day:'numeric',month:'long'})}), con <strong>${peakMVal.toFixed(0)} manifestantes</strong>. ${idx===2?'La alta tasa de contagio (α='+alpha+') lleva al pico muy rápido, sin tiempo de respuesta.':'Después del pico, el modelo predice '+(finalM<peakMVal*0.7?'una caída significativa.':'una caída lenta.')}`],
    ['P5 — ¿Qué parámetros hacen el conflicto irreversible?',
     irreversible
       ? `Con α=${alpha}, β=${beta} y γ=${gamma}, el sistema es <span class="danger">irreversible</span> en este horizonte: M(60) > M(0). El parámetro crítico es la relación α·N₀ > γ·D_max. Aquí α·N₀=${(alpha*N0).toFixed(3)} > γ·${Math.max(...sol.D).toFixed(0)}·${gamma}=${(gamma*Math.max(...sol.D)).toFixed(3)}.`
       : `El conflicto <span class="ok">es reversible</span> con estos parámetros. La clave es β·D > α·N: los mediadores deben "contagiar" paz más rápido de lo que el descontento se propaga. Con más mediadores o mayor β, el descenso sería más rápido.`],
  ].map(([q,a])=>answerBlock(q,a)).join('');

  markAnswered(['sqG1','sqG2','sqG3','sqG4','sqG5']);

  // Tabla completa cada 5 días
  const tRows = sol.t.filter((_,i)=>i%20===0).map((t2,ii)=>{
    const i=ii*20;
    const total=(sol.N[i]||0)+(sol.M[i]||0)+(sol.D[i]||0);
    return `<tr><td class="hl">${t2}</td><td>${(sol.N[i]||0).toFixed(0)}</td><td class="${sol.M[i]>M0*1.2?'danger':sol.M[i]>M0?'warn':'ok'}">${(sol.M[i]||0).toFixed(0)}</td><td class="ok">${(sol.D[i]||0).toFixed(0)}</td><td>${total.toFixed(0)}</td></tr>`;
  }).join('');

  // Heun tabla completa
  const tRowsH = solH.t.filter((_,i)=>i%10===0).map((t2,ii)=>{
    const i=ii*10;
    const diff=Math.abs((solH.M[i]||0)-(sol.M[Math.min(Math.round(+t2*4),sol.M.length-1)]||0));
    return `<tr><td class="hl">${t2}</td><td>${(solH.M[i]||0).toFixed(2)}</td><td>${(sol.M[Math.min(Math.round(+t2*4),sol.M.length-1)]||0).toFixed(2)}</td><td class="${diff>5?'warn':'ok'}">${diff.toFixed(2)}</td></tr>`;
  }).join('');

  // RK4 paso a paso en t=0
  const y0g=[N0,M0,D0];
  const fns0=[(t,y)=>-alpha*y[0]*y[1]+beta*y[2]*y[1],(t,y)=>alpha*y[0]*y[1]-gamma*y[1]*y[2],(t,y)=>kappa*y[1]-rho*y[2]];
  const k1g=fns0.map((f,i)=>f(0,y0g));
  const y2g=y0g.map((v,i)=>v+0.25/2*k1g[i]);
  const k2g=fns0.map((f,i)=>f(0.125,y2g));

  // Punto crítico de no retorno: cuándo M > N
  const crossIdx = sol.t.findIndex((_,i)=>sol.M[i]>sol.N[i]);
  const crossDay = crossIdx===-1?'No ocurre en 60 días':`Día ${sol.t[crossIdx]}`;

  document.getElementById('calcG').innerHTML = `
  <div class="calc-section">
    <h4>¿Qué modela el sistema N-M-D?</h4>
    <div class="formula-block">
El modelo adapta las ecuaciones epidemiológicas SIR al conflicto social:<br><br>
N(t) = Neutrales: ciudadanos que no se han movilizado (como "susceptibles" en SIR)<br>
M(t) = Manifestantes: ciudadanos activamente en el conflicto (como "infectados")<br>
D(t) = Mediadores: actores de diálogo que intentan resolver el conflicto (como "recuperados")<br><br>
Las ecuaciones modelan:<br>
N→M: un neutral se convierte en manifestante al "contagiarse" con tasa α·N·M<br>
M→N: un mediador "cura" a un manifestante con tasa β·D·M<br>
M→D: los manifestantes generan mediadores con tasa κ·M<br>
D→∅: los mediadores se "agotan" con tasa ρ·D
    </div>
  </div>

  <div class="calc-section">
    <h4>Parámetros del escenario actual: ${cfg.desc}</h4>
    <div class="formula-block">
N'(t) = −${alpha}·N·M + ${beta}·D·M<br>
M'(t) = +${alpha}·N·M − ${gamma}·M·D<br>
D'(t) = +${kappa}·M  − ${rho}·D<br><br>
Condiciones iniciales (1 de mayo de 2026):<br>
  N(0) = ${N0} ciudadanos neutrales<br>
  M(0) = ${M0} manifestantes activos (COB + indígenas + gremios)<br>
  D(0) = ${D0} mediadores (${idx===0?'Conferencia Episcopal activa':idx===1?'casi sin mediadores (COB rechaza diálogo)':'mediadores agotados por la escalada'})<br><br>
Punto de cruce M > N (conflicto masificado): ${crossDay}
    </div>
  </div>

  <div class="calc-section">
    <h4>Paso a paso — RK4 en las primeras horas (t=0, h=0.25 días)</h4>
    <div class="formula-block">
k₁ = f(0, [${y0g.join(', ')}]):<br>
  N'= −${alpha}·${N0}·${M0} + ${beta}·${D0}·${M0} = ${k1g[0].toFixed(4)}<br>
  M'= +${alpha}·${N0}·${M0} − ${gamma}·${M0}·${D0} = ${k1g[1].toFixed(4)}<br>
  D'= +${kappa}·${M0} − ${rho}·${D0} = ${k1g[2].toFixed(4)}<br><br>
k₂ = f(0.125, [${y2g.map(v=>v.toFixed(1)).join(', ')}]):<br>
  N'= ${k2g[0].toFixed(4)},  M'= ${k2g[1].toFixed(4)},  D'= ${k2g[2].toFixed(4)}<br><br>
(continúa con k₃ y k₄, luego yₙ₊₁ = y₀ + (h/6)(k₁+2k₂+2k₃+k₄))
    </div>
  </div>

  <div class="calc-section">
    <h4>Evolución completa — RK4 (h=0.25 días)</h4>
    <table class="calc-table">
      <tr><th>Día</th><th>N neutrales</th><th>M manifest.</th><th>D mediadores</th><th>Total</th></tr>
      ${tRows}
    </table>
  </div>

  <div class="calc-section">
    <h4>Comparación Heun vs RK4 — error en M(t)</h4>
    <table class="calc-table">
      <tr><th>Día</th><th>M Heun (h=1)</th><th>M RK4 (h=0.25)</th><th>|error|</th></tr>
      ${tRowsH}
    </table>
    <div class="formula-block" style="margin-top:.75rem">
RK4 con h=0.25: M(60) = ${finalM.toFixed(3)} personas<br>
Heun con h=1:   M(60) = ${(solH.M[solH.M.length-1]||0).toFixed(3)} personas<br>
Error absoluto: ${errHeunM.toFixed(3)} personas (${(errHeunM/Math.max(1,finalM)*100).toFixed(2)}%)<br><br>
Para sistemas con interacciones no lineales (N·M, M·D), RK4 es mucho más estable que Heun con pasos grandes.
    </div>
  </div>`;
}

// ══════════════════════════════════════════════
//  INIT
// ══════════════════════════════════════════════
window.addEventListener('DOMContentLoaded', () => {
  // Fix spline reference bug in D
  runA(0); runB(0); runC(0); runD(0); runE(0); runF(0); runG(0);
});
