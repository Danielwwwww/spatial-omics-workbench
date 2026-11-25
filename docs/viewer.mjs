import {regionSummary, associations, tilePairs} from './statistics.mjs';
import {regionDisplay} from './region-names.mjs';
import {organAtlases} from './organ-atlases.mjs';

const el = id => document.getElementById(id);
const pct = v => `${(v*100).toFixed(1)}%`;
const val = n => Number(n).toLocaleString('zh-CN',{maximumSignificantDigits:4});
const intensityText = (mode,value) => mode==='norm'?`${val(value*100)}%`:val(value);
const count = n => Number(n).toLocaleString('zh-CN',{maximumFractionDigits:0});
const node = (tag, text, cls) => {
  const element=document.createElement(tag);
  if (text!==undefined) element.textContent=text;
  if (cls) element.className=cls;
  return element;
};
const palettes = {ion:['#47348e','#7342a5','#d25182','#fcac65','#fff077'],
                  gene:['#9af7ff','#49dcff','#04bfff','#008be5','#004fbd']};
let data=null, image=null, localCompute=false;
const canvasGeometry={};
const decodedGenes = new Map();
let source=null, draft=[], editingBusy=false;
let computeWorker;

function notify(msg='',error=false){el('message').textContent=msg;el('message').classList.toggle('error',error);}

function sourceLink(label, url){
  const link=node('a',label);
  link.href=url;link.target='_blank';link.rel='noopener noreferrer';
  return link;
}

function showOrganReference(container, reference, selectName){
  container.replaceChildren();
  container.hidden=!reference;
  if(!reference)return;
  const figure=node('figure'), img=node('img'), caption=node('figcaption',reference.caption);
  img.src=reference.image;img.alt=`${reference.title}参考切片`;img.loading='lazy';img.referrerPolicy='no-referrer';
  figure.append(img,caption);
  const content=node('div'), list=node('ul');
  content.append(node('strong',`${reference.title} · 组织分区参考`),
    node('p','依据参考图，在自己的 H&E 上手绘边界。边界坐标由实际切片确定。','sub'));
  for(const [english,chinese,hint] of reference.structures){
    const item=node('li'), label=`${english}（${chinese}）`;
    if(selectName){
      const button=node('button',label,'atlas-label');
      button.type='button';button.onclick=()=>selectName(english);
      item.append(button);
    }else item.append(node('strong',label));
    item.append(node('span',hint,'sub'));list.append(item);
  }
  content.append(list,sourceLink('查看 NTP 啮齿类组织图谱',reference.url));
  container.append(figure,content);
}

function referenceForResult(){
  return organAtlases.get(source?.config?.organ_reference);
}
function options(id,names,withAll){
  const s=el(id);s.replaceChildren();
  names.forEach((name,i)=>{if(i===0&&withAll===true)return;s.add(new Option(name,String(i)));});
  if(withAll)s.add(new Option('全部组织','0'),0);
  if(!s.length)s.add(new Option('当前无基因','-1'));
  s.selectedIndex=0;
}

function validate(bundle){
  if(bundle?.schema!==1||!bundle.image?.jpeg||!Array.isArray(bundle.ionNames)||!bundle.ionNames.length||!bundle.grid?.xy||!bundle.msi?.xy){
    throw new Error('配准结果缺少可显示的组织图像、离子或空间网格。');
  }
}

async function load(bundle, origin=null){
  validate(bundle);
  const next=new Image(); next.src=`data:image/jpeg;base64,${bundle.image.jpeg}`;
  await next.decode(); data=bundle; image=next; source=origin; decodedGenes.clear();
  options('ion',data.ionNames,false);options('region',data.regions.map((_,id)=>displayRegion(id).english),true);options('gene',data.geneNames,false);
  if(data.stMap?.length){
    const clearest=data.stMap.reduce((best,layer,i)=>layer.n>data.stMap[best].n?i:best,0);
    el('gene').value=String(clearest);
  }
  el('title').textContent=data.name;el('origin').textContent=data.origin;
  el('coverage').textContent=pct(data.qa.tissueCommonFraction);el('tiles').textContent=val(data.qa.overlapGridCells);
  el('features').textContent=`${data.ionNames.length} / ${data.geneNames.length}`;
  el('results').hidden=false;
  el('roi-open').hidden=!source;
  const reference=referenceForResult();
  el('region-name-header').textContent=data.atlas?'Brain region':'Region';
  el('roi-reference').hidden=!reference;
  showOrganReference(el('roi-atlas-reference'),reference,name=>{
    el('roi-name').value=name;drawROI();el('roi-name').focus();
  });
  const resolutions=Object.keys(data.stResolutions||{});
  el('st-resolution').replaceChildren(...resolutions.map(size=>new Option(`bin${size}`,size)));
  el('st-resolution-control').hidden=!resolutions.length;
  el('atlas-only').checked=false;
  showQuality();render();
}

function coordinates(canvas){
  const w=canvas.width,h=canvas.height;
  const scale=Math.min((w-20)/data.image.width,(h-20)/data.image.height);
  return {scale,dx:(w-data.image.width*scale)/2,dy:(h-data.image.height*scale)/2};
}

function regionColor(id){
  if(id===1)return '#aaa';
  return data.manualRegions?.some(r=>r.name===data.regions[id])?'#151515':'#ffe44e';
}

function regionClass(id){
  return id===1?'region-neutral':regionColor(id)==='#151515'?'region-manual':'region-atlas';
}

function displayRegion(id){
  const name=data.regions[id];
  const pair=referenceForResult()?.structures.find(([english])=>english===name);
  return pair?{english:pair[0],chinese:pair[1],short:pair[0]}:
    regionDisplay(name,regionClass(id)==='region-manual');
}

function drawRegions(ctx,geom,region,showAll){
  const {scale,dx,dy}=geom;
  const occupied=[];
  ctx.font=showAll?'bold 13px sans-serif':'bold 11px sans-serif';
  ctx.textAlign='center';ctx.textBaseline='middle';
  for(const r of data.boundaries){
    if(region&&r.id!==region)continue;
    ctx.fillStyle=regionColor(r.id);ctx.globalAlpha=.76;
    const width=showAll?1.6:1.2;
    for(const [x,y] of r.xy)ctx.fillRect(dx+x*scale,dy+y*scale,width,width);
  }
  const sorted=[...data.boundaries].sort((a,b)=>(b.width||0)-(a.width||0));
  ctx.globalAlpha=1;
  for(const r of sorted){
    if((region&&r.id!==region)||r.id===1||(!showAll&&!region&&(r.width||0)<13))continue;
    const anchor=r.labelXY??r.xy[Math.floor(r.xy.length/2)];
    if(!anchor)continue;
    const name=displayRegion(r.id).map,half=ctx.measureText(name).width/2+3;
    const x=Math.min(dx+data.image.width*scale-half,Math.max(dx+half,dx+anchor[0]*scale));
    const y=dy+anchor[1]*scale;
    const bounds=[x-half,y-8,x+half,y+8];
    if(!region&&occupied.some(([l,t,rr,b])=>bounds[0]<rr&&bounds[2]>l&&bounds[1]<b&&bounds[3]>t))continue;
    occupied.push(bounds);
    ctx.fillStyle='#fff';
    ctx.fillRect(bounds[0],bounds[1],bounds[2]-bounds[0],bounds[3]-bounds[1]);
    ctx.fillStyle=regionClass(r.id)==='region-manual'?'#151515':'#806000';
    ctx.fillText(name,x,y);
  }
}

function map(kind){
  const canvas=el(`${kind}-map`),ctx=canvas.getContext('2d');
  const geom=coordinates(canvas),{scale,dx,dy}=geom;canvasGeometry[kind]=geom;
  ctx.fillStyle='#f7f7fb';ctx.fillRect(0,0,canvas.width,canvas.height);
  ctx.drawImage(image,dx,dy,data.image.width*scale,data.image.height*scale);
  const ion=Number(el('ion').value), gene=Number(el('gene').value), region=Number(el('region').value);
  const mode=el('mode-ion').value, atlasOnly=el('atlas-only').checked;
  const alpha=Number(el('opacity').value);
  ctx.save();ctx.beginPath();ctx.rect(dx,dy,data.image.width*scale,data.image.height*scale);ctx.clip();
  const paint = (x,y,fraction,radius) => {
    ctx.fillStyle=palettes[kind][Math.min(4,Math.floor(fraction*5))];
    ctx.globalAlpha=kind==='gene'?alpha:alpha*(.21+.79*fraction);
    ctx.fillRect(dx+x*scale-radius,dy+y*scale-radius,radius*2,radius*2);
  };
  let top=0;
  if(!atlasOnly && kind==='ion'){
    const m=data.msi, values=m[mode][ion];
    const positive=values.filter((v,i)=>v>0&&m.region[i]>0).sort((a,b)=>a-b);
    top=positive.length?positive[Math.floor(.98*(positive.length-1))]:0;
    const dot=Math.max(1,Math.min(5,Math.sqrt(data.image.width*data.image.height/m.xy.length)*.48*scale));
    for(let i=0;i<m.xy.length;i+=1){
      const v=values[i],r=m.region[i];
      if(!(v>0)||!r||(region&&r!==region))continue;
      paint(...m.xy[i],Math.min(1,Math.sqrt(v/top)),dot);
    }
  }else if(!atlasOnly && kind==='gene'){
    const size=el('st-resolution').value, key=`${size}:${gene}`;
    const layer=(data.stResolutions?.[size]?.maps??data.stMap)?.[gene];
    if(layer?.n){
      if(!decodedGenes.has(key)){
        const raw=Uint8Array.from(atob(layer.points),char=>char.charCodeAt(0));
        decodedGenes.set(key,new DataView(raw.buffer));
      }
      const pixels=decodedGenes.get(key), dot=Math.max(.85,scale*.95);
      top=layer.p98;
      for(let i=0;i+6<=pixels.byteLength;i+=6){
        const roi=pixels.getUint8(i+5);
        if(region&&region!==roi)continue;
        const pos=pixels.getUint32(i,true);
        const stored=pixels.getUint8(i+4)/255;
        const level=layer.scale==='fixed-log1p-cp10k'?stored:stored*layer.p98/Math.log1p(10000);
        paint(pos%data.image.width,Math.floor(pos/data.image.width),level,dot);
      }
    }
  }
  if(atlasOnly||el('boundaries').checked)drawRegions(ctx,geom,region,atlasOnly);
  ctx.restore();
  el(`${kind}-max`).textContent=atlasOnly?'H&E 分区':top?(kind==='gene'?'每万 UMI · 对数色标':`P98 ${val(top)}`):
    kind==='gene'&&!data.stMap?'此配准包只有网格汇总':'无阳性信号';
}

function showRegion(ion,mode){
  const container=el('region-stats');container.replaceChildren();
  const rows=regionSummary(data,ion,mode);
  for(const r of rows){
    const row=node('tr'), header=node('td');
    const names=displayRegion(r.id),select=node('button',names.english);
    select.type='button';
    select.setAttribute('aria-pressed',String(Number(el('region').value)===r.id));
    select.onclick=()=>{el('region').value=String(r.id);render();};
    select.classList.add(regionClass(r.id));
    header.append(select);
    row.append(header,node('td',names.chinese?`（${names.chinese}）`:'—'),
               node('td',intensityText(mode,r.mean)),node('td',count(r.n)),node('td',pct(r.detect)));
    container.append(row);
  }
  if(!rows.length){
    const row=node('tr'), cell=node('td','所选离子没有组织内测量。');
    cell.colSpan=5;row.append(cell);container.append(row);
  }
  const note=el('tissue-note');
  note.replaceChildren();
  if(data.atlas){
    note.append(data.atlas.method+' ');
    if(data.atlas.ap_mm_from_ccf_origin!==undefined){
      const a=data.atlas;
      const named=a.named_fraction??a.atlas_assigned_he_fraction;
      const rotation=a.rotation===undefined?'':` · 旋转 ${a.rotation}°`;
      const matched=named===undefined?'待复核':`图谱命名覆盖 ${pct(named)}`;
      note.append(node('span',`冠状层面 ${val(a.ap_mm_from_ccf_origin)} mm（CCF 原点）${rotation} · ${matched}。`),' ');
    }
    const source=node('a','Allen 冠状参考图谱');source.href=data.atlas.source;
    source.target='_blank';source.rel='noopener';note.append(source);
  }else{
    const reference=referenceForResult();
    note.textContent=reference?`${reference.title}参考图谱支持对照手绘，统计按切片上的手绘边界计算。`:
      data.regions.length===2?'可通过手绘添加组织分区。':'区域边界来自 H&E 标注。';
  }
}

function showMetabolites(region,mode,selected){
  const body=el('metabolites');
  body.replaceChildren();
  el('ion-table-note').textContent=`${region?displayRegion(region).english:'全部组织'} · ${mode==='norm'?'离子 / TIC':'原始强度'}`;
  const evidence=el('ion-evidence');
  evidence.replaceChildren();
  const selectedAnnotation=data.ionAnnotations?.[data.ionNames[selected]];
  if(selectedAnnotation){
    evidence.append(`${selectedAnnotation.name}${selectedAnnotation.ion?' · '+selectedAnnotation.ion:''}。`);
    if(Number.isFinite(selectedAnnotation.observedMedian))evidence.append(`窗内峰质心的中位 m/z ${selectedAnnotation.observedMedian.toFixed(6)}。 `);
    if(selectedAnnotation.evidence)evidence.append(selectedAnnotation.evidence+' ');
    if(typeof selectedAnnotation.source==='string'&&selectedAnnotation.source.startsWith('https://')){
      const link=node('a','参考证据');link.href=selectedAnnotation.source;
      link.target='_blank';link.rel='noopener';evidence.append(link);
    }
  }else evidence.textContent='本离子暂无候选注释。';
  data.ionNames.forEach((name, index)=>{
    const areas=regionSummary(data,index,mode);
    const chosen=region?areas.filter(a=>a.id===region):areas;
    const n=chosen.reduce((sum,r)=>sum+r.n,0);
    const total=chosen.reduce((sum,r)=>sum+r.sum,0);
    const detected=chosen.reduce((sum,r)=>sum+r.detected,0);
    const row=node('tr'), label=node('td'), button=node('button',name);
    button.type='button';
    button.setAttribute('aria-pressed',String(index===selected));
    button.onclick=()=>{
      el('ion').value=String(index);el('atlas-only').checked=false;render();
      el('ion-map').scrollIntoView({block:'nearest'});
    };
    label.append(button);
    row.append(label,node('td',data.ionAnnotations?.[name]?.name??'无候选注释'),node('td',n?intensityText(mode,total/n):'无测量'),
               node('td',n?pct(detected/n):'—'),node('td',count(n)),node('td',areas[0]?displayRegion(areas[0].id).english:'—'));
    body.append(row);
  });
}

function showExpression(){
  const summary = data.regionExpression;
  const gene = Number(el('gene').value), region = Number(el('region').value);
  const rows = el('expression-rows');
  rows.replaceChildren();
  el('expression').hidden = !summary || gene < 0;
  if (el('expression').hidden) return;
  el('expression-label').textContent = `${data.geneNames[gene]} · 每万 UMI 基因计数 · 三模态共同覆盖区${data.stResolutions?' · bin10':''}`;
  const selected=summary.filter(row => !region || row.id === region);
  if(!selected.length)rows.append(node('p','此区无通过质控的共同覆盖 ST 点。','sub'));
  for (const item of selected) {
    const caption = node('div',undefined,'region-caption');
    const names=displayRegion(item.id);
    caption.append(node('span',names.chinese?`${names.english}（${names.chinese}）`:names.english,regionClass(item.id)),node('strong',val(item.cp10k[gene])));
    rows.append(caption,node('p',`${count(item.stN)} 个 bin · ${count(item.umi)} UMI`,'sub'));
  }
}

function showCorrelations(ion,region){
  const body=el('correlations');body.replaceChildren();
  const method=el('correlation-method').value;
  const {n,rows}=associations(data,ion,region,8,method);
  const mapping=data.qa.gridMapping==='footprint'?'像素重叠加权':'按像素中心分格';
  el('correlation-note').textContent=`${n} 个共同网格 · ${mapping}`;
  for(const row of rows.slice(0,30)){
    const tr=node('tr'), name=node('td'), button=node('button',row.gene);button.type='button';
    button.onclick=()=>{el('gene').value=String(row.index);map('gene');showExpression();showGridPairs();showPairPlot();};name.append(button);
    tr.append(name,node('td',row.rho.toFixed(3)),node('td',String(n)));body.append(tr);
  }
  if(!rows.length){
    const tr=node('tr'); const msg=node('td',n<8?'至少需要 8 个共同网格。':data.geneNames.length?'没有可计算的变异信号。':'该 ST 输入未包含基因矩阵。');
    msg.colSpan=3;tr.append(msg);body.append(tr);
  }
}

function showPairPlot(){
  const ion=Number(el('ion').value), region=Number(el('region').value), gene=Number(el('gene').value);
  const canvas=el('pair-plot'), caption=el('pair-caption');
  const pairs=gene<0?[]:tilePairs(data,ion,region,gene);
  el('pair-panel').hidden=!pairs.length;
  if(!pairs.length)return;
  caption.textContent=`${data.ionNames[ion]} / ${data.geneNames[gene]} · 每个点对应一个共同网格`;
  const ctx=canvas.getContext('2d'),w=canvas.width,h=canvas.height,pad=44;
  const right=w-20, top=18, bottom=h-36;
  const maxX=Math.max(...pairs.map(p=>p.ion),1e-12),maxY=Math.max(...pairs.map(p=>p.log),1e-12);
  ctx.clearRect(0,0,w,h);
  ctx.strokeStyle='#aab6cd';ctx.beginPath();ctx.moveTo(pad,top);ctx.lineTo(pad,bottom);ctx.lineTo(right,bottom);ctx.stroke();
  ctx.fillStyle='#526276';ctx.font='12px sans-serif';
  ctx.fillText('ln(1 + CP10k)',pad+4,top-4);
  ctx.fillText('离子 / TIC',w-90,h-8);
  ctx.fillText('0',pad-10,bottom+15);
  if(maxX>1e-12)ctx.fillText(intensityText('norm',maxX),right-45,bottom+15);
  if(maxY>1e-12)ctx.fillText(val(maxY),3,top+7);
  ctx.fillStyle='#294fb79e';
  for(const pair of pairs){
    ctx.beginPath();ctx.arc(pad+pair.ion/maxX*(right-pad),bottom-pair.log/maxY*(bottom-top),3,0,2*Math.PI);ctx.fill();
  }
}

function showGridPairs(){
  if(!data||!el('grid-detail').open)return;
  const body=el('grid-pairs'),region=Number(el('region').value),ion=Number(el('ion').value);
  const gene=Number(el('gene').value);
  body.replaceChildren();
  if(gene<0){el('grid-detail-note').textContent='需要 ST 基因计数才能查看成对数据。';return;}
  const pairs=tilePairs(data,ion,region,gene);
  el('grid-detail-note').textContent=`${data.ionNames[ion]} 与 ${data.geneNames[gene]} · ${pairs.length} 个格子`;
  for(const entry of pairs){
    const row=node('tr');
    row.append(node('td',entry.xy?entry.xy.join(', '):'—'),node('td',count(entry.msi)),
               node('td',count(entry.st)),node('td',intensityText('norm',entry.ion)),
               node('td',val(entry.cp10k)),node('td',val(entry.log)));
    body.append(row);
  }
}

function showQuality(){
  const out=el('quality'),grid=node('div',undefined,'qa-grid');out.replaceChildren(grid);
  const labels={affine:'仿射',pilot_affine:'初始仿射',similarity:'相似变换',nonlinear:'平滑非线性',
                image_refined:'图像非线性细化',stalign:'STalign LDDMM'};
  for(const name of ['msi','st']){
    const qa=data.qa.alignment[name], result=qa.models[qa.selected];
    const item=node('div',undefined,'qa-item');item.append(node('b',name==='msi'?'代谢组 → H&E':'转录组 → H&E'));
    item.append(node('p',`当前图层使用 ${labels[qa.selected] ?? qa.selected}`));
    if(Number.isFinite(result.dice))item.append(node('p',`Dice ${result.dice.toFixed(3)}`));
    const baseline=qa.models.affine;
    if(!Number.isFinite(result.mean_tre))item.append(node('p','初步对齐 · 需用解剖点复核','warning'));
    else{
      item.append(node('p',`平均 / 最大 TRE ${val(result.mean_tre)} / ${val(result.max_tre)} ${qa.tre_units}`),node('p',`${qa.check_pairs} 个留出点`,'sub'));
      if(Number.isFinite(baseline?.max_tre)&&result.max_tre>baseline.max_tre)item.append(node('p','局部最大误差增加','warning'));
    }
    if(qa.stalign_refinement?.reason)item.append(node('p','STalign 形变未通过数值检查。','sub'));
    if(qa.stalign_refinement&&qa.tre_priority==='max')item.append(node('p','模型选择优先比较最差点误差。','sub'));
    const models=Object.entries(qa.models);
    if(models.length>1){
      const list=node('table',undefined,'qa-models');
      list.setAttribute('aria-label','各模型的留出点平均和最大目标配准误差');
      const header=node('tr');
      for(const text of ['模型','平均 TRE','最大 TRE'])header.append(node('th',text));
      list.append(header);
      for(const [key,model] of models){
        const row=node('tr');
        const format=x=>Number.isFinite(x)?val(x):'—';
        row.append(node('td',labels[key]??key),node('td',format(model.mean_tre)),node('td',format(model.max_tre)));
        list.append(row);
      }
      item.append(list);
    }
    if(qa.orientation_ambiguity)item.append(node('p','多个方向的轮廓重合接近','warning'));
    grid.append(item);
  }
}

function render(){
  if(!data)return;
  const ion=Number(el('ion').value),region=Number(el('region').value),mode=el('mode-ion').value;
  el('ion-label').textContent=data.ionNames[ion]+(data.qa.autoMz?' · 自动候选':'');
  el('map-note').textContent=`${count(data.qa.msiPixels)} 个 MSI 像素 · ${count(data.qa.stPoints)} 个 ST 点`;
  el('st-map-note').textContent=data.stMap?'各分辨率共用蓝青色标。未检出位置显示 H&E 底图。':'此配准包仅含共同网格统计。';
  const depth=data.stResolutions?.[el('st-resolution').value]?.depth;
  el('st-depth').textContent=depth?`裁剪后中位数 ${count(depth.umiQuartiles[1])} UMI / ${count(depth.geneQuartiles[1])} 个基因每 bin`:'';
  map('ion');map('gene');showRegion(ion,mode);showMetabolites(region,mode,ion);showCorrelations(ion,region);showExpression();showGridPairs();showPairPlot();
}

for(const id of ['ion','region','mode-ion'])el(id).addEventListener('change',render);
for(const id of ['gene','boundaries','opacity','atlas-only'])el(id).addEventListener('input',()=>{if(data){map('ion');map('gene');}});
el('gene').addEventListener('change',showExpression);
el('gene').addEventListener('change',showGridPairs);
el('gene').addEventListener('change',showPairPlot);
el('correlation-method').addEventListener('change',render);
function updateAtlasChoice(){
  const choice=el('atlas').value, active=choice==='allen-ccfv3';
  el('atlas-ap').disabled=!active;
  el('atlas-rotation').disabled=!active;
  el('atlas-input-note').hidden=!active;
  showOrganReference(el('atlas-reference'),organAtlases.get(choice));
}
el('atlas').addEventListener('change',updateAtlasChoice);
updateAtlasChoice();
el('grid-detail').addEventListener('toggle',showGridPairs);
el('st-resolution').addEventListener('change',render);
for(const kind of ['ion','gene'])el(`${kind}-map`).addEventListener('pointermove',event=>{
  if(!data)return;
  const canvas=event.currentTarget,r=canvas.getBoundingClientRect(),{scale,dx,dy}=canvasGeometry[kind];
  const x=((event.clientX-r.left)*canvas.width/r.width-dx)/scale;
  const y=((event.clientY-r.top)*canvas.height/r.height-dy)/scale;
  if(x<0||y<0||x>data.image.width||y>data.image.height)return;
  const point=toOriginal(x,y);
  el('map-hover').textContent=`H&E 原图坐标 x ${Math.round(point[0])}，y ${Math.round(point[1])}`;
});

function toOriginal(x,y){
  const [w,h]=data.image.originalSize;
  return data.image.toOriginal?.slice(0,2).map(([a,b,c])=>a*x+b*y+c) ??
    [(x+.5)*w/data.image.width-.5,(y+.5)*h/data.image.height-.5];
}

function drawROI(){
  const canvas=el('roi-canvas'),ctx=canvas.getContext('2d');
  const {scale,dx,dy}=coordinates(canvas);
  ctx.fillStyle='#f7f7fb';ctx.fillRect(0,0,canvas.width,canvas.height);
  ctx.drawImage(image,dx,dy,data.image.width*scale,data.image.height*scale);
  drawRegions(ctx,{scale,dx,dy},0,true);
  if(draft.length){
    ctx.beginPath();
    draft.forEach(([x,y],i)=>ctx[i?'lineTo':'moveTo'](dx+x*scale,dy+y*scale));
    if(draft.length>2){
      ctx.closePath();ctx.fillStyle='#ffe67244';ctx.fill();
    }
    ctx.strokeStyle='#1a1a1a';ctx.lineWidth=3;ctx.stroke();
    ctx.fillStyle='#fff8ab';
    for(const [x,y] of draft){
      ctx.beginPath();ctx.arc(dx+x*scale,dy+y*scale,4,0,2*Math.PI);ctx.fill();ctx.stroke();
    }
  }
  el('roi-undo').disabled=editingBusy||!draft.length;
  el('roi-clear').disabled=editingBusy||!draft.length;
  el('roi-save').disabled=editingBusy||draft.length<3||!el('roi-name').value.trim();
  el('roi-close').disabled=editingBusy;
}

async function compute(files, config, onProgress) {
  if (localCompute) {
    const body = new FormData();
    Object.entries(files).forEach(([key,file]) => body.set(key, file));
    body.set('config', JSON.stringify(config));
    const response = await fetch('./api/analyze', {method:'POST', body});
    const result = await response.json();
    if (!response.ok) throw new Error(typeof result.detail === 'string' ? result.detail : '计算失败，请核对文件格式。');
    return result;
  }
  if (!globalThis.Worker) throw new Error('当前浏览器缺少后台计算支持，请使用新版浏览器或运行本地服务。');
  computeWorker ??= new Worker(new URL('./compute/worker.js', import.meta.url));
  return new Promise((resolve, reject) => {
    computeWorker.onmessage = ({data: reply}) => {
      if (reply.message) onProgress(reply.message);
      else if (reply.error) reject(new Error(reply.error));
      else resolve(reply.result);
    };
    computeWorker.onerror = () => {
      computeWorker.terminate();
      computeWorker = undefined;
      reject(new Error('浏览器计算中断。请刷新重试，或使用本地计算服务处理较大的文件。'));
    };
    computeWorker.postMessage({files, config});
  });
}

function renderSavedRegions(){
  const saved=el('roi-saved');
  saved.replaceChildren();
  for(const roi of data.manualRegions||[]){
    const item=node('div',undefined,'roi-saved-row'),button=node('button','移除');
    button.type='button';button.className='btn';button.disabled=editingBusy;
    button.setAttribute('aria-label',`移除手绘区 ${roi.name}`);
    button.onclick=()=>applyRegions(data.manualRegions.filter(r=>r.name!==roi.name));
    item.append(node('span',roi.name),button);saved.append(item);
  }
  el('roi-export').disabled=editingBusy||!(data.manualRegions||[]).length;
}

async function applyRegions(regions, selected){
  if(editingBusy)return;
  editingBusy=true;
  const active=source, ion=el('ion').value, gene=el('gene').value, size=el('st-resolution').value;
  el('roi-status').textContent='正在重算区域的离子、表达和共同网格。';
  drawROI();renderSavedRegions();
  try{
    const config = {...active.config, manual_regions:regions};
    const result = await compute(active.files, config, message => {el('roi-status').textContent=message;});
    await load(result,{files:active.files,config});
    if(el('ion').options[Number(ion)])el('ion').value=ion;
    if(el('gene').options[Number(gene)])el('gene').value=gene;
    if(data.stResolutions?.[size])el('st-resolution').value=size;
    if(selected)el('region').value=String(data.regions.indexOf(selected));
    render();
    draft=[];el('roi-name').value='';
    el('roi-status').textContent=selected?`已保存 ${selected} 并更新统计。`:'已移除手绘区并更新统计。';
    if(selected){el('roi-dialog').close();notify(`已保存 ${selected} 并更新统计。`);}
  }catch(error){
    el('roi-status').textContent=error.message;
  }finally{
    editingBusy=false;drawROI();renderSavedRegions();
  }
}

el('roi-open').onclick=()=>{
  draft=[];el('roi-name').value='';el('roi-status').textContent='';
  el('roi-dialog').showModal();drawROI();renderSavedRegions();
};
el('roi-close').onclick=()=>el('roi-dialog').close();
el('roi-dialog').addEventListener('cancel',event=>{if(editingBusy)event.preventDefault();});
el('roi-name').addEventListener('input',drawROI);
el('roi-canvas').addEventListener('click',event=>{
  if(editingBusy||draft.length>=512)return;
  const canvas=event.currentTarget,r=canvas.getBoundingClientRect(),{scale,dx,dy}=coordinates(canvas);
  const x=((event.clientX-r.left)*canvas.width/r.width-dx)/scale;
  const y=((event.clientY-r.top)*canvas.height/r.height-dy)/scale;
  if(x<0||y<0||x>=data.image.width||y>=data.image.height)return;
  draft.push([x,y]);el('roi-status').textContent=`已标 ${draft.length} 个顶点。`;drawROI();
});
el('roi-undo').onclick=()=>{draft.pop();drawROI();};
el('roi-clear').onclick=()=>{draft=[];drawROI();};
el('roi-save').onclick=()=>{
  const name=el('roi-name').value.trim();
  if(data.regions.includes(name)){el('roi-status').textContent='这个名称已使用。';return;}
  const points=draft.map(([x,y])=>toOriginal(x,y).map((value,axis)=>
    Math.max(0,Math.min(data.image.originalSize[axis]-1,Math.round(value)))));
  applyRegions([...(data.manualRegions||[]),{name,points}],name);
};
el('roi-export').onclick=()=>{
  const link=node('a');
  const url=URL.createObjectURL(new Blob([JSON.stringify({manual_regions:data.manualRegions},null,2)],
    {type:'application/json'}));
  link.href=url;link.download='manual-regions.json';link.click();
  setTimeout(()=>URL.revokeObjectURL(url),1000);
};

function dataFiles(form){
  const msi=[...el('msi-files').files];
  const ibd=msi.find(file=>/\.ibd$/i.test(file.name));
  const msiMain=msi.find(file=>file!==ibd);
  if(msi.length>2||!msiMain||msi.length!==(ibd?2:1)||/\.imzml$/i.test(msiMain.name)!==Boolean(ibd))
    throw new Error('空间代谢组请选择一个强度表，或一对 imzML 和 IBD。');
  const st=[...el('st-files').files];
  const summary=st.find(file=>/summary/i.test(file.name));
  if(!st.length||st.length>2||(st.length===2&&!summary))
    throw new Error('空间转录组请选择一个 H5AD、坐标计数表，或一对 summary 和 DGE。');
  const files={he:form.elements.he.files[0],msi:msiMain,st:summary||st[0]};
  if(ibd)files.ibd=ibd;
  if(st.length===2)files.dge=st.find(file=>file!==summary);
  const mask=form.elements.he_mask.files[0];
  if(mask)files.he_mask=mask;
  return files;
}

const exampleButton=el('reference-example');
const submitButton=el('upload-form').querySelector('button[type=submit]');
exampleButton.onclick=async()=>{
  if(exampleButton.disabled||submitButton.disabled||editingBusy)return;
  exampleButton.disabled=true;submitButton.disabled=true;
  notify('正在加载参考样例。');
  try{
    const response=await fetch(new URL('./examples/brain-section.json',import.meta.url));
    if(!response.ok)throw new Error('参考样例暂不可用。');
    await load(await response.json());
    const annotated=data.ionNames.findIndex(name=>data.ionAnnotations?.[name]?.source);
    if(annotated>=0){el('ion').value=String(annotated);render();}
    el('upload').open=false;
    el('results').scrollIntoView({block:'start'});
    notify('');
  }catch(error){notify(error.message,true);}
  finally{exampleButton.disabled=false;submitButton.disabled=false;}
};

el('upload-form').onsubmit=async event=>{
  event.preventDefault();
  const form=event.currentTarget, submit=form.querySelector('button[type=submit]');
  if(submit.disabled||editingBusy)return;
  submit.disabled=true; exampleButton.disabled=true; el('roi-open').disabled=true; notify('正在读取文件和计算配准。大文件可能需要数分钟。');
  try{
    const files=dataFiles(form);
    const cfgFile=el('settings').files[0]; const cfg=cfgFile?JSON.parse(await cfgFile.text()):{};
    const title=el('dataset-name').value.trim(),mz=el('masses').value.trim(),genes=el('genes').value.trim();
    if(title)cfg.name=title;
    if(mz){cfg.mz=mz.split(/[,，\s]+/).map(Number);if(cfg.mz.some(x=>!Number.isFinite(x)||x<=0))throw new Error('请填写有效的正数 m/z');}
    if(genes)cfg.genes=genes.split(/[,，\s]+/);
    cfg.refine=el('refine').checked;
    const atlas=el('atlas').value;
    cfg.atlas=atlas==='allen-ccfv3'?atlas:'none';
    cfg.organ_reference=organAtlases.has(atlas)?atlas:undefined;
    const atlasAP=el('atlas-ap').value;
    if(cfg.atlas==='allen-ccfv3'&&atlasAP)cfg.atlas_ap_mm=Number(atlasAP);
    const atlasRotation=el('atlas-rotation').value;
    if(cfg.atlas==='allen-ccfv3'&&atlasRotation!=='auto')cfg.atlas_rotation=atlasRotation;
    cfg.stalign=localCompute && !el('stalign-control').hidden && el('stalign').checked;
    if(cfg.stalign)cfg.tre_priority=el('tre-priority').value;
    const result=await compute(files,cfg,notify);
    await load(result,{files,config:cfg}); el('upload').open=false;notify('配准已完成。请先复核质控指标。');
  }catch(error){notify(error.message,true);}finally{submit.disabled=false;exampleButton.disabled=false;el('roi-open').disabled=false;}
};

if (['localhost','127.0.0.1','[::1]'].includes(location.hostname)) {
  fetch('./api/health').then(r=>r.ok?r.json():null).then(result=>{
    localCompute=result?.upload===true;
    if(localCompute){
      el('mode').textContent='本地 Python 计算';
      el('stalign-control').hidden=result?.stalign!==true;
      el('compute-info').textContent='直接选择三类原始文件。配准由本机 Python 服务执行。';
    }
  }).catch(()=>{});
}
