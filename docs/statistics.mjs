export function ranks(values) {
  const order = values.map((_, i) => i).sort((a, b) => values[a] - values[b]);
  const result = new Array(values.length);
  for (let start = 0; start < values.length;) {
    let end = start + 1;
    while (end < values.length && values[order[end]] === values[order[start]]) end += 1;
    const rank = (start + end - 1) / 2;
    for (let i = start; i < end; i += 1) result[order[i]] = rank;
    start = end;
  }
  return result;
}

export function spearman(x, y) {
  return pearson(ranks(x), ranks(y));
}

export function pearson(x, y) {
  if (x.length !== y.length || x.length < 3) return null;
  const mean = values => values.reduce((sum,x) => sum+x,0) / values.length;
  const cx = mean(x), cy = mean(y);
  let product = 0, varA = 0, varB = 0;
  for (let i = 0; i < x.length; i += 1) {
    const u = x[i] - cx, v = y[i] - cy;
    product += u*v; varA += u*u; varB += v*v;
  }
  return varA && varB ? product / Math.sqrt(varA*varB) : null;
}

export function regionSummary(bundle, ion, mode) {
  const intensity = bundle.msi[mode][ion];
  const result = bundle.regions.map((name, id) => ({name, id, sum:0, n:0, detected:0}));
  for (let i = 0; i < intensity.length; i += 1) {
    if (mode === 'norm' && bundle.msi.tic[i] <= 0) continue;
    const region = result[bundle.msi.region[i]];
    if (!region || !region.id) continue;
    region.n += 1; region.sum += intensity[i]; region.detected += intensity[i] > 0 ? 1 : 0;
  }
  return result.filter(r => r.id && r.n).map(r => ({...r, mean:r.sum/r.n, detect:r.detected/r.n})).sort((a,b) => b.mean-a.mean);
}

function pairedTiles(bundle, region) {
  const g = !region && bundle.globalGrid ? bundle.globalGrid : bundle.grid;
  const idx = g.xy.map((_,i) => i).filter(i => (!region || g.region[i]===region) && g.stN[i]>0 && g.msiTicN[i]>0);
  const groups = new Map();
  for (const i of idx) {
    const key = g.tile_xy?.[i]?.join(',') ?? String(i);
    if (!groups.has(key)) groups.set(key, []);
    groups.get(key).push(i);
  }
  return {g, tiles:[...groups.values()], ionWeights:g.msiWeight ?? g.msiTicN};
}

const weightedMean = (indices, weights, read) =>
  indices.reduce((sum,i) => sum + weights[i]*read(i), 0) /
  indices.reduce((sum,i) => sum + weights[i], 0);

export function tilePairs(bundle, ion, region, gene) {
  const {g,tiles,ionWeights}=pairedTiles(bundle,region);
  return tiles.map(indices=>{
    const total = values => indices.reduce((sum,i)=>sum+values[i],0);
    const cp10k = weightedMean(indices,g.umi,i=>Math.expm1(g.genes[gene][i]));
    return {xy:g.tile_xy?.[indices[0]], msi:total(g.msiTicN),st:total(g.stN),umi:total(g.umi),
            ion:weightedMean(indices,ionWeights,i=>g.ions[ion][i]),cp10k,log:Math.log1p(cp10k)};
  });
}

export function associations(bundle, ion, region, minTiles=8, method='spearman') {
  const {g,tiles,ionWeights}=pairedTiles(bundle,region);
  const n = tiles.length;
  if (n < minTiles) return {n, rows:[]};
  const x = tiles.map(indices => weightedMean(indices,ionWeights,i => g.ions[ion][i]));
  const correlate = method === 'pearson' ? pearson : spearman;
  const rows = bundle.geneNames.flatMap((gene,j) => {
    const y = tiles.map(indices => Math.log1p(weightedMean(indices,g.umi,i => Math.expm1(g.genes[j][i]))));
    const rho = correlate(x,y);
    return rho === null ? [] : [{gene,index:j,rho}];
  }).sort((a,b) => Math.abs(b.rho)-Math.abs(a.rho) || a.gene.localeCompare(b.gene));
  return {n, rows};
}
