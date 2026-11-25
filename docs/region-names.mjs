const atlasNames = {
  '组织外':'Outside tissue', '待复核区域':'Unassigned tissue',
  '初级运动皮层':'Primary motor area', '次级运动皮层':'Secondary motor area',
  '初级躯体感觉皮层':'Primary somatosensory area', '补充躯体感觉皮层':'Supplemental somatosensory area',
  '听觉皮层':'Auditory areas', '视觉皮层':'Visual areas', '压后皮层':'Retrosplenial area',
  '后顶叶联络皮层':'Posterior parietal association areas', '其他新皮层':'Other isocortex',
  '嗅觉区':'Olfactory areas', '嗅皮层':'Olfactory areas',
  '海马 CA1':'Hippocampal CA1', '海马 CA2':'Hippocampal CA2', '海马 CA3':'Hippocampal CA3',
  '齿状回':'Dentate gyrus', '海马':'Hippocampus', '下托':'Subiculum',
  '海马旁区':'Parahippocampal region', '其他海马结构':'Other hippocampal formation',
  '丘脑':'Thalamus', '下丘脑':'Hypothalamus', '中脑':'Midbrain',
  '纹状体':'Striatum', '苍白球区':'Pallidum', '杏仁核复合区':'Amygdalar complex',
  '皮层下板':'Cortical subplate', '小脑':'Cerebellum', '脑桥':'Pons', '延髓':'Medulla',
  '胼胝体':'Corpus callosum', '外囊':'External capsule', '内囊':'Internal capsule',
  '海马伞':'Fimbria', '其他纤维束':'Other fiber tracts', '脑室系统':'Ventricular system',
};
const chineseNames = {'Outside H&E':'组织外','Other tissue':'其他组织','Tissue':'组织'};
const compactNames = {
  'Primary somatosensory area':'Primary somatosensory',
  'Supplemental somatosensory area':'Supplemental somatosensory',
  'Posterior parietal association areas':'Posterior parietal',
  'Other hippocampal formation':'Hippocampal formation',
};

export function regionDisplay(name, manual=false) {
  const chinese = chineseNames[name] ?? (/[\u3400-\u9fff]/u.test(name) ? name : '');
  if (manual) return {english:name, chinese:'' , map:name};
  const side = name.match(/^(?:图谱|图像)?([左右])(?:侧)?/);
  const key = side ? name.slice(side[0].length).replace(/（粗标注）$/, '') : name;
  const translation = atlasNames[key];
  const position = translation && side ? (side[1]==='左'?'Left ':'Right ') : '';
  const english = position+(translation ?? name);
  return {english,chinese:chinese===english?'':chinese,map:position+(compactNames[translation] ?? translation ?? name)};
}
