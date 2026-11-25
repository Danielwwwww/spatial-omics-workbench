const ntp = 'https://ntp.niehs.nih.gov';
const atlas = (path, title, image, caption, structures) => ({
  url: `${ntp}/atlas/us/${path}`,
  image: `${ntp}/sites/default/files/2023-08/${image}_large.jpg`,
  title, caption, structures
});

export const organAtlases = new Map([
  ['liver', atlas('hepatobiliary-system/liver', '肝脏',
    'usatlas-liver-figure-002_45473-25x_liver_semithin',
    'NTP 图 2 · 肝小叶甲苯胺蓝切片。1 中央静脉，2 肝血窦。门管区对照原图谱图 1。',
    [
      ['Central vein', '中央静脉', '肝板和肝血窦向静脉汇聚。'],
      ['Portal tract', '门管区', '寻找成簇的门静脉分支、小动脉和胆管。'],
      ['Bile duct', '胆管', '管腔由立方上皮围成，常伴随门管区血管。'],
      ['Hepatic parenchyma', '肝实质', '由肝细胞板及其间的肝血窦构成。'],
    ])],
  ['pancreas', atlas('alimentary-system/pancreas', '胰腺',
    'usatlas-pancreas-figure-001-pancreas-45841-18-3d-40x',
    'NTP 图 1 · 胰腺甲苯胺蓝切片。1 胰岛，2 外分泌腺泡，3 小叶间隔，4 小叶间导管。',
    [
      ['Islet of Langerhans', '胰岛', '识别与周围腺泡染色不同的内分泌细胞团。'],
      ['Exocrine acini', '外分泌腺泡', '细胞围绕细小腔隙呈腺泡排列。'],
      ['Interlobular duct', '小叶间导管', '在小叶之间寻找由上皮围成的管腔。'],
      ['Interlobular septum', '小叶间隔', '以分隔外分泌小叶的结缔组织为界。'],
    ])],
  ['kidney', atlas('urinary-system/kidney', '肾脏',
    'usatlas-kidney-figure-002-transverse-kidney-he',
    'NTP 图 2 · 啮齿类肾脏 H&E 横切。OSOM 为外髓外带，ISOM 为外髓内带。',
    [
      ['Renal cortex', '肾皮质', '从含肾小球与曲小管的外侧区域辨认。'],
      ['Outer medulla', '肾外髓', '从皮髓交界向内，参考外带和内带。'],
      ['Inner medulla', '肾内髓', '位于更内侧，含集合管与亨利袢细段。'],
      ['Renal papilla', '肾乳头', '在切面覆盖乳头时标记向肾盂突出的内髓末端。'],
      ['Glomerulus', '肾小球', '在皮质内识别位于肾小囊中的毛细血管团。'],
    ])],
  ['lung', atlas('respiratory-system/lung', '肺',
    'usatlas-lung-figure-002-44960-16-46b1-terminal-bronchiole-instilled-25x',
    'NTP 图 2 · 肺半薄切片。1 终末细支气管，2 肺泡管，3 肺泡，4 血管。',
    [
      ['Terminal bronchiole', '终末细支气管', '连续上皮围成气道，管壁缺少软骨。'],
      ['Alveolar duct', '肺泡管', '管腔沿途连接多个肺泡开口。'],
      ['Pulmonary alveoli', '肺泡区', '识别由薄肺泡隔分开的气腔。'],
      ['Pulmonary vessel', '肺血管', '结合内皮、管壁和腔内血细胞与气道区分。'],
    ])],
  ['small-intestine', atlas('alimentary-system/small-intestine', '小肠',
    'usatlas-sm_intestine-figure-001-45543-16-12-duodenum-15x',
    'NTP 图 1 · 十二指肠甲苯胺蓝切片。2 绒毛，3 十二指肠腺，6 黏膜下层，7–8 肌层。',
    [
      ['Villus', '肠绒毛', '从伸向肠腔的指状突起辨认。'],
      ['Intestinal crypt', '肠隐窝', '在相邻绒毛基底之间辨认向内延伸的腺体。'],
      ['Submucosa', '黏膜下层', '标记黏膜肌层与外肌层之间的区域。'],
      ['Muscularis externa', '外肌层', '通常可见内环行与外纵行两层平滑肌。'],
    ])],
]);
