const runtimeURL = 'https://cdn.jsdelivr.net/pyodide/v0.29.4/full/';
const h5adWheels = [
  'https://files.pythonhosted.org/packages/e0/b1/0542e0cab6f49f151a2d7a42400f84f706fc0b64e85dc1f56708b2e9fd37/array_api_compat-1.12.0-py3-none-any.whl',
  'https://files.pythonhosted.org/packages/ef/82/7a9d0550484a62c6da82858ee9419f3dd1ccc9aa1c26a1e43da3ecd20b0d/natsort-8.4.0-py3-none-any.whl',
  'https://files.pythonhosted.org/packages/0a/4b/ab615fea52e34579d5c6c7dba86b4f9d7f3cdb6a170b348ec49f34cf4355/anndata-0.11.4-py3-none-any.whl',
];
const modules = ['__init__', 'pipeline', 'registration', 'manual_regions', 'st_display'];
const mount = '/input';
let python, h5adReady = false;

function progress(message) {
  self.postMessage({message});
}

async function prepare() {
  if (python) return;
  progress('正在加载浏览器计算组件，首次使用需要联网。');
  importScripts(runtimeURL + 'pyodide.js');
  const runtime = await loadPyodide({indexURL: runtimeURL});
  await runtime.loadPackage(['numpy', 'scipy', 'pandas', 'pillow']);
  runtime.FS.mkdir('/workbench');
  runtime.FS.mkdir('/workbench/compute');
  runtime.FS.mkdir(mount);
  await Promise.all(modules.map(async name => {
    const file = name + '.py';
    const response = await fetch(new URL(file, self.location.href));
    if (!response.ok) throw new Error('配准程序加载失败');
    runtime.FS.writeFile('/workbench/compute/' + file, new Uint8Array(await response.arrayBuffer()));
  }));
  runtime.runPython(`import sys, json
from pathlib import Path
sys.path.insert(0, "/workbench")
from compute.pipeline import analyze`);
  python = runtime;
}

async function prepareH5AD() {
  if (h5adReady) return;
  progress('正在加载 H5AD 读取组件。');
  await python.loadPackage(['h5py', 'packaging', 'micropip']);
  await python.loadPackage(h5adWheels);
  python.runPython('import anndata');
  h5adReady = true;
}

self.onmessage = async ({data: {files, config}}) => {
  try {
    await prepare();
    if (files.st.name.toLowerCase().endsWith('.h5ad')) await prepareH5AD();
    const paths = {};
    const blobs = Object.entries(files).map(([key, file]) => {
      const suffix = file.name.toLowerCase().match(/(?:\.[a-z0-9]+)?\.gz$|\.[a-z0-9]+$/)?.[0] || '';
      const name = key + suffix;
      paths[key] = mount + '/' + name;
      return {name, data: file};
    });
    python.FS.mount(python.FS.filesystems.WORKERFS, {blobs}, mount);
    try {
      progress('正在浏览器内读取数据、配准和计算空间统计。');
      python.globals.set('request_json', JSON.stringify({paths, config}));
      const result = python.runPython(`def run_request(request):
    paths = {key: Path(value) for key, value in request["paths"].items()}
    result = analyze(config=request["config"], **paths)
    return json.dumps(result, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
run_request(json.loads(request_json))`);
      self.postMessage({result: JSON.parse(result)});
    } finally {
      python.globals.delete('request_json');
      python.FS.unmount(mount);
    }
  } catch (error) {
    self.postMessage({error: error.message.split('\n').filter(Boolean).at(-1) || '浏览器计算失败。'});
  }
};
