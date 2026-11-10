// Viewer entry point: toolbar state, file loading (picker, drag and drop, ?a= / ?b= URL
// parameters), error display and the synchronization of the two views (design D8, D9).
import { View } from './view.js';
import { parseMultiPatch } from './gismo.js';

const $ = (sel) => document.querySelector(sel);
const panes = { a: $('#pane-a'), b: $('#pane-b') };
const views = { a: new View(panes.a), b: null };

function settings() {
  return {
    subdiv: Number($('#subdiv').value),
    elementLines: $('#lines-elements').checked,
    allLines: $('#lines-all').checked,
    colorMode: $('#color-mode').value,
    markers: $('#markers').checked,
    interfaces: $('#interfaces').checked,
    clipAxis: $('#clip-axis').value,
    clipPos: Number($('#clip-pos').value) / 1000,
    clipFlip: $('#clip-flip').checked,
    netMode: $('#net-mode').value,
  };
}

function applyAll() {
  const s = settings();
  $('#subdiv-out').textContent = s.subdiv;
  $('#clip-pos').disabled = $('#clip-flip').disabled = s.clipAxis === 'none';
  for (const v of Object.values(views)) if (v) v.apply(s);
}

for (const id of ['subdiv', 'lines-elements', 'lines-all', 'color-mode', 'markers', 'interfaces', 'clip-axis', 'clip-pos', 'clip-flip', 'net-mode'])
  $('#' + id).addEventListener('input', applyAll);

// ---- errors --------------------------------------------------------------------------------

let errorTimer = null;
function showError(msg) {
  const el = $('#error');
  el.textContent = msg;
  el.hidden = false;
  clearTimeout(errorTimer);
  errorTimer = setTimeout(() => { el.hidden = true; }, 15000);
}
$('#error').addEventListener('click', () => { $('#error').hidden = true; });

// ---- views and camera sync -----------------------------------------------------------------

let syncing = false;
function onCameraChange(src) {
  if (syncing || !$('#sync').checked) return;
  const dst = src === views.a ? views.b : views.a;
  if (!dst || !dst.mp) return;
  syncing = true;
  dst.copyCamera(src);
  syncing = false;
}
views.a.onCameraChange = onCameraChange;

function ensureView(key) {
  if (key === 'b' && !views.b) {
    panes.b.hidden = false;
    views.b = new View(panes.b);
    views.b.onCameraChange = onCameraChange;
    views.a.resize();
  }
  return views[key];
}

$('#clear-b').addEventListener('click', () => {
  if (!views.b) return;
  views.b.clear();
  views.b.renderer.dispose();
  views.b.host.innerHTML = '';
  views.b = null;
  panes.b.hidden = true;
  $('#file-b').value = '';
  views.a.resize();
});

$('#reset-view').addEventListener('click', () => {
  views.a.resetCamera();
  if (views.b && views.b.mp) views.b.copyCamera(views.a);
});

// ---- loading -------------------------------------------------------------------------------

// Parses `text` and shows it in view `key`; on error nothing is shown in that view.
function show(key, text, name, timings) {
  const view = ensureView(key);
  let t = performance.now();
  let mp;
  try {
    mp = parseMultiPatch(text);
  } catch (e) {
    view.clear();
    showError(`无法加载 ${name}：${e.message}`);
    return false;
  }
  timings.parse = performance.now() - t;
  view.setModel(mp, name, timings);
  view.apply(settings());
  const other = key === 'a' ? views.b : views.a;
  if (other && other.mp && $('#sync').checked) view.copyCamera(other);
  $('#error').hidden = true;
  window.dispatchEvent(new CustomEvent('viewer-loaded', { detail: { key, name } }));
  return true;
}

async function loadFile(key, file) {
  const t = performance.now();
  const text = await file.text();
  return show(key, text, file.name, { read: performance.now() - t });
}

async function loadUrl(key, url) {
  const t = performance.now();
  let text;
  try {
    const r = await fetch(url);
    if (!r.ok) throw new Error(`HTTP ${r.status}`);
    text = await r.text();
  } catch (e) {
    showError(`无法读取 ${url}：${e.message}`);
    return false;
  }
  return show(key, text, url, { read: performance.now() - t });
}

$('#file-a').addEventListener('change', (e) => e.target.files[0] && loadFile('a', e.target.files[0]));
$('#file-b').addEventListener('change', (e) => e.target.files[0] && loadFile('b', e.target.files[0]));

for (const [key, pane] of Object.entries(panes)) {
  pane.addEventListener('dragover', (e) => { e.preventDefault(); pane.classList.add('drag-over'); });
  pane.addEventListener('dragleave', () => pane.classList.remove('drag-over'));
  pane.addEventListener('drop', (e) => {
    e.preventDefault();
    pane.classList.remove('drag-over');
    const file = e.dataTransfer.files[0];
    if (file) loadFile(key, file);
  });
}
// Dropping a second file next to the first pane opens it on the right when B is closed.
document.body.addEventListener('dragover', (e) => e.preventDefault());
document.body.addEventListener('drop', (e) => e.preventDefault());

const params = new URLSearchParams(location.search);
(async () => {
  if (params.get('a')) await loadUrl('a', params.get('a'));
  if (params.get('b')) await loadUrl('b', params.get('b'));
  window.viewerReady = true;
})();

// Exposed for automated checks (headless browser) and the browser console.
window.viewer = { views, settings, applyAll, show, loadUrl };
