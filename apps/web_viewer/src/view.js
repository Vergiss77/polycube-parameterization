// One 3D view of a multipatch (design D4-D8): boundary mesh with element / subdivision lines,
// coloring, non-positive Jacobian markers, interfaces, clipping, picking and control net.
// Geometry is stored relative to the model center (Float32 precision near the origin).
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { boundarySides, interfaceSides, tessellateSides } from './tessellate.js';
import { jacobianQuality } from './quality.js';
import { controlNet } from './controlnet.js';
import { BAD, PLAIN, blockColoring, paletteColor, viridis, viridisGradient } from './colors.js';
import { elementCount } from './bspline.js';

const now = () => performance.now();

export class View {
  constructor(pane) {
    this.pane = pane;
    this.host = pane.querySelector('.canvas-host');
    this.infoEl = pane.querySelector('.info');
    this.legendEl = pane.querySelector('.legend');
    this.renderer = new THREE.WebGLRenderer({ antialias: true });
    this.renderer.setPixelRatio(window.devicePixelRatio);
    this.renderer.setClearColor(0x1d2026);
    this.renderer.localClippingEnabled = true;
    this.host.appendChild(this.renderer.domElement);

    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(35, 1, 0.01, 1000);
    this.camera.up.set(0, 0, 1);
    this.camera.position.set(3, -4, 2.5);
    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = false;
    this.scene.add(new THREE.HemisphereLight(0xffffff, 0x50505a, 1.5));
    const sun = new THREE.DirectionalLight(0xffffff, 1.6);
    sun.position.set(0.6, 0.8, 1.5);
    this.camera.add(sun);
    this.scene.add(this.camera);

    this.clipPlane = new THREE.Plane(new THREE.Vector3(0, 0, -1), 0);
    this.clipping = [];  // empty or [clipPlane], shared by all materials
    this.group = new THREE.Group();
    this.scene.add(this.group);

    this.onCameraChange = null;  // set by the app for camera synchronization
    this.controls.addEventListener('change', () => {
      this.render();
      if (this.onCameraChange) this.onCameraChange(this);
    });
    this._installPicking();
    new ResizeObserver(() => this.resize()).observe(this.host);
    this.resize();
  }

  // ---- model -------------------------------------------------------------------------------

  setModel(mp, name, timings) {
    this.clear();
    this.mp = mp;
    this.name = name;
    this.timings = { ...timings };
    let t = now();
    this.quality = jacobianQuality(mp);
    this.timings.quality = now() - t;
    this.colorIndex = blockColoring(mp);
    // Unique control points: shared control points have identical coordinates.
    const seen = new Set();
    let total = 0;
    for (const p of mp.patches)
      for (let i = 0; i < p.coefs.length; i += 3, ++total) seen.add(`${p.coefs[i]},${p.coefs[i + 1]},${p.coefs[i + 2]}`);
    this.controlPoints = { total, unique: seen.size };
    // Bounding box from the control points (convex hull property).
    const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
    for (const p of mp.patches)
      for (let i = 0; i < p.coefs.length; i += 3)
        for (let d = 0; d < 3; ++d) { lo[d] = Math.min(lo[d], p.coefs[i + d]); hi[d] = Math.max(hi[d], p.coefs[i + d]); }
    this.origin = lo.map((l, d) => 0.5 * (l + hi[d]));
    this.box = { lo: lo.map((l, d) => l - this.origin[d]), hi: hi.map((h, d) => h - this.origin[d]) };
    this.radius = 0.5 * Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]);
    this.selected = -1;
    this.pane.classList.add('loaded');
    this.resetCamera();
  }

  clear() {
    for (const o of [...this.group.children]) {
      this.group.remove(o);
      o.geometry?.dispose();
      o.material?.dispose();
    }
    this.layers = {};
    this.mp = null;
    this.infoEl.innerHTML = '';
    this.legendEl.hidden = true;
    this.pane.classList.remove('loaded');
    this.render();
  }

  // Applies all display settings; rebuilds only what changed.
  apply(settings) {
    if (!this.mp) return;
    const prev = this.settings ?? {};
    this.settings = { ...settings };
    if (prev.subdiv !== settings.subdiv || !this.layers.surface) this._buildSurfaces(settings.subdiv);
    if (prev.netMode !== settings.netMode || this.netFor !== this._netKey()) this._buildNet();
    if (!this.layers.markers) this._buildMarkers();
    this._applyColors();
    this.layers.elementLines.visible = settings.elementLines && !settings.allLines;
    this.layers.allLines.visible = settings.allLines;
    this.layers.markers.visible = settings.markers;
    this.layers.interfaces.visible = settings.interfaces;
    this._applyClip();
    this._updateInfo();
    this._updateLegend();
    this.render();
  }

  // ---- geometry layers ---------------------------------------------------------------------

  _positions32(pos64) {
    const p = new Float32Array(pos64.length);
    for (let i = 0; i < pos64.length; i += 3)
      for (let d = 0; d < 3; ++d) p[i + d] = pos64[i + d] - this.origin[d];
    return p;
  }

  _replace(key, object) {
    const old = this.layers[key];
    if (old) {
      this.group.remove(old);
      old.geometry?.dispose();
      old.material?.dispose();
    }
    this.layers[key] = object;
    if (object) this.group.add(object);
  }

  _buildSurfaces(subdiv) {
    let t = now();
    const B = tessellateSides(this.mp, boundarySides(this.mp), subdiv);
    const I = tessellateSides(this.mp, interfaceSides(this.mp), subdiv);
    this.timings.tessellate = now() - t;
    this.boundary = B;
    this.interfaceMesh = I;
    this.triangles = B.index.length / 3;

    const geo = new THREE.BufferGeometry();
    const pos = new THREE.BufferAttribute(this._positions32(B.positions), 3);
    geo.setAttribute('position', pos);
    geo.setAttribute('normal', new THREE.BufferAttribute(B.normals, 3));
    geo.setAttribute('color', new THREE.BufferAttribute(new Float32Array(B.positions.length), 3));
    geo.setIndex(new THREE.BufferAttribute(B.index, 1));
    geo.computeBoundingSphere();
    const mat = new THREE.MeshLambertMaterial({
      vertexColors: true, side: THREE.DoubleSide, clippingPlanes: this.clipping,
      polygonOffset: true, polygonOffsetFactor: 1, polygonOffsetUnits: 1,
    });
    this._replace('surface', new THREE.Mesh(geo, mat));

    const lines = (index, color, opacity) => {
      const g = new THREE.BufferGeometry();
      g.setAttribute('position', pos);
      g.setIndex(new THREE.BufferAttribute(index, 1));
      return new THREE.LineSegments(g, new THREE.LineBasicMaterial({
        color, transparent: opacity < 1, opacity, clippingPlanes: this.clipping }));
    };
    this._replace('elementLines', lines(B.elementLines, 0x15181d, 0.9));
    this._replace('allLines', lines(B.allLines, 0x15181d, 0.55));

    const ig = new THREE.BufferGeometry();
    ig.setAttribute('position', new THREE.BufferAttribute(this._positions32(I.positions), 3));
    ig.setAttribute('normal', new THREE.BufferAttribute(I.normals, 3));
    ig.setAttribute('color', new THREE.BufferAttribute(new Float32Array(I.positions.length), 3));
    ig.setIndex(new THREE.BufferAttribute(I.index, 1));
    this._replace('interfaces', new THREE.Mesh(ig, new THREE.MeshLambertMaterial({
      vertexColors: true, side: THREE.DoubleSide, transparent: true, opacity: 0.35, depthWrite: false,
      clippingPlanes: this.clipping })));
    this._buildHighlight();
  }

  _buildMarkers() {
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.BufferAttribute(this._positions32(this.quality.bad), 3));
    this._replace('markers', new THREE.Points(geo, new THREE.PointsMaterial({
      color: new THREE.Color(...BAD), size: 7, sizeAttenuation: false, depthTest: false, clippingPlanes: this.clipping })));
    this.layers.markers.renderOrder = 10;
  }

  _netKey() {
    const mode = this.settings?.netMode ?? 'none';
    return mode === 'selected' ? `selected:${this.selected}` : mode;
  }

  _buildNet() {
    this.netFor = this._netKey();
    const mode = this.settings.netMode;
    if (mode === 'none' || (mode === 'selected' && this.selected < 0)) {
      this._replace('net', null);
      this._replace('netPoints', null);
      this.netInfo = null;
      return;
    }
    const net = controlNet(this.mp, mode === 'selected' ? [this.selected] : undefined);
    this.netInfo = net;
    const pos = new THREE.BufferAttribute(this._positions32(net.positions), 3);
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', pos);
    g.setIndex(new THREE.BufferAttribute(net.index, 1));
    this._replace('net', new THREE.LineSegments(g, new THREE.LineBasicMaterial({
      color: 0xffc94d, transparent: true, opacity: 0.8, depthTest: false, clippingPlanes: this.clipping })));
    const pg = new THREE.BufferGeometry();
    pg.setAttribute('position', pos);
    this._replace('netPoints', new THREE.Points(pg, new THREE.PointsMaterial({
      color: 0xffe08a, size: 4, sizeAttenuation: false, depthTest: false, clippingPlanes: this.clipping })));
    this.layers.net.renderOrder = 8;
    this.layers.netPoints.renderOrder = 9;
  }

  _buildHighlight() {
    if (this.selected < 0) {
      this._replace('highlight', null);
      return;
    }
    const sides = this.mp.boundary.filter(([b]) => b === this.selected);
    const H = tessellateSides(this.mp, sides, this.settings?.subdiv ?? 8);
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(this._positions32(H.positions), 3));
    g.setIndex(new THREE.BufferAttribute(H.index, 1));
    this._replace('highlight', new THREE.Mesh(g, new THREE.MeshBasicMaterial({
      color: 0xffffff, transparent: true, opacity: 0.35, depthWrite: false, side: THREE.DoubleSide,
      clippingPlanes: this.clipping, polygonOffset: true, polygonOffsetFactor: -1, polygonOffsetUnits: -1 })));
  }

  // ---- coloring ----------------------------------------------------------------------------

  jacobianRange() {
    const s = this.settings;
    const lo = Number.isFinite(s?.sjMin) ? s.sjMin : Math.max(0, Math.min(this.quality.min, 0.999));
    return [lo, 1];
  }

  _colorFor(patch, sj, mode, range) {
    if (mode === 'block') return paletteColor(this.colorIndex[patch]);
    if (mode === 'jacobian') {
      if (sj <= 0) return BAD;
      return viridis((sj - range[0]) / Math.max(1e-12, range[1] - range[0]));
    }
    return PLAIN;
  }

  _applyColors() {
    const mode = this.settings.colorMode, range = this.jacobianRange();
    for (const [key, data] of [['surface', this.boundary], ['interfaces', this.interfaceMesh]]) {
      const attr = this.layers[key].geometry.getAttribute('color');
      for (let v = 0; v < data.patch.length; ++v) attr.array.set(this._colorFor(data.patch[v], data.scaledJacobian[v], mode, range), 3 * v);
      attr.needsUpdate = true;
    }
  }

  // ---- clipping ----------------------------------------------------------------------------

  // axis: 'none' | 0 | 1 | 2; pos in [0, 1] of the model box; flip keeps the upper side.
  _applyClip() {
    const { clipAxis, clipPos, clipFlip } = this.settings;
    this.clipping.length = 0;
    if (clipAxis === 'none' || clipAxis === undefined) return;
    const a = Number(clipAxis);
    const c = this.box.lo[a] + (this.box.hi[a] - this.box.lo[a]) * clipPos;
    const n = new THREE.Vector3();
    n.setComponent(a, clipFlip ? 1 : -1);
    this.clipPlane.set(n, clipFlip ? -c : c);
    this.clipping.push(this.clipPlane);
  }

  // ---- picking -----------------------------------------------------------------------------

  _installPicking() {
    const ray = new THREE.Raycaster();
    let down = null;
    const el = this.renderer.domElement;
    el.addEventListener('pointerdown', (e) => { down = [e.clientX, e.clientY]; });
    el.addEventListener('pointerup', (e) => {
      if (!down || Math.hypot(e.clientX - down[0], e.clientY - down[1]) > 4 || !this.mp) return;
      const r = el.getBoundingClientRect();
      const ndc = new THREE.Vector2(((e.clientX - r.left) / r.width) * 2 - 1, -((e.clientY - r.top) / r.height) * 2 + 1);
      this.pickAt(ndc);
    });
  }

  // Picks the block under normalized device coordinates (also used by tests); returns its id.
  pickAt(ndc) {
    const ray = new THREE.Raycaster();
    ray.setFromCamera(ndc, this.camera);
    const hits = ray.intersectObject(this.layers.surface).filter((h) =>
      this.clipping.every((p) => p.distanceToPoint(h.point) >= 0));
    this.select(hits.length ? this.boundary.patch[hits[0].face.a] : -1);
    return this.selected;
  }

  select(block) {
    this.selected = block;
    this._buildHighlight();
    if (this.settings) {
      if (this.settings.netMode === 'selected') this._buildNet();
      this._applyClip();
      this._updateInfo();
    }
    this.render();
  }

  // ---- panels ------------------------------------------------------------------------------

  _updateInfo() {
    const q = this.quality, mp = this.mp, f = (x, d = 3) => Number(x).toFixed(d);
    const deg = [...new Set(mp.patches.flatMap((p) => p.degree))].join(', ');
    const ms = (x) => (x === undefined ? '-' : `${Math.round(x)} ms`);
    let html = `<h3>${escapeHtml(this.name)}</h3><table>
      <tr><td>块 / 接口 / 边界面</td><td>${mp.patches.length} / ${mp.interfaces.length} / ${mp.boundary.length}</td></tr>
      <tr><td>控制点（去重 / 合计）</td><td>${this.controlPoints.unique} / ${this.controlPoints.total}</td></tr>
      <tr><td>次数</td><td>${deg}</td></tr>
      <tr><td>缩放 Jacobian 最小 / 平均</td><td class="${q.min <= 0 ? 'bad' : ''}">${f(q.min)} / ${f(q.mean)}</td></tr>
      <tr><td>非正点 / 评估点</td><td class="${q.nonpositive ? 'bad' : ''}">${q.nonpositive} / ${q.evalPoints}（${q.blocksWithNonpositive} 块）</td></tr>
      <tr><td>三角形</td><td>${this.triangles}</td></tr>
      <tr><td>耗时</td><td>读取 ${ms(this.timings.read)}，解析 ${ms(this.timings.parse)}，Jacobian ${ms(this.timings.quality)}，离散 ${ms(this.timings.tessellate)}</td></tr>
    </table>`;
    if (this.selected >= 0) {
      const p = mp.patches[this.selected], b = q.perBlock[this.selected];
      const el = p.knots.map(elementCount);
      html += `<div class="section"><table>
        <tr><td>选中块</td><td>${this.selected}</td></tr>
        <tr><td>次数 / 单元数</td><td>${p.degree.join(' × ')} / ${el.join(' × ')}</td></tr>
        <tr><td>控制点</td><td>${p.n.join(' × ')} = ${p.n[0] * p.n[1] * p.n[2]}</td></tr>
        <tr><td>缩放 Jacobian 最小</td><td class="${b.min <= 0 ? 'bad' : ''}">${f(b.min)}</td></tr>
        <tr><td>非正点</td><td class="${b.nonpositive ? 'bad' : ''}">${b.nonpositive}</td></tr>
      </table></div>`;
    }
    if (this.netInfo) html += `<div class="section">控制网：${this.netInfo.pointCount} 个点，${this.netInfo.segmentCount} 条线</div>`;
    this.infoEl.innerHTML = html;
  }

  _updateLegend() {
    if (this.settings.colorMode !== 'jacobian') {
      this.legendEl.hidden = true;
      return;
    }
    const [lo, hi] = this.jacobianRange();
    this.legendEl.hidden = false;
    this.legendEl.innerHTML = `缩放 Jacobian<div class="bar" style="background:${viridisGradient()}"></div>
      <div class="ticks"><span>${lo.toFixed(2)}</span><span>${((lo + hi) / 2).toFixed(2)}</span><span>${hi.toFixed(2)}</span></div>
      <div><span style="color:rgb(${BAD.map((c) => Math.round(255 * c)).join(',')})">■</span> ≤ 0</div>`;
  }

  // ---- camera ------------------------------------------------------------------------------

  resetCamera() {
    const r = this.radius || 1;
    this.camera.near = r / 100;
    this.camera.far = r * 100;
    this.camera.position.set(1.6 * r, -2.2 * r, 1.4 * r);
    this.controls.target.set(0, 0, 0);
    this.camera.updateProjectionMatrix();
    this.controls.update();
  }

  copyCamera(other) {
    this.camera.position.copy(other.camera.position);
    this.camera.quaternion.copy(other.camera.quaternion);
    this.camera.near = other.camera.near;
    this.camera.far = other.camera.far;
    this.camera.updateProjectionMatrix();
    this.controls.target.copy(other.controls.target);
    this.controls.update();
  }

  resize() {
    const w = Math.max(1, this.host.clientWidth), h = Math.max(1, this.host.clientHeight);
    this.renderer.setSize(w, h, false);
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
    this.render();
  }

  render() {
    this.renderer.render(this.scene, this.camera);
  }
}

function escapeHtml(s) {
  return String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' })[c]);
}
