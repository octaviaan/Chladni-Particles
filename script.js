import * as THREE from "three";
import { Pane } from "tweakpane";
import * as EssentialsPlugin from "tweakpane-plugin-essentials";

// --- Configuration ---
const CONFIG = {
  particleCount: 400000,
  gridSize: 256, // Resolution of the underlying energy field
  settleStrength: 1.5, // How fast particles move to nodal lines
  jitter: 0.1, // Brownian motion to keep them "alive"
  drag: 0.85, // Friction
  speedLimit: 2.0,
  viewScale: 600, // Size of the world
  color: "#72dcff",
  particleSize: 2,
  particleOpacity: 1.0,
  cameraControl: { x: 0, y: 0, z: 1.9 },
  exportBaseSize: 2048,
  exportOpaque: false,

  modeCount: 4,
  mRange: { min: 2, max: 12 },
  nRange: { min: 2, max: 12 },

  rectAspect: 1.9,

  waveTypeA: "Cartesian",
  waveTypeB: "Radial",
  waveMix: 0.5,

  // Less organic controls
  integerModes: true, // round m,n to integers for cleaner straight lines
  randomPhase: true,
  randomRotation: true,
  rotationMaxDeg: 30,

  // Wave types:
  // "Cartesian", "Radial", "Spiral", "Hexagonal", "Quasicrystal", "Lissajous",
  // "Moire", "Chevron", "Ring", "Flower", "GridWarp", "Diamond", "Square",
  // "Lattice", "Kaleidoscope", "Voronoi", "Perlin", "Superellipse", "Bessel",
  // "TriLattice", "Rose", "Astroid", "Checker"
};

// --- Global Variables ---
let scene, camera, renderer, geometry, points;
let positions, velocities; // Float32Arrays for high performance
let refreshUI = null;
let suppressUIEvents = false;
let fpsGraph = null;

// Field Data (Pre-computed grid)
let fieldSize = 0;
let energy;
let gradX;
let gradY;

// Initial Chladni Modes (m, n)
let modes = [];

// --- Initialization ---
init();
allocateField();
initModes();
rebuildField();
animate();

function init() {
  // 1. Scene & Camera
  scene = new THREE.Scene();
  const aspect = window.innerWidth / window.innerHeight;
  const frustumSize = CONFIG.viewScale * 2;
  camera = new THREE.OrthographicCamera(
    (-frustumSize * aspect) / 2,
    (frustumSize * aspect) / 2,
    frustumSize / 2,
    -frustumSize / 2,
    1,
    3000,
  );
  camera.position.set(0, 0, 500); // True top-down view
  applyCameraFromControl(CONFIG.cameraControl);

  // 2. Renderer
  renderer = new THREE.WebGLRenderer({
    antialias: true,
    powerPreference: "high-performance",
    preserveDrawingBuffer: true,
  });
  renderer.setSize(window.innerWidth, window.innerHeight);
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  renderer.setClearColor(0x000000, 1);
  document.body.appendChild(renderer.domElement);

  // 3. Particle System Setup
  buildParticles();

  // 5. Material
  const sprite = new THREE.TextureLoader().load(
    "https://threejs.org/examples/textures/sprites/disc.png",
  );

  const material = new THREE.PointsMaterial({
    color: CONFIG.color,
    vertexColors: false,
    size: CONFIG.particleSize,
    map: sprite,
    alphaTest: 0.0,
    transparent: true,
    opacity: CONFIG.particleOpacity,
    blending: THREE.AdditiveBlending,
    depthWrite: false,
    sizeAttenuation: false,
  });

  points = new THREE.Points(geometry, material);
  scene.add(points);

  // 6. Event Listeners
  window.addEventListener("resize", onWindowResize);

  setupGUI();
}

// --- Physics & Math ---

function allocateField() {
  fieldSize = CONFIG.gridSize * CONFIG.gridSize;
  energy = new Float32Array(fieldSize);
  gradX = new Float32Array(fieldSize);
  gradY = new Float32Array(fieldSize);
}

function buildParticles() {
  geometry = new THREE.BufferGeometry();

  positions = new Float32Array(CONFIG.particleCount * 3);
  velocities = new Float32Array(CONFIG.particleCount * 2);

  for (let i = 0; i < CONFIG.particleCount; i++) {
    const [x, y] = randomPointInShape();
    positions[i * 3 + 0] = x;
    positions[i * 3 + 1] = y;
    positions[i * 3 + 2] = 0;

    velocities[i * 2 + 0] = 0;
    velocities[i * 2 + 1] = 0;
  }

  geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
}

function rebuildParticles() {
  const oldGeometry = points.geometry;
  buildParticles();
  points.geometry = geometry;
  oldGeometry.dispose();
}

function rebuildField() {
  const G = CONFIG.gridSize;
  const typeA = CONFIG.waveTypeA;
  const typeB = CONFIG.waveTypeB;
  const mix = CONFIG.waveMix;

  const lerp = (a, b, t) => a + (b - a) * t;
  const smoothstep = (t) => t * t * (3 - 2 * t);
  const hash2 = (x, y) => {
    const s = Math.sin(x * 127.1 + y * 311.7) * 43758.5453;
    return s - Math.floor(s);
  };
  const valueNoise = (x, y) => {
    const x0 = Math.floor(x);
    const y0 = Math.floor(y);
    const x1 = x0 + 1;
    const y1 = y0 + 1;
    const sx = smoothstep(x - x0);
    const sy = smoothstep(y - y0);
    const n00 = hash2(x0, y0) * 2 - 1;
    const n10 = hash2(x1, y0) * 2 - 1;
    const n01 = hash2(x0, y1) * 2 - 1;
    const n11 = hash2(x1, y1) * 2 - 1;
    const ix0 = lerp(n00, n10, sx);
    const ix1 = lerp(n01, n11, sx);
    return lerp(ix0, ix1, sy);
  };

  const baseCartesian = (cx, cy, mode) => {
    const cosR = mode.cos;
    const sinR = mode.sin;

    // rotated coordinates
    let rx = cx * cosR - cy * sinR;
    let ry = cx * sinR + cy * cosR;

    // controlled "bend" that stays structured (not organic noise)
    const u = rx + 0.5;
    const v = ry + 0.5;

    return (
      Math.sin(mode.m * Math.PI * u + mode.px) *
      Math.sin(mode.n * Math.PI * v + mode.py)
    );
  };

  const waveValue = (type, cx, cy, idx, mode) => {
    if (type === "Radial") {
      const r = Math.sqrt(cx * cx + cy * cy) * 2.0;
      const theta = Math.atan2(cy, cx);
      const mInt = Math.round(mode.m);
      return (
        Math.sin(mode.n * Math.PI * r + mode.px) *
        Math.cos(mInt * theta + mode.py)
      );
    }
    if (type === "Spiral") {
      const r = Math.sqrt(cx * cx + cy * cy) * 2.0;
      const theta = Math.atan2(cy, cx);
      const mInt = Math.round(mode.m);
      return Math.sin(mode.n * Math.PI * r + mInt * theta + mode.px);
    }
    if (type === "HyperSpiral") {
      const r = Math.sqrt(cx * cx + cy * cy);
      const theta = Math.atan2(cy, cx);
      return Math.sin(
        mode.n * Math.log(r + 0.001) * 5 + mode.m * theta + mode.px,
      );
    }
    if (type === "Parabolic") {
      const u = cx * cx - cy * cy;
      const v = 2 * cx * cy;
      return Math.sin(mode.m * Math.PI * u) * Math.sin(mode.n * Math.PI * v);
    }
    if (type === "Hexagonal") {
      const cosR = mode.cos;
      const sinR = mode.sin;
      const rx = cx * cosR - cy * sinR;
      const ry = cx * sinR + cy * cosR;
      const kFreq = mode.m * Math.PI;
      const v1 = rx;
      const v2 = -0.5 * rx + 0.866 * ry;
      const v3 = -0.5 * rx - 0.866 * ry;
      return (
        Math.sin(kFreq * v1 + mode.px) +
        Math.sin(kFreq * v2 + mode.py) +
        Math.sin(kFreq * v3)
      );
    }
    if (type === "Quasicrystal") {
      const cosR = mode.cos;
      const sinR = mode.sin;
      const rx = cx * cosR - cy * sinR;
      const ry = cx * sinR + cy * cosR;
      const kFreq = mode.m * Math.PI;
      let sum = 0;
      for (let j = 0; j < 5; j++) {
        const theta = (Math.PI * 2 * j) / 5;
        const v = rx * Math.cos(theta) + ry * Math.sin(theta);
        sum += Math.cos(kFreq * v + mode.px);
      }
      return sum;
    }
    if (type === "Gyroid") {
      const scale = mode.m * 10;
      return (
        Math.sin(cx * scale) * Math.cos(cy * scale) +
        Math.sin(cy * scale) * Math.cos(mode.px)
      );
    }
    if (type === "Lissajous") {
      const cosR = mode.cos;
      const sinR = mode.sin;
      const rx = cx * cosR - cy * sinR;
      const ry = cx * sinR + cy * cosR;
      return (
        Math.sin(mode.m * Math.PI * rx + mode.px) *
        Math.sin(mode.n * Math.PI * ry + mode.py)
      );
    }
    if (type === "Moire") {
      const cosR = mode.cos;
      const sinR = mode.sin;
      const rx = cx * cosR - cy * sinR;
      const ry = cx * sinR + cy * cosR;
      const rot2 = mode.px * 0.25;
      const cos2 = Math.cos(rot2);
      const sin2 = Math.sin(rot2);
      const rx2 = rx * cos2 - ry * sin2;
      const ry2 = rx * sin2 + ry * cos2;
      return (
        Math.sin(mode.m * Math.PI * rx + mode.px) *
          Math.sin(mode.n * Math.PI * ry + mode.py) +
        0.8 *
          Math.sin(mode.m * 1.05 * Math.PI * rx2 + mode.px) *
          Math.sin(mode.n * 0.95 * Math.PI * ry2 + mode.py)
      );
    }
    if (type === "Chevron") {
      const cosR = mode.cos;
      const sinR = mode.sin;
      const rx = cx * cosR - cy * sinR;
      const ry = cx * sinR + cy * cosR;
      const d1 = rx + ry;
      const d2 = rx - ry;
      return (
        Math.sin(mode.m * Math.PI * d1 + mode.px) *
        Math.sin(mode.n * Math.PI * d2 + mode.py)
      );
    }
    if (type === "Ring") {
      const r = Math.sqrt(cx * cx + cy * cy) * 2.0;
      return Math.sin(mode.m * Math.PI * r + mode.px);
    }
    if (type === "Cassini") {
      const a = mode.m * 0.1;
      const d1 = (cx - a) * (cx - a) + cy * cy;
      const d2 = (cx + a) * (cx + a) + cy * cy;
      const val = Math.sqrt(d1 * d2);
      return Math.sin(val * 20 - mode.n * 5);
    }
    if (type === "Flower") {
      const r = Math.sqrt(cx * cx + cy * cy) * 2.0;
      const theta = Math.atan2(cy, cx);
      const nInt = Math.round(mode.n);
      return (
        Math.sin(mode.m * Math.PI * r + mode.px) *
        Math.cos(nInt * theta + mode.py)
      );
    }
    if (type === "Interference") {
      const d1 = Math.sqrt((cx + 0.3) ** 2 + cy ** 2);
      const d2 = Math.sqrt((cx - 0.3) ** 2 + cy ** 2);
      return Math.sin(mode.m * 20 * d1) + Math.sin(mode.m * 20 * d2 + mode.px);
    }
    if (type === "GridWarp") {
      const cosR = mode.cos;
      const sinR = mode.sin;
      const rx = cx * cosR - cy * sinR;
      const ry = cx * sinR + cy * cosR;
      const warp = 0.35;
      const u = rx + warp * Math.sin(mode.n * Math.PI * ry + mode.py);
      const v = ry + warp * Math.sin(mode.m * Math.PI * rx + mode.px);
      return (
        Math.sin(mode.m * Math.PI * u + mode.px) *
        Math.sin(mode.n * Math.PI * v + mode.py)
      );
    }
    if (type === "Diamond") {
      const r = (Math.abs(cx) + Math.abs(cy)) * 2.0;
      return Math.sin(mode.m * Math.PI * r + mode.px);
    }
    if (type === "Square") {
      const r = Math.max(Math.abs(cx), Math.abs(cy)) * 2.0;
      return Math.sin(mode.m * Math.PI * r + mode.px);
    }
    if (type === "Lattice") {
      const cosR = mode.cos;
      const sinR = mode.sin;
      const rx = cx * cosR - cy * sinR;
      const ry = cx * sinR + cy * cosR;
      return (
        Math.cos(mode.m * Math.PI * rx + mode.px) *
        Math.cos(mode.n * Math.PI * ry + mode.py)
      );
    }
    if (type === "Kaleidoscope") {
      const r = Math.sqrt(cx * cx + cy * cy) * 2.0;
      const baseTheta = Math.atan2(cy, cx) + mode.py;
      const slices = Math.max(3, Math.round(mode.n));
      const wedge = Math.PI / slices;
      let theta = (baseTheta + Math.PI) % (2 * wedge);
      if (theta > wedge) theta = 2 * wedge - theta;
      return (
        Math.sin(mode.m * Math.PI * r + mode.px) * Math.cos(theta * slices)
      );
    }
    if (type === "Voronoi") {
      let minD = 1e9;
      for (let j = 0; j < 6; j++) {
        const sx = hash2(j + mode.m * 1.3, mode.n * 2.1) - 0.5;
        const sy = hash2(j + mode.n * 1.7, mode.m * 2.5) - 0.5;
        const dx = cx - sx;
        const dy = cy - sy;
        const d = Math.sqrt(dx * dx + dy * dy);
        if (d < minD) minD = d;
      }
      return Math.sin(mode.m * Math.PI * minD * 2 + mode.px);
    }
    if (type === "Perlin") {
      const cosR = mode.cos;
      const sinR = mode.sin;
      const rx = cx * cosR - cy * sinR;
      const ry = cx * sinR + cy * cosR;
      const scale = Math.max(1.5, mode.m);
      const n = valueNoise(rx * scale + mode.px, ry * scale + mode.py);
      return n;
    }
    if (type === "Superellipse") {
      const p = 0.8 + (Math.abs(mode.m) % 4) * 0.4;
      const rx = Math.abs(cx) * 2.0;
      const ry = Math.abs(cy) * 2.0;
      const r = Math.pow(Math.pow(rx, p) + Math.pow(ry, p), 1 / p) * 0.5;
      return Math.sin(mode.n * Math.PI * r + mode.px);
    }
    if (type === "Bessel") {
      const r = Math.sqrt(cx * cx + cy * cy) * 2.0;
      const theta = Math.atan2(cy, cx);
      const nInt = Math.max(0, Math.round(mode.n));
      const x = Math.max(1e-4, mode.m * Math.PI * r + 1e-4);
      const j = Math.sin(x - (nInt * Math.PI) / 2) / Math.sqrt(x);
      return j * Math.cos(nInt * theta + mode.py);
    }
    if (type === "TriLattice") {
      const cosR = mode.cos;
      const sinR = mode.sin;
      const rx = cx * cosR - cy * sinR;
      const ry = cx * sinR + cy * cosR;
      const k = mode.m * Math.PI;
      const a0 = k * rx + mode.px;
      const a1 = k * (rx * -0.5 + ry * 0.866) + mode.py;
      const a2 = k * (rx * -0.5 - ry * 0.866);
      return Math.sin(a0) + Math.sin(a1) + Math.sin(a2);
    }
    if (type === "Astroid") {
      const rx = Math.abs(cx) * 2.0;
      const ry = Math.abs(cy) * 2.0;
      const r = Math.pow(Math.pow(rx, 2 / 3) + Math.pow(ry, 2 / 3), 3 / 2);
      return Math.sin(mode.m * Math.PI * r + mode.px);
    }
    if (type === "Checker") {
      const cosR = mode.cos;
      const sinR = mode.sin;
      const rx = cx * cosR - cy * sinR;
      const ry = cx * sinR + cy * cosR;
      const a = mode.m * Math.PI;
      const b = mode.n * Math.PI;
      return (
        Math.sin(a * rx + mode.px) * Math.sin(b * ry + mode.py) +
        Math.sin(b * rx + mode.px) * Math.sin(a * ry + mode.py)
      );
    }
    // Cartesian (default)
    return baseCartesian(cx, cy, mode);
  };

  // 1) Compute Energy Field
  const gamma = 1.0;
  for (let y = 0; y < G; y++) {
    for (let x = 0; x < G; x++) {
      const tx = x / (G - 1);
      const ty = y / (G - 1);
      const cx = tx - 0.5;
      const cy = ty - 0.5;

      const idx = y * G + x;
      let phi = 0;
      for (let k = 0; k < modes.length; k++) {
        const mode = modes[k];
        const waveA = waveValue(typeA, cx, cy, idx, mode);
        const waveB = waveValue(typeB, cx, cy, idx, mode);
        const wave = Math.tanh(waveA * (1 - mix) + waveB * mix);
        phi += mode.a * wave;
      }
      let e = phi * phi;
      if (gamma !== 1.0) e = Math.pow(e, gamma);
      energy[idx] = e;
    }
  }

  // 2) Compute Gradients
  const range = CONFIG.viewScale * 2;
  const cellSize = range / (G - 1);

  for (let y = 0; y < G; y++) {
    const y0 = Math.max(0, y - 1);
    const y1 = Math.min(G - 1, y + 1);

    for (let x = 0; x < G; x++) {
      const x0 = Math.max(0, x - 1);
      const x1 = Math.min(G - 1, x + 1);

      const idx = y * G + x;

      const eL = energy[y * G + x0];
      const eR = energy[y * G + x1];
      const eU = energy[y0 * G + x];
      const eD = energy[y1 * G + x];

      gradX[idx] = (eR - eL) / (2 * cellSize);
      gradY[idx] = (eD - eU) / (2 * cellSize);
    }
  }
}

function updateParticles() {
  const G = CONFIG.gridSize;
  const range = CONFIG.viewScale;
  const fullRange = range * 2;
  const pos = positions;
  const vel = velocities;
  const count = CONFIG.particleCount;

  const gridScale = (G - 1) / fullRange;
  const rectRy = Math.max(1, Math.round(range / CONFIG.rectAspect));

  for (let i = 0; i < count; i++) {
    const i3 = i * 3;
    const i2 = i * 2;

    let x = pos[i3];
    let y = pos[i3 + 1];
    let vx = vel[i2];
    let vy = vel[i2 + 1];

    // 1. Boundary Check
    let inside = true;
    if (x < -range || x > range || y < -rectRy || y > rectRy) inside = false;

    if (!inside) {
      const respawn = randomPointInShape();
      x = respawn[0];
      y = respawn[1];
      vx = 0;
      vy = 0;
    }

    // 2. Bilinear Interpolation of Gradients
    let gx_pos = (x + range) * gridScale;
    let gy_pos = (y + range) * gridScale;

    gx_pos = Math.max(0, Math.min(G - 1.001, gx_pos));
    gy_pos = Math.max(0, Math.min(G - 1.001, gy_pos));

    const g_x0 = Math.floor(gx_pos);
    const g_y0 = Math.floor(gy_pos);
    const tx = gx_pos - g_x0;
    const ty = gy_pos - g_y0;

    const idx00 = g_y0 * G + g_x0;
    const idx10 = idx00 + 1;
    const idx01 = idx00 + G;
    const idx11 = idx01 + 1;

    const gxVal =
      (gradX[idx00] * (1 - tx) + gradX[idx10] * tx) * (1 - ty) +
      (gradX[idx01] * (1 - tx) + gradX[idx11] * tx) * ty;

    const gyVal =
      (gradY[idx00] * (1 - tx) + gradY[idx10] * tx) * (1 - ty) +
      (gradY[idx01] * (1 - tx) + gradY[idx11] * tx) * ty;

    const eVal =
      (energy[idx00] * (1 - tx) + energy[idx10] * tx) * (1 - ty) +
      (energy[idx01] * (1 - tx) + energy[idx11] * tx) * ty;

    // 3. Update Velocity (Gradient Descent)
    vx -= gxVal * CONFIG.settleStrength;
    vy -= gyVal * CONFIG.settleStrength;

    vx += (Math.random() - 0.5) * CONFIG.jitter;
    vy += (Math.random() - 0.5) * CONFIG.jitter;

    vx *= CONFIG.drag;
    vy *= CONFIG.drag;

    const speed = Math.sqrt(vx * vx + vy * vy);
    if (speed > CONFIG.speedLimit) {
      const scale = CONFIG.speedLimit / speed;
      vx *= scale;
      vy *= scale;
    }

    // 4. Update Position
    x += vx;
    y += vy;

    pos[i3] = x;
    pos[i3 + 1] = y;
    pos[i3 + 2] = eVal * 2;

    vel[i2] = vx;
    vel[i2 + 1] = vy;
  }

  geometry.attributes.position.needsUpdate = true;
}

function onWindowResize() {
  const aspect = window.innerWidth / window.innerHeight;
  const frustumSize = CONFIG.viewScale * 2;
  camera.left = (-frustumSize * aspect) / 2;
  camera.right = (frustumSize * aspect) / 2;
  camera.top = frustumSize / 2;
  camera.bottom = -frustumSize / 2;
  camera.updateProjectionMatrix();
  renderer.setSize(window.innerWidth, window.innerHeight);
}

function randomPointInShape() {
  const range = CONFIG.viewScale;
  const rectRy = Math.max(1, Math.round(range / CONFIG.rectAspect));

  return [
    (Math.random() - 0.5) * 2 * range,
    (Math.random() - 0.5) * 2 * rectRy,
  ];
}

function quantizeModeValue(v) {
  if (!CONFIG.integerModes) return v;
  return Math.max(1, Math.round(v));
}

function modeValueAt(index, count, min, max, flip = false) {
  if (count <= 1) return (min + max) * 0.5;
  const t = index / (count - 1);
  const u = flip ? 1 - t : t;
  return min + (max - min) * u;
}

function randomizeModes() {
  const deviceMemory = navigator.deviceMemory || 8;
  const cores = navigator.hardwareConcurrency || 8;
  const perfScale = Math.min(
    1,
    Math.max(0.5, (deviceMemory / 8) * (cores / 8)),
  );
  const maxModes = Math.max(4, Math.round(10 * perfScale));
  CONFIG.modeCount = Math.floor(1 + Math.random() * maxModes);
  CONFIG.mRange.min = Math.random() * 20;
  CONFIG.mRange.max =
    CONFIG.mRange.min + Math.random() * (40 - CONFIG.mRange.min);
  CONFIG.nRange.min = Math.random() * 20;
  CONFIG.nRange.max =
    CONFIG.nRange.min + Math.random() * (40 - CONFIG.nRange.min);

  CONFIG.integerModes = Math.random() < 0.5;
  CONFIG.randomPhase = Math.random() < 0.5;
  CONFIG.randomRotation = Math.random() < 0.5;
  CONFIG.rotationMaxDeg = Math.floor(Math.random() * 91);

  modes = [];
  for (let i = 0; i < CONFIG.modeCount; i++) {
    const rot =
      CONFIG.randomRotation && CONFIG.rotationMaxDeg > 0
        ? ((Math.random() * 2 - 1) * CONFIG.rotationMaxDeg * Math.PI) / 180
        : 0;

    const mRaw = modeValueAt(
      i,
      CONFIG.modeCount,
      CONFIG.mRange.min,
      CONFIG.mRange.max,
    );
    const nRaw = modeValueAt(
      i,
      CONFIG.modeCount,
      CONFIG.nRange.min,
      CONFIG.nRange.max,
      true,
    );

    const m = quantizeModeValue(mRaw);
    const n = quantizeModeValue(nRaw);

    modes.push({
      m,
      n,
      a: 1 - i * 0.2,
      px: CONFIG.randomPhase ? Math.random() * Math.PI * 2 : 0,
      py: CONFIG.randomPhase ? Math.random() * Math.PI * 2 : 0,
      cos: Math.cos(rot),
      sin: Math.sin(rot),
    });
  }

  normalizeModeAmplitudes();
  rebuildField();

  if (refreshUI) {
    suppressUIEvents = true;
    refreshUI();
    suppressUIEvents = false;
  }
}

function initModes() {
  modes = [];
  for (let i = 0; i < CONFIG.modeCount; i++) {
    const rot =
      CONFIG.randomRotation && CONFIG.rotationMaxDeg > 0
        ? ((Math.random() * 2 - 1) * CONFIG.rotationMaxDeg * Math.PI) / 180
        : 0;

    const mRaw = modeValueAt(
      i,
      CONFIG.modeCount,
      CONFIG.mRange.min,
      CONFIG.mRange.max,
    );
    const nRaw = modeValueAt(
      i,
      CONFIG.modeCount,
      CONFIG.nRange.min,
      CONFIG.nRange.max,
      true,
    );

    const m = quantizeModeValue(mRaw);
    const n = quantizeModeValue(nRaw);

    modes.push({
      m,
      n,
      a: 1 - i * 0.2,
      px: CONFIG.randomPhase ? Math.random() * Math.PI * 2 : 0,
      py: CONFIG.randomPhase ? Math.random() * Math.PI * 2 : 0,
      cos: Math.cos(rot),
      sin: Math.sin(rot),
    });
  }
  normalizeModeAmplitudes();
}

function rebuildModesFromConfig() {
  initModes();
  rebuildField();
}

function normalizeModeAmplitudes() {
  let sum = 0;
  for (let i = 0; i < modes.length; i++) sum += Math.abs(modes[i].a);
  if (sum <= 0) return;
  const inv = 1 / sum;
  for (let i = 0; i < modes.length; i++) modes[i].a *= inv;
}

function animate() {
  requestAnimationFrame(animate);
  if (fpsGraph) fpsGraph.begin();
  updateParticles();
  renderer.render(scene, camera);
  if (fpsGraph) fpsGraph.end();
}

function applyCameraZoom(value) {
  camera.zoom = value;
  camera.updateProjectionMatrix();
}

function applyCameraPan(x, y) {
  camera.position.x = x;
  camera.position.y = y;
  camera.lookAt(x, y, 0);
  camera.updateProjectionMatrix();
}

function applyCameraFromControl(control) {
  applyCameraZoom(control.z);
  applyCameraPan(control.x, control.y);
}

function setupGUI() {
  const pane = new Pane({ title: "Controls" });
  pane.element.classList.add("tp-minimal");
  pane.registerPlugin(EssentialsPlugin);
  const inputs = [];

  fpsGraph = pane.addBlade({
    view: "fpsgraph",
    label: "FPS",
    rows: 2,
  });

  const waveTypeOptions = {
    Cartesian: "Cartesian",
    Radial: "Radial",
    Spiral: "Spiral",
    HyperSpiral: "HyperSpiral",
    Parabolic: "Parabolic",
    Hexagonal: "Hexagonal",
    Quasicrystal: "Quasicrystal",
    Gyroid: "Gyroid",
    Lissajous: "Lissajous",
    Moire: "Moire",
    Chevron: "Chevron",
    Ring: "Ring",
    Cassini: "Cassini",
    Flower: "Flower",
    Interference: "Interference",
    GridWarp: "GridWarp",
    Diamond: "Diamond",
    Square: "Square",
    Lattice: "Lattice",
    Kaleidoscope: "Kaleidoscope",
    Voronoi: "Voronoi",
    Perlin: "Perlin",
    Superellipse: "Superellipse",
    Bessel: "Bessel",
    TriLattice: "TriLattice",
    Astroid: "Astroid",
    Checker: "Checker",
  };

  const particlesFolder = pane.addFolder({ title: "Particles" });
  inputs.push(
    particlesFolder
      .addBinding(CONFIG, "particleCount", {
        min: 13370,
        max: 1337000,
        step: 10000,
        label: "Count",
      })
      .on("change", () => rebuildParticles()),
  );
  inputs.push(
    particlesFolder
      .addBinding(CONFIG, "particleSize", {
        min: 1,
        max: 5,
        step: 0.5,
        label: "Size",
      })
      .on("change", (ev) => {
        points.material.size = ev.value;
      }),
  );
  inputs.push(
    particlesFolder
      .addBinding(CONFIG, "color", { label: "Color" })
      .on("change", (ev) => {
        points.material.color.set(ev.value);
      }),
  );
  inputs.push(
    particlesFolder
      .addBinding(CONFIG, "viewScale", {
        options: {
          200: 200,
          300: 300,
          400: 400,
          500: 500,
          600: 600,
          700: 700,
          800: 800,
        },
        label: "View scale",
      })
      .on("change", () => {
        onWindowResize();
        rebuildParticles();
        allocateField();
        rebuildField();
      }),
  );
  inputs.push(
    particlesFolder
      .addBinding(CONFIG, "gridSize", {
        options: {
          64: 64,
          128: 128,
          256: 256,
          384: 384,
          512: 512,
        },
        label: "Grid size",
      })
      .on("change", () => {
        allocateField();
        rebuildField();
      }),
  );

  const motionFolder = pane.addFolder({ title: "Motion" });
  inputs.push(
    motionFolder.addBinding(CONFIG, "settleStrength", {
      min: 0.1,
      max: 4.0,
      step: 0.01,
      label: "Settle strength",
    }),
  );
  inputs.push(
    motionFolder.addBinding(CONFIG, "jitter", {
      min: 0.0,
      max: 0.25,
      step: 0.01,
      label: "Jitter",
    }),
  );
  inputs.push(
    motionFolder.addBinding(CONFIG, "drag", {
      min: 0.7,
      max: 0.9,
      step: 0.01,
      label: "Drag",
    }),
  );
  inputs.push(
    motionFolder.addBinding(CONFIG, "speedLimit", {
      min: 0.5,
      max: 4.0,
      step: 0.1,
      label: "Speed limit",
    }),
  );

  const wavesFolder = pane.addFolder({ title: "Waves" });
  const waveAInput = wavesFolder
    .addBinding(CONFIG, "waveTypeA", {
      options: waveTypeOptions,
      label: "Wave type A",
    })
    .on("change", () => rebuildField());
  const waveBInput = wavesFolder
    .addBinding(CONFIG, "waveTypeB", {
      options: waveTypeOptions,
      label: "Wave type B",
    })
    .on("change", () => rebuildField());
  inputs.push(waveAInput, waveBInput);
  inputs.push(
    wavesFolder
      .addBinding(CONFIG, "waveMix", {
        min: 0,
        max: 1,
        step: 0.01,
        label: "Wave mix",
      })
      .on("change", () => rebuildField()),
  );

  const modesFolder = pane.addFolder({ title: "Modes" });

  inputs.push(
    modesFolder
      .addBinding(CONFIG, "modeCount", {
        min: 1,
        max: 20,
        step: 1,
        label: "Mode count",
      })
      .on("change", () => {
        if (suppressUIEvents) return;
        rebuildModesFromConfig();
      }),
  );
  inputs.push(
    modesFolder
      .addBinding(CONFIG, "mRange", {
        min: 0,
        max: 40,
        step: 0.1,
        label: "m range",
      })
      .on("change", () => {
        if (suppressUIEvents) return;
        rebuildModesFromConfig();
      }),
  );
  inputs.push(
    modesFolder
      .addBinding(CONFIG, "nRange", {
        min: 0,
        max: 40,
        step: 0.1,
        label: "n range",
      })
      .on("change", () => {
        if (suppressUIEvents) return;
        rebuildModesFromConfig();
      }),
  );

  inputs.push(
    modesFolder
      .addBinding(CONFIG, "integerModes", { label: "Integer modes" })
      .on("change", () => {
        if (suppressUIEvents) return;
        rebuildModesFromConfig();
      }),
  );
  inputs.push(
    modesFolder
      .addBinding(CONFIG, "randomPhase", { label: "Random phase" })
      .on("change", () => {
        if (suppressUIEvents) return;
        rebuildModesFromConfig();
      }),
  );
  inputs.push(
    modesFolder
      .addBinding(CONFIG, "randomRotation", { label: "Random rotation" })
      .on("change", () => {
        if (suppressUIEvents) return;
        rebuildModesFromConfig();
      }),
  );
  inputs.push(
    modesFolder
      .addBinding(CONFIG, "rotationMaxDeg", {
        min: 0,
        max: 90,
        step: 1,
        label: "Rotation max (deg)",
      })
      .on("change", () => {
        if (suppressUIEvents) return;
        rebuildModesFromConfig();
      }),
  );

  const captureFolder = pane.addFolder({ title: "Capture" });
  inputs.push(
    captureFolder.addBinding(CONFIG, "exportBaseSize", {
      min: 1024,
      max: 8192,
      step: 256,
      label: "Export size (px)",
    }),
  );
  inputs.push(
    captureFolder.addBinding(CONFIG, "exportOpaque", {
      label: "Transparent",
    }),
  );
  inputs.push(
    captureFolder
      .addBinding(CONFIG, "cameraControl", {
        label: "Pan/Zoom",
        x: { min: -CONFIG.viewScale, max: CONFIG.viewScale, step: 1 },
        y: { min: -CONFIG.viewScale, max: CONFIG.viewScale, step: 1 },
        z: { min: 0.5, max: 10, step: 0.1 },
      })
      .on("change", (ev) => {
        applyCameraFromControl(ev.value);
      }),
  );
  captureFolder
    .addButton({ title: "Save image" })
    .on("click", () => saveImage());

  if (typeof captureFolder.addSeparator === "function") {
    captureFolder.addSeparator();
  } else if (typeof captureFolder.addBlade === "function") {
    captureFolder.addBlade({ view: "separator" });
  }
  captureFolder
    .addButton({ title: "Randomize all" })
    .on("click", () => randomizeAll());
  captureFolder
    .addButton({ title: "Randomize modes" })
    .on("click", () => randomizeModes());

  refreshUI = () => {
    inputs.forEach((input) => input.refresh());
  };

  particlesFolder.expanded = true;
  captureFolder.expanded = true;
}

function saveImage() {
  const oldSize = new THREE.Vector2();
  renderer.getSize(oldSize);
  const oldPixelRatio = renderer.getPixelRatio();
  const oldClearAlpha = renderer.getClearAlpha();
  const oldClearColor = renderer.getClearColor(new THREE.Color());
  const oldPointSize = points.material.size;

  const baseSize = Math.max(256, Math.round(CONFIG.exportBaseSize));

  const aspect = oldSize.x / oldSize.y;
  const exportWidth = Math.max(1, Math.round(baseSize));
  const exportHeight = Math.max(1, Math.round(exportWidth / aspect));

  renderer.setPixelRatio(1);
  renderer.setSize(exportWidth, exportHeight, false);
  const exportAlpha = CONFIG.exportOpaque ? 0 : 1;
  renderer.setClearColor(oldClearColor, exportAlpha);

  const frustumSize = CONFIG.viewScale * 2;
  camera.left = (-frustumSize * aspect) / 2;
  camera.right = (frustumSize * aspect) / 2;
  camera.top = frustumSize / 2;
  camera.bottom = -frustumSize / 2;
  camera.updateProjectionMatrix();

  const sizeScale = exportWidth / oldSize.x;
  points.material.size = oldPointSize * sizeScale;

  renderer.render(scene, camera);

  const link = document.createElement("a");
  link.download = `chladni-${Date.now()}-${exportWidth}x${exportHeight}.png`;
  link.href = renderer.domElement.toDataURL("image/png");
  link.click();

  points.material.size = oldPointSize;
  renderer.setClearColor(oldClearColor, oldClearAlpha);
  renderer.setPixelRatio(oldPixelRatio);
  renderer.setSize(oldSize.x, oldSize.y, false);
  onWindowResize();
}

function randomizeAll() {
  const waveTypes = [
    "Cartesian",
    "Radial",
    "Spiral",
    "HyperSpiral",
    "Parabolic",
    "Hexagonal",
    "Quasicrystal",
    "Gyroid",
    "Lissajous",
    "Moire",
    "Chevron",
    "Ring",
    "Cassini",
    "Flower",
    "Interference",
    "GridWarp",
    "Diamond",
    "Square",
    "Lattice",
    "Kaleidoscope",
    "Voronoi",
    "Perlin",
    "Superellipse",
    "Bessel",
    "TriLattice",
    "Astroid",
    "Checker",
  ];

  const deviceMemory = navigator.deviceMemory || 8;
  const cores = navigator.hardwareConcurrency || 8;
  const perfScale = Math.min(
    1,
    Math.max(0.5, (deviceMemory / 8) * (cores / 8)),
  );
  const minParticles = 100000;
  const maxParticles = Math.max(minParticles, Math.round(400000 * perfScale));
  const particleStep = 20000;
  const minGrid = 64;
  const maxGrid = Math.max(minGrid, Math.round(192 * perfScale));
  const maxModes = Math.max(4, Math.round(10 * perfScale));

  const particleOptions = [];
  for (let v = minParticles; v <= maxParticles; v += particleStep) {
    particleOptions.push(v);
  }
  CONFIG.particleCount =
    particleOptions[Math.floor(Math.random() * particleOptions.length)];
  CONFIG.settleStrength = 0.1 + Math.random() * 3.9;
  CONFIG.jitter = Math.random() * 0.25;
  CONFIG.drag = 0.7 + Math.random() * 0.2;
  CONFIG.speedLimit = 0.5 + Math.random() * 9.5;

  CONFIG.waveTypeA = waveTypes[Math.floor(Math.random() * waveTypes.length)];
  CONFIG.waveTypeB = waveTypes[Math.floor(Math.random() * waveTypes.length)];
  CONFIG.waveMix = Math.random();

  CONFIG.particleSize = 1 + Math.random() * 2;

  CONFIG.modeCount = Math.floor(1 + Math.random() * maxModes);
  CONFIG.mRange.min = Math.random() * 20;
  CONFIG.mRange.max =
    CONFIG.mRange.min + Math.random() * (40 - CONFIG.mRange.min);
  CONFIG.nRange.min = Math.random() * 20;
  CONFIG.nRange.max =
    CONFIG.nRange.min + Math.random() * (40 - CONFIG.nRange.min);

  CONFIG.integerModes = Math.random() < 0.7;
  CONFIG.randomPhase = Math.random() < 0.5;
  CONFIG.randomRotation = Math.random() < 0.5;
  CONFIG.rotationMaxDeg = Math.floor(Math.random() * 91);

  rebuildParticles();
  allocateField();
  initModes();
  rebuildField();

  points.material.color.set(CONFIG.color);
  points.material.size = CONFIG.particleSize;

  applyCameraFromControl(CONFIG.cameraControl);
  onWindowResize();

  if (refreshUI) refreshUI();
}
