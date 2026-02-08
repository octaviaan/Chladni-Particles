import * as THREE from "three";
import { Pane } from "tweakpane";
import * as EssentialsPlugin from "tweakpane-plugin-essentials";

// --- Configuration ---
const CONFIG = {
  particleCount: 300000,
  gridSize: 256,
  settleStrength: 2.0,
  jitter: 0.1,
  drag: 0.85,
  speedLimit: 2.0,
  viewScale: 600,
  color: "#b3a79b",
  particleSize: 2.0,
  particleOpacity: 1.0,
  cameraControl: { x: 0, y: 0, z: 1 },
  cameraDefault: { x: 0, y: 0, z: 1 },
  exportBaseSize: 3000,
  exportOpaque: false,

  modeCount: 4,
  mRange: { min: 2, max: 4 },
  nRange: { min: 4, max: 8 },

  rectAspect: 1,
  fieldScale: 1.0,

  waveTypeA: "Cartesian",
  waveTypeB: "Radial",
  waveMix: 0.5,
  spatialMixNoise: 0.35,

  integerModes: true,
};

// --- Math Helpers ---
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

// --- Wave Strategy Definitions ---
const WAVE_FUNCTIONS = {
  Beat: (cx, cy, mode) => {
    const ax = 0.35;
    const ay = 0.1;
    const d1 = Math.sqrt((cx + ax) ** 2 + (cy + ay) ** 2);
    const d2 = Math.sqrt((cx - ax) ** 2 + (cy - ay) ** 2);
    const k1 = mode.m * 16;
    const k2 = (mode.n + 0.35) * 16;
    return Math.sin(k1 * d1 + mode.px) + Math.sin(k2 * d2 + mode.py);
  },
  Biharmonic: (cx, cy, mode) => {
    const rx = cx * mode.cos - cy * mode.sin;
    const ry = cx * mode.sin + cy * mode.cos;
    const a = Math.max(1, mode.m);
    const b = Math.max(1, mode.n);
    return (
      Math.sin(a * Math.PI * rx + mode.px) *
        Math.sin(b * Math.PI * ry + mode.py) +
      0.5 *
        Math.sin((a + 1) * Math.PI * rx + mode.py) *
        Math.sin((b + 2) * Math.PI * ry + mode.px)
    );
  },
  Cartesian: (cx, cy, mode) => {
    const rx = cx * mode.cos - cy * mode.sin;
    const ry = cx * mode.sin + cy * mode.cos;
    return (
      Math.sin(mode.m * Math.PI * (rx + 0.5) + mode.px) *
      Math.sin(mode.n * Math.PI * (ry + 0.5) + mode.py)
    );
  },
  Cassini: (cx, cy, mode) => {
    const a = mode.m * 0.1;
    const d1 = (cx - a) * (cx - a) + cy * cy;
    const d2 = (cx + a) * (cx + a) + cy * cy;
    return Math.sin(Math.sqrt(d1 * d2) * 20 - mode.n * 5);
  },
  Chebyshev: (cx, cy, mode) => {
    const r = Math.sqrt(cx * cx + cy * cy) * 2.0;
    const n = Math.max(2, Math.round(mode.m));
    return Math.cos(n * Math.acos(Math.max(-1, Math.min(1, r))) + mode.px);
  },
  Chevron: (cx, cy, mode) => {
    const rx = cx * mode.cos - cy * mode.sin;
    const ry = cx * mode.sin + cy * mode.cos;
    return (
      Math.sin(mode.m * Math.PI * (rx + ry) + mode.px) *
      Math.sin(mode.n * Math.PI * (rx - ry) + mode.py)
    );
  },
  Crosshatch: (cx, cy, mode) => {
    const rx = cx * mode.cos - cy * mode.sin;
    const ry = cx * mode.sin + cy * mode.cos;
    return (
      Math.sin(mode.m * Math.PI * rx + mode.px) +
      Math.sin(mode.n * Math.PI * ry + mode.py)
    );
  },
  Diamond: (cx, cy, mode) => {
    return Math.sin(
      mode.m * Math.PI * (Math.abs(cx) + Math.abs(cy)) * 2.0 + mode.px,
    );
  },
  GridWarp: (cx, cy, mode) => {
    const rx = cx * mode.cos - cy * mode.sin;
    const ry = cx * mode.sin + cy * mode.cos;
    const warp = 0.35;
    const u = rx + warp * Math.sin(mode.n * Math.PI * ry + mode.py);
    const v = ry + warp * Math.sin(mode.m * Math.PI * rx + mode.px);
    return (
      Math.sin(mode.m * Math.PI * u + mode.px) *
      Math.sin(mode.n * Math.PI * v + mode.py)
    );
  },
  Gyroid: (cx, cy, mode) => {
    const scaleX = mode.m * 10;
    const scaleY = Math.max(0.5, mode.n) * 10;
    return (
      Math.sin(cx * scaleX) * Math.cos(cy * scaleY) +
      Math.sin(cy * scaleY) * Math.cos(mode.px)
    );
  },
  Hexagonal: (cx, cy, mode) => {
    const rx = cx * mode.cos - cy * mode.sin;
    const ry = cx * mode.sin + cy * mode.cos;
    const kFreq = mode.m * Math.PI;
    const v1 = rx;
    const v2 = -0.5 * rx + 0.866 * ry;
    const v3 = -0.5 * rx - 0.866 * ry;
    return (
      Math.sin(kFreq * v1 + mode.px) +
      Math.sin(kFreq * v2 + mode.py) +
      Math.sin(kFreq * v3)
    );
  },
  Interference: (cx, cy, mode) => {
    const d1 = Math.sqrt((cx + 0.3) ** 2 + cy ** 2);
    const d2 = Math.sqrt((cx - 0.3) ** 2 + cy ** 2);
    return Math.sin(mode.m * 20 * d1) + Math.sin(mode.m * 20 * d2 + mode.px);
  },
  Lattice: (cx, cy, mode) => {
    const rx = cx * mode.cos - cy * mode.sin;
    const ry = cx * mode.sin + cy * mode.cos;
    const k = mode.m * Math.PI;
    const v1 = rx;
    const v2 = -0.5 * rx + 0.866 * ry;
    const v3 = -0.5 * rx - 0.866 * ry;
    return (
      Math.cos(k * v1 + mode.px) * Math.cos(k * v2 + mode.py) * Math.cos(k * v3)
    );
  },
  Lissajous: (cx, cy, mode) => {
    const rx = cx * mode.cos - cy * mode.sin;
    const ry = cx * mode.sin + cy * mode.cos;
    return (
      Math.sin(mode.m * Math.PI * rx + mode.px) *
      Math.sin(mode.n * Math.PI * ry + mode.py)
    );
  },
  MultiSource: (cx, cy, mode) => {
    const k = mode.m * 14;
    const sources = [
      [0.35, 0.0],
      [-0.35, 0.0],
      [0.0, 0.3],
      [0.0, -0.3],
    ];
    let sum = 0;
    for (let i = 0; i < sources.length; i++) {
      const sx = sources[i][0];
      const sy = sources[i][1];
      const d = Math.sqrt((cx - sx) ** 2 + (cy - sy) ** 2);
      sum += Math.sin(k * d + mode.px + i * 0.7);
    }
    return sum / sources.length;
  },
  Octagon: (cx, cy, mode) => {
    const x = Math.abs(cx);
    const y = Math.abs(cy);
    const d1 = (x + y) * 0.7071;
    const d2 = Math.abs(x - y) * 0.7071;
    const r = Math.max(x, y, d1, d2) * 2.0;
    return Math.sin(mode.m * Math.PI * r + mode.px);
  },
  Parabolic: (cx, cy, mode) => {
    const u = cx * cx - cy * cy;
    const v = 2 * cx * cy;
    return Math.sin(mode.m * Math.PI * u) * Math.sin(mode.n * Math.PI * v);
  },
  Quasicrystal: (cx, cy, mode) => {
    const rx = cx * mode.cos - cy * mode.sin;
    const ry = cx * mode.sin + cy * mode.cos;
    const kFreq = mode.m * Math.PI;
    let sum = 0;
    for (let j = 0; j < 5; j++) {
      const theta = (Math.PI * 2 * j) / 5;
      const v = rx * Math.cos(theta) + ry * Math.sin(theta);
      sum += Math.cos(kFreq * v + mode.px);
    }
    return sum;
  },
  Radial: (cx, cy, mode) => {
    const r = Math.sqrt(cx * cx + cy * cy) * 2.0;
    const theta = Math.atan2(cy, cx);
    return (
      Math.sin(mode.n * Math.PI * r + mode.px) *
      Math.cos(Math.round(mode.m) * theta + mode.py)
    );
  },
  RingInterference: (cx, cy, mode) => {
    const r = Math.sqrt(cx * cx + cy * cy) * 2.0;
    const ax = 0.35;
    const d1 = Math.sqrt((cx + ax) ** 2 + cy * cy) * 2.0;
    const d2 = Math.sqrt((cx - ax) ** 2 + cy * cy) * 2.0;
    const k = mode.m * Math.PI * 2;
    return (
      Math.sin(k * r + mode.px) +
      0.5 * Math.sin(k * d1 + mode.py) +
      0.5 * Math.sin(k * d2 + mode.py)
    );
  },
  Rose: (cx, cy, mode) => {
    const r = Math.sqrt(cx * cx + cy * cy) * 2.0;
    const theta = Math.atan2(cy, cx);
    const k = Math.max(2, Math.round(mode.m));
    return Math.sin(k * theta + mode.px) * Math.cos(mode.n * Math.PI * r);
  },
  Rosette: (cx, cy, mode) => {
    const N = Math.max(3, Math.round(mode.m));
    let sum = 0;
    for (let i = 0; i < N; i++) {
      const theta = (Math.PI * 2 * i) / N;
      const rx = cx * Math.cos(theta) - cy * Math.sin(theta);
      sum += Math.cos(rx * mode.n * 20 + mode.px);
    }
    return sum / N;
  },
  Spiral: (cx, cy, mode) => {
    const r = Math.sqrt(cx * cx + cy * cy) * 2.0;
    const theta = Math.atan2(cy, cx);
    return Math.sin(
      mode.n * Math.PI * r + Math.round(mode.m) * theta + mode.px,
    );
  },
  Square: (cx, cy, mode) => {
    return Math.sin(
      mode.m * Math.PI * Math.max(Math.abs(cx), Math.abs(cy)) * 2.0 + mode.px,
    );
  },
  Stripe: (cx, cy, mode) => {
    const rx = cx * mode.cos - cy * mode.sin;
    return Math.sin(mode.m * Math.PI * rx + mode.px);
  },
  Superellipse: (cx, cy, mode) => {
    const p = 0.8 + (Math.abs(mode.m) % 4) * 0.4;
    const r =
      Math.pow(
        Math.pow(Math.abs(cx) * 2, p) + Math.pow(Math.abs(cy) * 2, p),
        1 / p,
      ) * 0.5;
    return Math.sin(mode.n * Math.PI * r + mode.px);
  },
  Voronoi: (cx, cy, mode) => {
    let minD = 1e9;
    for (let j = 0; j < 6; j++) {
      const sx = hash2(j + mode.m * 1.3, mode.n * 2.1) - 0.5;
      const sy = hash2(j + mode.n * 1.7, mode.m * 2.5) - 0.5;
      const d = Math.sqrt((cx - sx) ** 2 + (cy - sy) ** 2);
      if (d < minD) minD = d;
    }
    return Math.sin(mode.m * Math.PI * minD * 2 + mode.px);
  },
  Wallpaper: (cx, cy, mode) => {
    const scale = 3 + (mode.m % 5);
    const wx = (cx * scale + 100) % 2.0;
    const wy = (cy * scale + 100) % 2.0;
    const fx = Math.abs(wx - 1.0);
    const fy = Math.abs(wy - 1.0);
    return (
      Math.sin(fx * Math.PI * mode.n + mode.px) *
      Math.cos(fy * Math.PI * mode.n + mode.py)
    );
  },
};

const WAVE_TYPE_KEYS = Object.keys(WAVE_FUNCTIONS);

// --- Global Variables ---
let scene, camera, renderer, geometry, points;
let positions, velocities, colors;
let refreshUI = null;
let suppressUIEvents = false;
let isPanning = false;
let lastPointer = { x: 0, y: 0 };

// Field Data
let energy;
let gradX;
let gradY;

// Modes
let modes = [];

// --- Initialization ---
init();
allocateField();
initModes();
rebuildField();
animate();

function init() {
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
  camera.position.set(0, 0, 500);
  applyCameraFromControl(CONFIG.cameraControl);

  renderer = new THREE.WebGLRenderer({
    antialias: true,
    powerPreference: "high-performance",
    preserveDrawingBuffer: true,
  });
  renderer.setSize(window.innerWidth, window.innerHeight);
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  renderer.setClearColor(0x000000, 1);
  document.body.appendChild(renderer.domElement);

  buildParticles();

  const sprite = new THREE.TextureLoader().load(
    "https://threejs.org/examples/textures/sprites/disc.png",
  );

  const material = new THREE.PointsMaterial({
    color: 0xffffff,
    vertexColors: true,
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

  window.addEventListener("resize", onWindowResize);
  setupMouseControls();
  window.addEventListener("keydown", onKeyDown);

  setupGUI();
}

// --- Physics & Math ---
function allocateField() {
  const size = CONFIG.gridSize * CONFIG.gridSize;
  energy = new Float32Array(size);
  gradX = new Float32Array(size);
  gradY = new Float32Array(size);
}

function buildParticles() {
  geometry = new THREE.BufferGeometry();

  positions = new Float32Array(CONFIG.particleCount * 3);
  velocities = new Float32Array(CONFIG.particleCount * 2);
  colors = new Float32Array(CONFIG.particleCount * 3);

  const baseColor = new THREE.Color(CONFIG.color);

  for (let i = 0; i < CONFIG.particleCount; i++) {
    const [x, y] = randomPointInShape();
    positions[i * 3 + 0] = x;
    positions[i * 3 + 1] = y;
    positions[i * 3 + 2] = 0;

    velocities[i * 2 + 0] = 0;
    velocities[i * 2 + 1] = 0;

    colors[i * 3 + 0] = baseColor.r;
    colors[i * 3 + 1] = baseColor.g;
    colors[i * 3 + 2] = baseColor.b;
  }

  geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
  geometry.setAttribute("color", new THREE.BufferAttribute(colors, 3));
}

function rebuildParticles() {
  const oldGeometry = points.geometry;
  buildParticles();
  points.geometry = geometry;
  oldGeometry.dispose();
}

function applyParticleColor(hex) {
  if (!colors || !geometry?.attributes?.color) return;
  const nextColor = new THREE.Color(hex);
  const count = CONFIG.particleCount;
  for (let i = 0; i < count; i++) {
    const i3 = i * 3;
    colors[i3] = nextColor.r;
    colors[i3 + 1] = nextColor.g;
    colors[i3 + 2] = nextColor.b;
  }
  geometry.attributes.color.needsUpdate = true;
}

function rebuildField() {
  const G = CONFIG.gridSize;
  const funcA = WAVE_FUNCTIONS[CONFIG.waveTypeA] || WAVE_FUNCTIONS["Cartesian"];
  const funcB = WAVE_FUNCTIONS[CONFIG.waveTypeB] || WAVE_FUNCTIONS["Cartesian"];
  const baseBias = CONFIG.waveMix - 0.5;
  const fieldScale = Math.max(0.05, CONFIG.fieldScale || 1);
  const maxR = Math.SQRT1_2 * fieldScale;
  const spatialMixScale = 3.0;

  for (let y = 0; y < G; y++) {
    for (let x = 0; x < G; x++) {
      const tx = x / (G - 1);
      const ty = y / (G - 1);
      const cx = (tx - 0.5) * fieldScale;
      const cy = (ty - 0.5) * fieldScale;

      const idx = y * G + x;
      let phi = 0;
      const r = Math.sqrt(cx * cx + cy * cy);
      const rn = clamp(r / maxR, 0, 1);
      const noiseWeight = smoothstep(CONFIG.spatialMixNoise);
      const noiseRange = lerp(1.0, 2.0, CONFIG.spatialMixNoise);
      const noiseRaw = valueNoise(cx * spatialMixScale, cy * spatialMixScale);
      const noise = clamp(0.5 + 0.5 * noiseRaw * noiseRange, 0, 1);
      const mixed = (1 - noiseWeight) * rn + noiseWeight * noise + baseBias;
      const spatialMix = smoothstep(clamp(mixed, 0, 1));

      for (let k = 0; k < modes.length; k++) {
        const mode = modes[k];
        const wa = Math.tanh(funcA(cx, cy, mode));
        const wb = Math.tanh(funcB(cx, cy, mode));
        const wave = wa * (1 - spatialMix) + wb * spatialMix;
        phi += mode.a * wave;
      }
      energy[idx] = phi * phi;
    }
  }

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

  const settle = CONFIG.settleStrength;
  const jitter = CONFIG.jitter;
  const drag = CONFIG.drag;
  const limit = CONFIG.speedLimit;

  const gridScale = (G - 1) / fullRange;
  const rectRy = Math.max(1, Math.round(range / CONFIG.rectAspect));

  for (let i = 0; i < count; i++) {
    const i3 = i * 3;
    const i2 = i * 2;

    let x = pos[i3];
    let y = pos[i3 + 1];
    let vx = vel[i2];
    let vy = vel[i2 + 1];

    let inside = true;
    if (x < -range || x > range || y < -rectRy || y > rectRy) inside = false;

    if (!inside) {
      const respawn = randomPointInShape(Math.random);
      x = respawn[0];
      y = respawn[1];
      vx = 0;
      vy = 0;
    }

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

    vx -= gxVal * settle;
    vy -= gyVal * settle;

    vx += (Math.random() - 0.5) * jitter;
    vy += (Math.random() - 0.5) * jitter;

    vx *= drag;
    vy *= drag;

    const speedSq = vx * vx + vy * vy;
    const speed = Math.sqrt(speedSq);

    if (speed > limit) {
      const scale = limit / speed;
      vx *= scale;
      vy *= scale;
    }

    x += vx;
    y += vy;

    pos[i3] = x;
    pos[i3 + 1] = y;
    pos[i3 + 2] = 0;

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

function randomPointInShape(randomFn = Math.random) {
  const range = CONFIG.viewScale;
  const rectRy = Math.max(1, Math.round(range / CONFIG.rectAspect));
  return [(randomFn() - 0.5) * 2 * range, (randomFn() - 0.5) * 2 * rectRy];
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
  const maxModes = 12;
  CONFIG.modeCount = Math.floor(1 + Math.random() * maxModes);
  CONFIG.mRange.min = Math.random() * 20;
  CONFIG.mRange.max =
    CONFIG.mRange.min + Math.random() * (15 - CONFIG.mRange.min);
  CONFIG.nRange.min = Math.random() * 20;
  CONFIG.nRange.max =
    CONFIG.nRange.min + Math.random() * (15 - CONFIG.nRange.min);

  CONFIG.integerModes = Math.random() < 0.5;

  initModes();
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
    modes.push({
      m: quantizeModeValue(mRaw),
      n: quantizeModeValue(nRaw),
      a: Math.exp(-i * 0.35),
      px: 0,
      py: 0,
      cos: 1,
      sin: 0,
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
  updateParticles();
  renderer.render(scene, camera);
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

function clamp(value, min, max) {
  return Math.max(min, Math.min(max, value));
}

function setupMouseControls() {
  const canvas = renderer.domElement;

  canvas.addEventListener(
    "wheel",
    (event) => {
      event.preventDefault();
      const zoomFactor = Math.exp(-event.deltaY * 0.0015);
      CONFIG.cameraControl.z = clamp(
        CONFIG.cameraControl.z * zoomFactor,
        0.5,
        10,
      );
      applyCameraZoom(CONFIG.cameraControl.z);
      if (refreshUI) refreshUI();
    },
    { passive: false },
  );

  canvas.addEventListener("pointerdown", (event) => {
    if (event.button !== 0) return;
    isPanning = true;
    lastPointer.x = event.clientX;
    lastPointer.y = event.clientY;
    canvas.setPointerCapture(event.pointerId);
  });

  canvas.addEventListener("dblclick", (event) => {
    event.preventDefault();
    CONFIG.cameraControl.x = CONFIG.cameraDefault.x;
    CONFIG.cameraControl.y = CONFIG.cameraDefault.y;
    CONFIG.cameraControl.z = CONFIG.cameraDefault.z;
    applyCameraFromControl(CONFIG.cameraControl);
    if (refreshUI) refreshUI();
  });

  canvas.addEventListener("pointermove", (event) => {
    if (!isPanning) return;
    const deltaX = event.clientX - lastPointer.x;
    const deltaY = event.clientY - lastPointer.y;
    lastPointer.x = event.clientX;
    lastPointer.y = event.clientY;

    const viewWidth = (camera.right - camera.left) / camera.zoom;
    const viewHeight = (camera.top - camera.bottom) / camera.zoom;
    const worldPerPixelX = viewWidth / canvas.clientWidth;
    const worldPerPixelY = viewHeight / canvas.clientHeight;

    CONFIG.cameraControl.x = clamp(
      CONFIG.cameraControl.x - deltaX * worldPerPixelX,
      -CONFIG.viewScale,
      CONFIG.viewScale,
    );
    CONFIG.cameraControl.y = clamp(
      CONFIG.cameraControl.y + deltaY * worldPerPixelY,
      -CONFIG.viewScale,
      CONFIG.viewScale,
    );
    applyCameraPan(CONFIG.cameraControl.x, CONFIG.cameraControl.y);
    if (refreshUI) refreshUI();
  });

  const stopPan = (event) => {
    if (!isPanning) return;
    isPanning = false;
    if (canvas.hasPointerCapture(event.pointerId))
      canvas.releasePointerCapture(event.pointerId);
  };
  canvas.addEventListener("pointerup", stopPan);
  canvas.addEventListener("pointercancel", stopPan);
  canvas.addEventListener("pointerleave", () => {
    isPanning = false;
  });
}

function onKeyDown(event) {
  const tag = event.target ? event.target.tagName.toLowerCase() : "";
  if (tag === "input" || tag === "textarea" || tag === "select") return;
  if (event.key.toLowerCase() === "r") randomizeAll();
}

function setupGUI() {
  const pane = new Pane({ title: "Controls" });
  pane.element.classList.add("tp-minimal");
  pane.registerPlugin(EssentialsPlugin);
  const inputs = [];

  const titleEl = pane.element.querySelector(".tp-rotv_t");

  const legendStyle = document.createElement("style");
  legendStyle.textContent = `
    .tp-legend{
      margin: 4px 14px 4px;
      padding: 2px 2px;
      color: #afcea1;
      font: 12px/1.2 ui-sans-serif, system-ui, -apple-system, Segoe UI, Roboto;
      display: grid;
      gap: 4px;
    }
    .tp-legend-row{ display:flex; align-items:center; gap:8px; white-space:nowrap; }
    .tp-legend-ico{ width:24px; height:24px; display:inline-flex; opacity:.95; }
    .tp-legend svg{ display:block; }
  `;
  document.head.appendChild(legendStyle);

  const svgs = {
    wheel: `
    <svg width="24" height="24" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
<path d="M12 21C8.68629 21 6 18.3137 6 15V9C6 5.68629 8.68629 3 12 3C15.3137 3 18 5.68629 18 9V15C18 18.3137 15.3137 21 12 21Z" stroke="#B8B8B8" stroke-linecap="round" stroke-linejoin="round"/>
<path d="M12 10V8" stroke="#AFCEA1" stroke-width="3" stroke-linecap="round" stroke-linejoin="round"/>
</svg>

`,
    drag: `
    <svg width="24" height="24" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
<path d="M12.0049 8.75V14.75" stroke="#B8B8B8"/>
<path d="M15.0049 12L9.00488 12" stroke="#B8B8B8"/>
<path d="M12 3L14.5981 6.75H9.40192L12 3Z" fill="#AFCEA1"/>
<path d="M12 21L9.40192 17.25L14.5981 17.25L12 21Z" fill="#AFCEA1"/>
<path d="M21 12L17.25 14.5981L17.25 9.40192L21 12Z" fill="#AFCEA1"/>
<path d="M3 12L6.75 9.40192L6.75 14.5981L3 12Z" fill="#AFCEA1"/>
</svg>

`,
    dbl: `
    <svg width="24" height="24" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
<path d="M12 21C8.68629 21 6 18.3137 6 15V9C6 5.68629 8.68629 3 12 3C15.3137 3 18 5.68629 18 9V15C18 18.3137 15.3137 21 12 21Z" stroke="#B8B8B8" stroke-linecap="round" stroke-linejoin="round"/>
<circle cx="8.5" cy="6.5" r="3.5" fill="#AFCEA1"/>
</svg>
`,
    keyR: `
    <svg width="24" height="24" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
<rect x="2.5" y="2.5" width="19" height="19" rx="2.5" stroke="#B8B8B8"/>
<path d="M10.7831 16H9V8H12.253C12.6627 8 13.0281 8.06113 13.3494 8.18338C13.6707 8.29799 13.9398 8.46609 14.1566 8.68768C14.3815 8.90926 14.5502 9.1767 14.6627 9.48997C14.7751 9.80325 14.8313 10.1547 14.8313 10.5444C14.8313 11.1098 14.7028 11.595 14.4458 12C14.1888 12.405 13.8233 12.6724 13.3494 12.8023L15 16H13.0241L11.6145 13.0315H10.7831V16ZM11.8193 11.702C12.245 11.702 12.5382 11.6256 12.6988 11.4728C12.8675 11.32 12.9518 11.0678 12.9518 10.7163V10.3725C12.9518 10.021 12.8675 9.76886 12.6988 9.61605C12.5382 9.46323 12.245 9.38682 11.8193 9.38682H10.7831V11.702H11.8193Z" fill="#AFCEA1"/>
</svg>

`,
  };

  const legend = document.createElement("div");
  legend.className = "tp-legend";
  legend.innerHTML = `
    <div class="tp-legend-row"><span class="tp-legend-ico">${svgs.wheel}</span><span>Wheel: zoom</span></div>
    <div class="tp-legend-row"><span class="tp-legend-ico">${svgs.drag}</span><span>Drag: pan</span></div>
    <div class="tp-legend-row"><span class="tp-legend-ico">${svgs.dbl}</span><span>Double click: reset postion</span></div>
    <div class="tp-legend-row"><span class="tp-legend-ico">${svgs.keyR}</span><span>R: randomize all</span></div>
  `;
  if (titleEl && titleEl.parentElement) {
    titleEl.parentElement.insertAdjacentElement("afterend", legend);
  } else {
    pane.element.prepend(legend);
  }

  const waveTypeOptions = WAVE_TYPE_KEYS.reduce((obj, key) => {
    obj[key] = key;
    return obj;
  }, {});

  const particlesFolder = pane.addFolder({ title: "Particles" });
  inputs.push(
    particlesFolder
      .addBinding(CONFIG, "particleCount", {
        min: 50000,
        max: 500000,
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
      .on("change", (ev) => applyParticleColor(ev.value)),
  );
  particlesFolder
    .addButton({ title: "Rebuild particles" })
    .on("click", () => rebuildParticles());

  const motionFolder = pane.addFolder({ title: "Motion" });
  inputs.push(
    motionFolder.addBinding(CONFIG, "settleStrength", {
      min: 0.5,
      max: 10.0,
      step: 0.1,
    }),
  );
  inputs.push(
    motionFolder.addBinding(CONFIG, "jitter", {
      min: 0.0,
      max: 0.3,
      step: 0.01,
    }),
  );
  inputs.push(
    motionFolder.addBinding(CONFIG, "drag", { min: 0.7, max: 0.9, step: 0.01 }),
  );
  inputs.push(
    motionFolder.addBinding(CONFIG, "speedLimit", {
      min: 0.2,
      max: 3.0,
      step: 0.1,
    }),
  );

  const wavesFolder = pane.addFolder({ title: "Waves" });
  inputs.push(
    wavesFolder
      .addBinding(CONFIG, "waveTypeA", {
        options: waveTypeOptions,
        label: "Wave A",
      })
      .on("change", () => rebuildField()),
  );
  inputs.push(
    wavesFolder
      .addBinding(CONFIG, "waveTypeB", {
        options: waveTypeOptions,
        label: "Wave B",
      })
      .on("change", () => rebuildField()),
  );
  inputs.push(
    wavesFolder
      .addBinding(CONFIG, "waveMix", {
        min: 0.0,
        max: 1.0,
        step: 0.01,
        label: "Mix",
      })
      .on("change", () => rebuildField()),
  );
  inputs.push(
    wavesFolder
      .addBinding(CONFIG, "spatialMixNoise", {
        min: 0.0,
        max: 1.0,
        step: 0.01,
        label: "Noise",
      })
      .on("change", () => rebuildField()),
  );
  inputs.push(
    wavesFolder
      .addBinding(CONFIG, "fieldScale", {
        min: 0.5,
        max: 5.0,
        step: 0.1,
        label: "Field scale",
      })
      .on("change", () => rebuildField()),
  );
  wavesFolder
    .addButton({ title: "Randomize waves" })
    .on("click", () => randomizeWaves());

  const modesFolder = pane.addFolder({ title: "Modes" });
  inputs.push(
    modesFolder
      .addBinding(CONFIG, "modeCount", { min: 1, max: 10, step: 1 })
      .on("change", () => {
        if (!suppressUIEvents) rebuildModesFromConfig();
      }),
  );
  inputs.push(
    modesFolder
      .addBinding(CONFIG, "mRange", { min: 1, max: 10, step: 0.1 })
      .on("change", () => {
        if (!suppressUIEvents) rebuildModesFromConfig();
      }),
  );
  inputs.push(
    modesFolder
      .addBinding(CONFIG, "nRange", { min: 1, max: 10, step: 0.1 })
      .on("change", () => {
        if (!suppressUIEvents) rebuildModesFromConfig();
      }),
  );

  const captureFolder = pane.addFolder({ title: "Capture" });
  inputs.push(
    captureFolder.addBinding(CONFIG, "exportBaseSize", {
      min: 1000,
      max: 9000,
      step: 500,
    }),
  );
  inputs.push(
    captureFolder
      .addBinding(CONFIG, "rectAspect", {
        view: "radiogrid",
        groupName: "rect-aspect",
        size: [2, 1],
        cells: (x) => {
          const options = [
            { title: "Square", value: 1 },
            { title: "Landscape", value: 1.78 },
          ];
          return options[x];
        },
        label: "Aspect",
      })
      .on("change", () => rebuildParticles()),
  );
  inputs.push(
    captureFolder.addBinding(CONFIG, "exportOpaque", {
      label: "Transparent background",
    }),
  );
  inputs.push(
    captureFolder
      .addBinding(CONFIG, "cameraControl", {
        x: { min: -CONFIG.viewScale, max: CONFIG.viewScale },
        y: { min: -CONFIG.viewScale, max: CONFIG.viewScale },
        z: { min: 0.5, max: 10 },
      })
      .on("change", (ev) => applyCameraFromControl(ev.value)),
  );
  captureFolder
    .addButton({ title: "Save image" })
    .on("click", () => saveImage({ useViewport: false }));
  captureFolder
    .addButton({ title: "Save current view" })
    .on("click", () => saveImage({ useViewport: true }));

  refreshUI = () => inputs.forEach((input) => input.refresh());
  particlesFolder.expanded = true;
  captureFolder.expanded = true;
}

function saveImage({ useViewport = false } = {}) {
  const oldSize = new THREE.Vector2();
  renderer.getSize(oldSize);
  const oldPixelRatio = renderer.getPixelRatio();
  const oldClearColor = renderer.getClearColor(new THREE.Color());
  const oldPointSize = points.material.size;
  const oldCameraControl = {
    x: CONFIG.cameraControl.x,
    y: CONFIG.cameraControl.y,
    z: CONFIG.cameraControl.z,
  };

  const baseSize = Math.max(256, Math.round(CONFIG.exportBaseSize));
  const viewAspect = oldSize.x / oldSize.y;
  const aspect = useViewport
    ? viewAspect
    : Math.max(0.1, CONFIG.rectAspect || viewAspect);
  const exportWidth = Math.max(1, Math.round(baseSize));
  const exportHeight = Math.max(1, Math.round(exportWidth / aspect));

  renderer.setPixelRatio(1);
  renderer.setSize(exportWidth, exportHeight, false);
  renderer.setClearColor(oldClearColor, CONFIG.exportOpaque ? 0 : 1);

  if (!useViewport) {
    const range = CONFIG.viewScale;
    const rectRy = Math.max(1, Math.round(range / CONFIG.rectAspect));
    let halfWidth = range;
    let halfHeight = rectRy;
    const boundsAspect = halfWidth / halfHeight;
    if (boundsAspect > aspect) {
      halfHeight = halfWidth / aspect;
    } else {
      halfWidth = halfHeight * aspect;
    }
    camera.zoom = 1;
    camera.position.x = 0;
    camera.position.y = 0;
    camera.lookAt(0, 0, 0);
    camera.left = -halfWidth;
    camera.right = halfWidth;
    camera.top = halfHeight;
    camera.bottom = -halfHeight;
    camera.updateProjectionMatrix();
  }

  points.material.size = oldPointSize * (exportWidth / oldSize.x);
  renderer.render(scene, camera);

  const link = document.createElement("a");
  link.download = `chladni-${Date.now()}-${exportWidth}x${exportHeight}.png`;
  link.href = renderer.domElement.toDataURL("image/png");
  link.click();

  points.material.size = oldPointSize;
  renderer.setClearColor(oldClearColor, 1);
  renderer.setPixelRatio(oldPixelRatio);
  renderer.setSize(oldSize.x, oldSize.y, false);

  if (!useViewport) {
    CONFIG.cameraControl.x = oldCameraControl.x;
    CONFIG.cameraControl.y = oldCameraControl.y;
    CONFIG.cameraControl.z = oldCameraControl.z;
    applyCameraFromControl(CONFIG.cameraControl);
  }
  onWindowResize();
}

function randomizeWaves() {
  const keys = WAVE_TYPE_KEYS;
  CONFIG.waveTypeA = keys[Math.floor(Math.random() * keys.length)];
  CONFIG.waveTypeB = keys[Math.floor(Math.random() * keys.length)];

  rebuildField();
  if (refreshUI) refreshUI();
}

function randomizeAll() {
  // Fixed: your original call used parameters that randomizeModes() doesn't accept.
  randomizeModes();
  randomizeWaves();
  if (refreshUI) refreshUI();
}
