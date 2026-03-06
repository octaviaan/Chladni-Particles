import * as THREE from "three";
import { Pane } from "tweakpane";
import * as EssentialsPlugin from "tweakpane-plugin-essentials";
import * as InfodumpPlugin from "tweakpane-plugin-infodump";

// --- Configuration ---
const CONFIG = {
  particleCount: 500000,
  gridSize: 256,
  settleStrength: 5.0,
  jitter: 0.07,
  drag: 0.85,
  speedLimit: 2.0,
  viewScale: 400,
  color: "#4f4030",
  backgroundColor: "#1C1A17",
  particleSize: 2.0,
  particleOpacity: 0.7,
  cameraControl: { x: 0, y: 0, z: 1.5 },
  cameraDefault: { x: 0, y: 0, z: 1.5 },
  exportBaseSize: 2048,
  videoLoopSeconds: 5,
  videoFps: 12,
  videoWidth: 1280,
  videoHeight: 720,
  gifScale: 0.7,
  gifQuality: 50,
  formAnimate: false,
  formAnimSpeed: 2.0,

  modeCount: 6,
  mRange: { min: 2, max: 6 },
  nRange: { min: 4, max: 10 },

  rectAspect: 1.78,
  fieldScale: 0.9,

  waveTypeA: "Cartesian",
  waveTypeB: "Gyroid",
  waveMix: 0.5,
  phaseJitter: 0.69,

  integerModes: true,

  imageMode: false,
  imageBlend: 1.0,
  imageStyle: "Fill",
  imageThreshold: 0.7,
  imageMaskMode: "Auto",
  imageInvert: false,
  imageScale: 0.7,
  iconOffsetX: 0,
  logoParticleCount: 50000,
  logoParticleSize: 3.0,
  logoParticleOpacity: 0.9,
  logoColor: "#ffffff",
  textValue: "Gitcoin",
  textAlign: "Center",
  textFontSize: 180,
  textLineWidth: 0.9,
  textOffsetX: 0,
  textParticleCount: 100000,
  textParticleSize: 4.0,
  textParticleOpacity: 0.85,
  textColor: "#ffffff",
};

// --- Math Helpers ---
const smoothstep = (t) => t * t * (3 - 2 * t);

const hash2 = (x, y) => {
  const s = Math.sin(x * 127.1 + y * 311.7) * 43758.5453;
  return s - Math.floor(s);
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
  TriSineX: (cx, cy, mode) => {
    const rx = cx * mode.cos - cy * mode.sin;
    const ry = cx * mode.sin + cy * mode.cos;
    const a = mode.m * Math.PI;
    const b = mode.n * Math.PI;
    const c = ((mode.m + mode.n) * 0.5 + 0.5) * Math.PI;
    return (
      Math.sin(a * rx + mode.px) +
      Math.sin(b * ry + mode.py) +
      Math.sin(c * rx + mode.px * 0.5 + mode.py * 0.5)
    );
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
const ACCENT_COLOR = new THREE.Color("#02E2AC");
const ACCENT_COLOR_RATIO = 0.05;
const IMAGE_ATTRACT_STRENGTH = 0.03;
const TEXT_FONT_FAMILY = '"BDO Grotesk", sans-serif';

// --- Global Variables ---
let scene, camera, renderer, geometry, points, logoGeometry, logoPoints;
let textGeometry, textPoints;
let positions, velocities, colors, logoPositions, logoVelocities;
let textPositions, textVelocities;
let refreshUI = null;
let suppressUIEvents = false;
let pane;
let paneSecondary = null;
let isPanning = false;
let lastPointer = { x: 0, y: 0 };
let imageFileInput = null;
let logoImageSource = null;
let textImageSource = null;
let imagePointCloud = null;
let imageTargets = null;
let textPointCloud = null;
let textTargets = null;
let gifRecorder = null;
let gifRecordingTimer = null;
let gifFrameTimer = null;
let gifLibraryPromise = null;
let gifIsRendering = false;
let webmRecorder = null;
let webmRecordingTimer = null;
let webmRecordingChunks = [];
let captureStatusEl = null;
let captureStatusTimer = null;

// Field Data
let energy;
let gradX;
let gradY;

// Modes
let modes = [];

// Formation animation
let formAnimTime = 0;
let lastAnimTime = null;
let formProgress = null; // null = inactive, 0-1 = formation progress

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
  renderer.setClearColor(CONFIG.backgroundColor, 1);
  document.body.appendChild(renderer.domElement);

  buildParticles();
  buildLogoParticles();
  buildTextParticles();

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

  const logoMaterial = new THREE.PointsMaterial({
    color: CONFIG.logoColor,
    size: CONFIG.logoParticleSize,
    map: sprite,
    alphaTest: 0.0,
    transparent: true,
    opacity: CONFIG.logoParticleOpacity,
    blending: THREE.AdditiveBlending,
    depthWrite: false,
    depthTest: false,
    sizeAttenuation: false,
  });

  logoPoints = new THREE.Points(logoGeometry, logoMaterial);
  logoPoints.visible = false;
  logoPoints.renderOrder = 10;
  scene.add(logoPoints);

  const textMaterial = new THREE.PointsMaterial({
    color: CONFIG.textColor,
    size: CONFIG.textParticleSize,
    map: sprite,
    alphaTest: 0.0,
    transparent: true,
    opacity: CONFIG.textParticleOpacity,
    blending: THREE.AdditiveBlending,
    depthWrite: false,
    depthTest: false,
    sizeAttenuation: false,
  });

  textPoints = new THREE.Points(textGeometry, textMaterial);
  textPoints.visible = false;
  textPoints.renderOrder = 11;
  scene.add(textPoints);

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

  const range = CONFIG.viewScale;
  const rectRy = Math.max(1, Math.round(range / CONFIG.rectAspect));

  for (let i = 0; i < CONFIG.particleCount; i++) {
    const x = (Math.random() - 0.5) * 2 * range;
    const y = (Math.random() - 0.5) * 2 * rectRy;
    positions[i * 3 + 0] = x;
    positions[i * 3 + 1] = y;
    positions[i * 3 + 2] = 0;

    velocities[i * 2 + 0] = 0;
    velocities[i * 2 + 1] = 0;

    const colorRoll = Math.random();
    const particleColor =
      colorRoll < ACCENT_COLOR_RATIO ? ACCENT_COLOR : baseColor;
    colors[i * 3 + 0] = particleColor.r;
    colors[i * 3 + 1] = particleColor.g;
    colors[i * 3 + 2] = particleColor.b;
  }

  geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
  geometry.setAttribute("color", new THREE.BufferAttribute(colors, 3));
}

function buildLogoParticles() {
  logoGeometry = new THREE.BufferGeometry();

  const logoCount = Math.max(1, Math.round(CONFIG.logoParticleCount));
  logoPositions = new Float32Array(logoCount * 3);
  logoVelocities = new Float32Array(logoCount * 2);

  const range = CONFIG.viewScale;
  const rectRy = Math.max(1, Math.round(range / CONFIG.rectAspect));

  for (let i = 0; i < logoCount; i++) {
    const x = (Math.random() - 0.5) * 2 * range;
    const y = (Math.random() - 0.5) * 2 * rectRy;
    logoPositions[i * 3 + 0] = x;
    logoPositions[i * 3 + 1] = y;
    logoPositions[i * 3 + 2] = 0;

    logoVelocities[i * 2 + 0] = 0;
    logoVelocities[i * 2 + 1] = 0;
  }

  logoGeometry.setAttribute(
    "position",
    new THREE.BufferAttribute(logoPositions, 3),
  );
}

function buildTextParticles() {
  textGeometry = new THREE.BufferGeometry();

  const textCount = Math.max(1, Math.round(CONFIG.textParticleCount));
  textPositions = new Float32Array(textCount * 3);
  textVelocities = new Float32Array(textCount * 2);

  const range = CONFIG.viewScale;
  const rectRy = Math.max(1, Math.round(range / CONFIG.rectAspect));

  for (let i = 0; i < textCount; i++) {
    const x = (Math.random() - 0.5) * 2 * range;
    const y = (Math.random() - 0.5) * 2 * rectRy;
    textPositions[i * 3 + 0] = x;
    textPositions[i * 3 + 1] = y;
    textPositions[i * 3 + 2] = 0;

    textVelocities[i * 2 + 0] = 0;
    textVelocities[i * 2 + 1] = 0;
  }

  textGeometry.setAttribute(
    "position",
    new THREE.BufferAttribute(textPositions, 3),
  );
}

function rebuildParticles() {
  const oldGeometry = points.geometry;
  buildParticles();
  points.geometry = geometry;
  oldGeometry.dispose();
  applyParticleColor(CONFIG.color);
  resetLogoParticles();
  resetTextParticles();
}

function resetParticles() {
  if (!positions || !velocities || !geometry?.attributes?.position) return;

  const count = Math.min(
    CONFIG.particleCount,
    Math.floor(positions.length / 3),
    Math.floor(velocities.length / 2),
  );

  const range = CONFIG.viewScale;
  const rectRy = Math.max(1, Math.round(range / CONFIG.rectAspect));

  for (let i = 0; i < count; i++) {
    const i3 = i * 3;
    const i2 = i * 2;
    const x = (Math.random() - 0.5) * 2 * range;
    const y = (Math.random() - 0.5) * 2 * rectRy;
    positions[i3] = x;
    positions[i3 + 1] = y;
    positions[i3 + 2] = 0;
    velocities[i2] = 0;
    velocities[i2 + 1] = 0;
  }

  geometry.attributes.position.needsUpdate = true;
  resetLogoParticles();
  resetTextParticles();
}

function resetLogoParticles() {
  if (
    !logoPositions ||
    !logoVelocities ||
    !logoGeometry?.attributes?.position
  ) {
    return;
  }

  const hasTargets =
    imageTargets &&
    imageTargets.length === logoVelocities.length &&
    imageTargets.length > 0;
  if (hasTargets) {
    snapParticlesToImage();
    return;
  }

  const count = Math.min(
    Math.floor(logoPositions.length / 3),
    Math.floor(logoVelocities.length / 2),
  );
  const range = CONFIG.viewScale;
  const rectRy = Math.max(1, Math.round(range / CONFIG.rectAspect));

  for (let i = 0; i < count; i++) {
    const i3 = i * 3;
    const i2 = i * 2;
    logoPositions[i3] = (Math.random() - 0.5) * 2 * range;
    logoPositions[i3 + 1] = (Math.random() - 0.5) * 2 * rectRy;
    logoPositions[i3 + 2] = 0;
    logoVelocities[i2] = 0;
    logoVelocities[i2 + 1] = 0;
  }

  logoGeometry.attributes.position.needsUpdate = true;
  syncLogoVisibility();
}

function resetTextParticles() {
  if (
    !textPositions ||
    !textVelocities ||
    !textGeometry?.attributes?.position
  ) {
    return;
  }

  const hasTargets =
    textTargets &&
    textTargets.length === textVelocities.length &&
    textTargets.length > 0;
  if (hasTargets) {
    snapParticlesToText();
    return;
  }

  const count = Math.min(
    Math.floor(textPositions.length / 3),
    Math.floor(textVelocities.length / 2),
  );
  const range = CONFIG.viewScale;
  const rectRy = Math.max(1, Math.round(range / CONFIG.rectAspect));

  for (let i = 0; i < count; i++) {
    const i3 = i * 3;
    const i2 = i * 2;
    textPositions[i3] = (Math.random() - 0.5) * 2 * range;
    textPositions[i3 + 1] = (Math.random() - 0.5) * 2 * rectRy;
    textPositions[i3 + 2] = 0;
    textVelocities[i2] = 0;
    textVelocities[i2 + 1] = 0;
  }

  textGeometry.attributes.position.needsUpdate = true;
  syncLogoVisibility();
}

function resetAllParticles() {
  resetParticles();
}

function rebuildLogoParticles() {
  if (!logoPoints) return;
  const oldGeometry = logoPoints.geometry;
  buildLogoParticles();
  logoPoints.geometry = logoGeometry;
  if (oldGeometry) oldGeometry.dispose();

  if (imagePointCloud) {
    rebuildImageTargets();
    if (CONFIG.imageMode) snapParticlesToImage();
  } else {
    imageTargets = null;
  }
  syncLogoVisibility();
}

function rebuildTextParticles() {
  if (!textPoints) return;
  const oldGeometry = textPoints.geometry;
  buildTextParticles();
  textPoints.geometry = textGeometry;
  if (oldGeometry) oldGeometry.dispose();

  if (textPointCloud) {
    rebuildTextTargets();
    if (CONFIG.imageMode) snapParticlesToText();
  } else {
    textTargets = null;
  }
  syncLogoVisibility();
}

function applyParticleColor(hex) {
  if (!colors || !geometry?.attributes?.color) return;
  const nextColor = new THREE.Color(hex);
  const count = CONFIG.particleCount;
  for (let i = 0; i < count; i++) {
    const i3 = i * 3;
    const colorRoll = Math.random();
    const particleColor =
      colorRoll < ACCENT_COLOR_RATIO ? ACCENT_COLOR : nextColor;
    colors[i3] = particleColor.r;
    colors[i3 + 1] = particleColor.g;
    colors[i3 + 2] = particleColor.b;
  }
  geometry.attributes.color.needsUpdate = true;
}

function applyLogoColor(hex) {
  if (!logoPoints?.material) return;
  logoPoints.material.color.set(hex);
}

function applyTextColor(hex) {
  if (!textPoints?.material) return;
  textPoints.material.color.set(hex);
}

function syncLogoVisibility() {
  if (!logoPoints || !textPoints) return;
  const logoCount = logoVelocities ? logoVelocities.length / 2 : 0;
  const textCount = textVelocities ? textVelocities.length / 2 : 0;
  const hasLogoTargets =
    imageTargets && logoCount > 0 && imageTargets.length === logoCount * 2;
  const hasTextTargets =
    textTargets && textCount > 0 && textTargets.length === textCount * 2;
  logoPoints.visible =
    !!CONFIG.imageMode && !!hasLogoTargets && CONFIG.logoParticleOpacity > 0;
  textPoints.visible =
    !!CONFIG.imageMode && !!hasTextTargets && CONFIG.textParticleOpacity > 0;
}

function refreshImageModeFromTargets() {
  CONFIG.imageMode = !!imageTargets || !!textTargets;
  syncLogoVisibility();
}

function rebuildField() {
  const G = CONFIG.gridSize;
  const funcA = WAVE_FUNCTIONS[CONFIG.waveTypeA] || WAVE_FUNCTIONS["Cartesian"];
  const funcB = WAVE_FUNCTIONS[CONFIG.waveTypeB] || WAVE_FUNCTIONS["Cartesian"];
  const baseBias = CONFIG.waveMix - 0.5;
  const fieldScale = Math.max(0.05, CONFIG.fieldScale || 1);
  const maxR = Math.SQRT1_2 * fieldScale;

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
      const mixed = rn + baseBias;
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

  resetLogoParticles();
  resetTextParticles();
}

function scatterParticles() {
  if (!positions || !velocities) return;
  const range = CONFIG.viewScale;
  const rectRy = Math.max(1, Math.round(range / CONFIG.rectAspect));
  const count = Math.min(
    CONFIG.particleCount,
    Math.floor(positions.length / 3),
  );
  for (let i = 0; i < count; i++) {
    positions[i * 3] = (Math.random() - 0.5) * 2 * range;
    positions[i * 3 + 1] = (Math.random() - 0.5) * 2 * rectRy;
    positions[i * 3 + 2] = 0;
    velocities[i * 2] = 0;
    velocities[i * 2 + 1] = 0;
  }
  if (geometry?.attributes?.position)
    geometry.attributes.position.needsUpdate = true;
}

function updateParticles() {
  if (!positions || !velocities || !geometry?.attributes?.position) return;

  const count = Math.min(
    CONFIG.particleCount,
    Math.floor(positions.length / 3),
    Math.floor(velocities.length / 2),
  );

  const targets =
    CONFIG.imageMode && imageTargets && imageTargets.length === count * 2
      ? imageTargets
      : null;

  stepSimulation(positions, velocities, count, targets);
  geometry.attributes.position.needsUpdate = true;
}

function stepSimulation(pos, vel, count, targets) {
  const G = CONFIG.gridSize;
  const range = CONFIG.viewScale;
  const fullRange = range * 2;
  const settle = CONFIG.settleStrength;
  const jitter = CONFIG.jitter;
  const drag = CONFIG.drag;
  const limit = CONFIG.speedLimit;
  const limitSq = limit * limit;

  const gridScale = (G - 1) / fullRange;
  const rectRy = Math.max(1, Math.round(range / CONFIG.rectAspect));

  let useTargets = false;
  let fieldMult = 1.0;
  let targetMult = 0.0;
  let jitterMult = 1.0;

  if (targets) {
    useTargets = true;
    const blend = clamp(CONFIG.imageBlend, 0, 1);
    fieldMult = 1.0 - blend;
    targetMult = blend * IMAGE_ATTRACT_STRENGTH;
    // When particles are driven to image/text targets, reduce random noise
    // so they can settle instead of endlessly re-attracting.
    jitterMult = 1.0 - blend;
  }

  if (formProgress !== null) {
    fieldMult *= formProgress;
    jitterMult *= 1.0 - formProgress * 0.95;
  }

  for (let i = 0; i < count; i++) {
    const i3 = i * 3;
    const i2 = i * 2;

    let x = pos[i3];
    let y = pos[i3 + 1];
    let vx = vel[i2];
    let vy = vel[i2 + 1];

    // Keep particles within the simulation bounds without random respawn.
    if (x < -range || x > range || y < -rectRy || y > rectRy) {
      x = clamp(x, -range, range);
      y = clamp(y, -rectRy, rectRy);
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

    if (fieldMult > 0.001) {
      vx -= gxVal * settle * fieldMult;
      vy -= gyVal * settle * fieldMult;
    }

    if (useTargets && targetMult > 0.001) {
      vx += (targets[i2] - x) * targetMult;
      vy += (targets[i2 + 1] - y) * targetMult;
    }

    vx += (Math.random() - 0.5) * jitter * jitterMult;
    vy += (Math.random() - 0.5) * jitter * jitterMult;

    vx *= drag;
    vy *= drag;

    const speedSq = vx * vx + vy * vy;

    if (speedSq > limitSq) {
      const scale = limit / Math.sqrt(speedSq);
      vx *= scale;
      vy *= scale;
    }

    x += vx;
    y += vy;

    pos[i3] = x;
    pos[i3 + 1] = y;

    vel[i2] = vx;
    vel[i2 + 1] = vy;
  }
}

function updateLogoParticles() {
  if (
    !logoPositions ||
    !logoVelocities ||
    !logoGeometry?.attributes?.position
  ) {
    return;
  }

  const count = Math.min(
    Math.floor(logoPositions.length / 3),
    Math.floor(logoVelocities.length / 2),
  );
  const hasTargets =
    CONFIG.imageMode &&
    CONFIG.logoParticleOpacity > 0 &&
    imageTargets &&
    imageTargets.length === count * 2;
  syncLogoVisibility();
  if (!hasTargets) return;

  stepSimulation(logoPositions, logoVelocities, count, imageTargets);
  logoGeometry.attributes.position.needsUpdate = true;
}

function updateTextParticles() {
  if (
    !textPositions ||
    !textVelocities ||
    !textGeometry?.attributes?.position
  ) {
    return;
  }

  const count = Math.min(
    Math.floor(textPositions.length / 3),
    Math.floor(textVelocities.length / 2),
  );
  const hasTargets =
    CONFIG.imageMode &&
    CONFIG.textParticleOpacity > 0 &&
    textTargets &&
    textTargets.length === count * 2;
  syncLogoVisibility();
  if (!hasTargets) return;

  stepSimulation(textPositions, textVelocities, count, textTargets);
  textGeometry.attributes.position.needsUpdate = true;
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

function randomPhaseOffset() {
  const phaseRange = Math.max(0, CONFIG.phaseJitter) * Math.PI * 2;
  return (Math.random() - 0.5) * phaseRange;
}

function randomizeModes(rebuild = true) {
  const maxModes = 12;
  CONFIG.modeCount = Math.floor(1 + Math.random() * maxModes);

  const rangeFloor = 1.0;
  const rangeCeil = 10.0;

  const biased = (power) => Math.pow(Math.random(), power);

  // Bias ranges toward low-to-mid frequencies while still allowing variation.
  const mMin = rangeFloor + biased(1.4) * 4.5; // ~1..5.5
  const mSpan = 1.2 + biased(1.8) * 3.3; // ~1.2..4.5
  CONFIG.mRange.min = Math.round(mMin * 10) / 10;
  CONFIG.mRange.max =
    Math.round(Math.min(rangeCeil, CONFIG.mRange.min + mSpan) * 10) / 10;

  const nMin = rangeFloor + biased(1.5) * 5.2; // ~1..6.2
  const nSpan = 1.4 + biased(1.9) * 3.6; // ~1.4..5.0
  CONFIG.nRange.min = Math.round(nMin * 10) / 10;
  CONFIG.nRange.max =
    Math.round(Math.min(rangeCeil, CONFIG.nRange.min + nSpan) * 10) / 10;

  CONFIG.integerModes = Math.random() < 0.5;

  initModes();

  if (rebuild) {
    rebuildField();
    if (refreshUI) {
      suppressUIEvents = true;
      refreshUI();
      suppressUIEvents = false;
    }
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
    const px = randomPhaseOffset();
    const py = randomPhaseOffset();
    const angle = (Math.random() - 0.5) * 2 * Math.PI;
    modes.push({
      m: quantizeModeValue(mRaw),
      n: quantizeModeValue(nRaw),
      a: Math.exp(-i * 0.35),
      px: px,
      py: py,
      cos: Math.cos(angle),
      sin: Math.sin(angle),
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
  const now = performance.now();

  if (CONFIG.formAnimate) {
    const dt = lastAnimTime !== null ? (now - lastAnimTime) / 1000 : 0;
    const loopSec = Math.max(0.1, CONFIG.videoLoopSeconds);
    formAnimTime = Math.min(formAnimTime + dt * CONFIG.formAnimSpeed, loopSec);
    const t = formAnimTime / loopSec; // 0..1
    formProgress = smoothstep(clamp(t * 1.4, 0, 1));
  } else {
    formProgress = null;
  }

  lastAnimTime = now;
  updateParticles();
  updateLogoParticles();
  updateTextParticles();
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

function hasImageSources() {
  return !!logoImageSource || !!textImageSource;
}

function ensureImageFileInput() {
  if (imageFileInput) return imageFileInput;

  imageFileInput = document.createElement("input");
  imageFileInput.type = "file";
  imageFileInput.accept = "image/*";
  imageFileInput.style.display = "none";
  imageFileInput.addEventListener("change", onImageFileSelected);
  document.body.appendChild(imageFileInput);

  return imageFileInput;
}

function promptImageUpload() {
  const input = ensureImageFileInput();
  input.click();
}

function onImageFileSelected(event) {
  const input = event.target;
  const file = input.files && input.files[0];
  if (!file) return;

  const imageUrl = URL.createObjectURL(file);
  const img = new Image();
  img.onload = () => {
    updateImageSource(img, file.name);
    const hasPoints = rebuildImagePointCloud();
    if (hasPoints) {
      CONFIG.imageMode = true;
      rebuildImageTargets();
      snapParticlesToImage();
      syncLogoVisibility();
      if (refreshUI) refreshUI();
    } else {
      imageTargets = null;
      if (!textTargets) CONFIG.imageMode = false;
      syncLogoVisibility();
    }
    URL.revokeObjectURL(imageUrl);
    input.value = "";
  };
  img.onerror = () => {
    URL.revokeObjectURL(imageUrl);
    input.value = "";
  };
  img.src = imageUrl;
}

function updateImageSource(img, name) {
  const maxDim = 512;
  const longest = Math.max(img.width, img.height);
  const scale = longest > maxDim ? maxDim / longest : 1;
  const width = Math.max(1, Math.round(img.width * scale));
  const height = Math.max(1, Math.round(img.height * scale));

  const canvas = document.createElement("canvas");
  canvas.width = width;
  canvas.height = height;
  const ctx = canvas.getContext("2d", { willReadFrequently: true });
  if (!ctx) {
    logoImageSource = null;
    imagePointCloud = null;
    imageTargets = null;
    return;
  }
  ctx.drawImage(img, 0, 0, width, height);
  const pixels = ctx.getImageData(0, 0, width, height).data;

  logoImageSource = {
    name,
    width,
    height,
    pixels: new Uint8ClampedArray(pixels),
  };
}

function updateTextSource() {
  const rawText = typeof CONFIG.textValue === "string" ? CONFIG.textValue : "";
  const text = rawText.trim();
  if (!text) {
    textImageSource = null;
    return false;
  }

  const width = 1024;
  const height = 512;
  const canvas = document.createElement("canvas");
  canvas.width = width;
  canvas.height = height;
  const ctx = canvas.getContext("2d", { willReadFrequently: true });
  if (!ctx) return false;

  const weight = 700;
  const maxLineWidthRatio = clamp(CONFIG.textLineWidth, 0.4, 0.95);
  const maxTextWidth = width * maxLineWidthRatio;
  const maxTextHeight = height * 0.72;
  let fontSize = clamp(Math.round(CONFIG.textFontSize), 16, 360);
  const lineHeightFactor = 1.15;
  const align = CONFIG.textAlign || "Center";
  const blockLeft = (width - maxTextWidth) * 0.5;
  const blockRight = width - blockLeft;
  const drawX =
    align === "Left" ? blockLeft : align === "Right" ? blockRight : width * 0.5;

  ctx.fillStyle = "#000000";
  ctx.textAlign =
    align === "Left" ? "left" : align === "Right" ? "right" : "center";
  ctx.textBaseline = "middle";

  const wrapText = (value) => {
    const words = value.split(/\s+/).filter(Boolean);
    if (words.length === 0) return [];

    const lines = [];
    let currentLine = "";

    for (let i = 0; i < words.length; i++) {
      const word = words[i];
      const candidate = currentLine ? `${currentLine} ${word}` : word;
      if (ctx.measureText(candidate).width <= maxTextWidth) {
        currentLine = candidate;
        continue;
      }

      if (currentLine) lines.push(currentLine);

      if (ctx.measureText(word).width <= maxTextWidth) {
        currentLine = word;
        continue;
      }

      let chunk = "";
      for (const char of word) {
        const next = chunk + char;
        if (ctx.measureText(next).width > maxTextWidth && chunk) {
          lines.push(chunk);
          chunk = char;
        } else {
          chunk = next;
        }
      }
      currentLine = chunk;
    }

    if (currentLine) lines.push(currentLine);
    return lines;
  };

  let lines = [];
  while (fontSize > 14) {
    ctx.font = `${weight} ${fontSize}px ${TEXT_FONT_FAMILY}`;
    lines = wrapText(text);
    const textHeight = lines.length * fontSize * lineHeightFactor;
    if (lines.length > 0 && textHeight <= maxTextHeight) break;
    fontSize -= 2;
  }

  ctx.clearRect(0, 0, width, height);
  ctx.font = `${weight} ${fontSize}px ${TEXT_FONT_FAMILY}`;
  lines = wrapText(text);

  const lineHeight = fontSize * lineHeightFactor;
  const startY = height * 0.5 - ((lines.length - 1) * lineHeight) / 2;
  for (let i = 0; i < lines.length; i++) {
    ctx.fillText(lines[i], drawX, startY + i * lineHeight);
  }

  const pixels = ctx.getImageData(0, 0, width, height).data;
  textImageSource = {
    name: `Text: ${text}`,
    width,
    height,
    pixels: new Uint8ClampedArray(pixels),
  };
  return true;
}

function applyTextSource() {
  updateTextSource();

  const hasPoints = rebuildTextPointCloud();
  if (hasPoints) {
    CONFIG.imageMode = true;
    rebuildTextTargets();
    snapParticlesToText();
    syncLogoVisibility();
    if (refreshUI) refreshUI();
  } else {
    textPointCloud = null;
    textTargets = null;
    if (!imageTargets) CONFIG.imageMode = false;
    syncLogoVisibility();
  }
}

function createPointCloudFromSource(source, offsetX) {
  if (!source) return null;

  const threshold = clamp(CONFIG.imageThreshold, 0, 1);
  const alphaCutoff = 0.02;
  const points = [];
  const { width, height, pixels } = source;
  const pixelCount = width * height;
  const baseMask = new Uint8Array(pixelCount);
  let transparentPixels = 0;
  let opaquePixels = 0;
  let minBrightness = 1;
  let maxBrightness = 0;

  for (let i = 0; i < pixelCount; i++) {
    const idx = i * 4;
    const alpha = pixels[idx + 3] / 255;
    if (alpha <= alphaCutoff) {
      transparentPixels++;
      continue;
    }
    const brightness = (pixels[idx] + pixels[idx + 1] + pixels[idx + 2]) / 765;
    if (brightness < minBrightness) minBrightness = brightness;
    if (brightness > maxBrightness) maxBrightness = brightness;
    opaquePixels++;
  }

  const hasMeaningfulTransparency =
    transparentPixels > 0 && transparentPixels / pixelCount > 0.01;
  const luminanceSpread = opaquePixels > 0 ? maxBrightness - minBrightness : 0;
  // Auto mode favors alpha only when transparency defines the silhouette and
  // opaque pixels do not have strong brightness contrast.
  const autoUseAlphaMask = hasMeaningfulTransparency && luminanceSpread < 0.2;
  const useAlphaMask =
    CONFIG.imageMaskMode === "Alpha"
      ? true
      : CONFIG.imageMaskMode === "Luminance"
        ? false
        : autoUseAlphaMask;
  const alphaThreshold = 1 - threshold;

  for (let i = 0; i < pixelCount; i++) {
    const idx = i * 4;
    const alpha = pixels[idx + 3] / 255;
    if (alpha <= alphaCutoff) continue;

    let selected = false;
    if (useAlphaMask) {
      selected = CONFIG.imageInvert
        ? alpha <= alphaThreshold
        : alpha >= alphaThreshold;
    } else {
      const brightness =
        (pixels[idx] + pixels[idx + 1] + pixels[idx + 2]) / 765;
      selected = CONFIG.imageInvert
        ? brightness >= threshold
        : brightness <= threshold;
    }

    if (selected) baseMask[i] = 1;
  }

  for (let y = 0; y < height; y++) {
    for (let x = 0; x < width; x++) {
      const i = y * width + x;
      if (!baseMask[i]) continue;

      if (CONFIG.imageStyle === "Outline") {
        const left = x > 0 ? baseMask[i - 1] : 0;
        const right = x + 1 < width ? baseMask[i + 1] : 0;
        const up = y > 0 ? baseMask[i - width] : 0;
        const down = y + 1 < height ? baseMask[i + width] : 0;
        if (left && right && up && down) continue;
      }

      const nx = (x + 0.5) / width - 0.5 + clamp(offsetX, -1.5, 1.5);
      const ny = 0.5 - (y + 0.5) / height;
      points.push(nx, ny);
    }
  }

  if (points.length === 0) return null;

  return {
    name: source.name,
    width,
    height,
    points: new Float32Array(points),
  };
}

function rebuildImagePointCloud() {
  imagePointCloud = createPointCloudFromSource(
    logoImageSource,
    CONFIG.iconOffsetX,
  );
  if (!imagePointCloud) {
    imageTargets = null;
    syncLogoVisibility();
    return false;
  }
  return true;
}

function rebuildTextPointCloud() {
  textPointCloud = createPointCloudFromSource(
    textImageSource,
    CONFIG.textOffsetX,
  );
  if (!textPointCloud) {
    textTargets = null;
    syncLogoVisibility();
    return false;
  }
  return true;
}

function createTargetsFromPointCloud(pointCloud, targetCount) {
  if (!pointCloud) return null;

  const pointCount = pointCloud.points.length / 2;
  if (pointCount === 0) return null;

  const targets = new Float32Array(targetCount * 2);

  const range = CONFIG.viewScale;
  const rectRy = Math.max(1, Math.round(range / CONFIG.rectAspect));
  const logoScale = clamp(CONFIG.imageScale, 0.1, 1.5);
  const boundHalfWidth = range * 0.9 * logoScale;
  const boundHalfHeight = rectRy * 0.9 * logoScale;
  const boundsAspect = boundHalfWidth / boundHalfHeight;
  const sourceAspect = pointCloud.width / pointCloud.height;

  let halfWidth = boundHalfWidth;
  let halfHeight = boundHalfHeight;
  if (sourceAspect > boundsAspect) {
    halfHeight = halfWidth / sourceAspect;
  } else {
    halfWidth = halfHeight * sourceAspect;
  }

  const points = pointCloud.points;
  const pixelJitterX = 0.9 / pointCloud.width;
  const pixelJitterY = 0.9 / pointCloud.height;

  for (let i = 0; i < targetCount; i++) {
    const i2 = i * 2;
    const src = Math.floor(Math.random() * pointCount);
    const src2 = src * 2;
    const nx = clamp(
      points[src2] + (Math.random() - 0.5) * pixelJitterX,
      -1.5,
      1.5,
    );
    const ny = clamp(
      points[src2 + 1] + (Math.random() - 0.5) * pixelJitterY,
      -1.5,
      1.5,
    );
    targets[i2] = nx * (halfWidth * 2);
    targets[i2 + 1] = ny * (halfHeight * 2);
  }

  return targets;
}

function rebuildImageTargets() {
  if (!imagePointCloud) {
    imageTargets = null;
    syncLogoVisibility();
    return;
  }

  const targetCount = Math.max(1, Math.round(CONFIG.logoParticleCount));
  imageTargets = createTargetsFromPointCloud(imagePointCloud, targetCount);
  syncLogoVisibility();
}

function rebuildTextTargets() {
  if (!textPointCloud) {
    textTargets = null;
    syncLogoVisibility();
    return;
  }

  const targetCount = Math.max(1, Math.round(CONFIG.textParticleCount));
  textTargets = createTargetsFromPointCloud(textPointCloud, targetCount);
  syncLogoVisibility();
}

function clearImageData() {
  logoImageSource = null;
  imagePointCloud = null;
  imageTargets = null;
  if (!textTargets) CONFIG.imageMode = false;
  syncLogoVisibility();
  if (refreshUI) refreshUI();
}

function clearTextData() {
  textImageSource = null;
  textPointCloud = null;
  textTargets = null;
  if (!imageTargets) CONFIG.imageMode = false;
  syncLogoVisibility();
  if (refreshUI) refreshUI();
}

function snapParticlesToImage() {
  if (
    !imageTargets ||
    !logoPositions ||
    !logoVelocities ||
    imageTargets.length !== logoVelocities.length
  ) {
    return;
  }

  const count = logoVelocities.length / 2;
  for (let i = 0; i < count; i++) {
    const i2 = i * 2;
    const i3 = i * 3;
    logoPositions[i3] = imageTargets[i2];
    logoPositions[i3 + 1] = imageTargets[i2 + 1];
    logoPositions[i3 + 2] = 0;
    logoVelocities[i2] = 0;
    logoVelocities[i2 + 1] = 0;
  }
  logoGeometry.attributes.position.needsUpdate = true;
  syncLogoVisibility();
}

function snapParticlesToText() {
  if (
    !textTargets ||
    !textPositions ||
    !textVelocities ||
    textTargets.length !== textVelocities.length
  ) {
    return;
  }

  const count = textVelocities.length / 2;
  for (let i = 0; i < count; i++) {
    const i2 = i * 2;
    const i3 = i * 3;
    textPositions[i3] = textTargets[i2];
    textPositions[i3 + 1] = textTargets[i2 + 1];
    textPositions[i3 + 2] = 0;
    textVelocities[i2] = 0;
    textVelocities[i2 + 1] = 0;
  }
  textGeometry.attributes.position.needsUpdate = true;
  syncLogoVisibility();
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
  const key = event.key.toLowerCase();
  if (key === "r") randomizeAll();
  if (key === "s") saveImage({ useViewport: false });
  if (key === "c") saveImage({ useViewport: true });
  if (key === "v") recordVideoLoop();
  if (key === "h" && pane) {
    const nextHidden = !pane.hidden;
    pane.hidden = nextHidden;
    if (paneSecondary) paneSecondary.hidden = nextHidden;
  }
}

function setupGUI() {
  const panesWrapper = document.createElement("div");
  panesWrapper.className = "tp-columns";
  document.body.appendChild(panesWrapper);

  pane = new Pane({ title: "Controls", container: panesWrapper });
  paneSecondary = new Pane({ title: "Logo & Text", container: panesWrapper });

  pane.element.classList.add("tp-minimal");
  paneSecondary.element.classList.add("tp-minimal", "tp-minimal-secondary");
  pane.registerPlugin(EssentialsPlugin);
  pane.registerPlugin(InfodumpPlugin);
  paneSecondary.registerPlugin(EssentialsPlugin);
  paneSecondary.registerPlugin(InfodumpPlugin);
  const inputs = [];

  const legendStyle = document.createElement("style");
  legendStyle.textContent = `
    .tp-legend{
      padding: 8px;
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
    keyS: `
    <svg width="24" height="24" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
<rect x="2.5" y="2.5" width="19" height="19" rx="2.5" stroke="#B8B8B8"/>
<text x="12" y="16.5" font-family="sans-serif" font-size="13" font-weight="600" fill="#AFCEA1" text-anchor="middle">S</text>
</svg>
`,
    keyC: `
    <svg width="24" height="24" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
<rect x="2.5" y="2.5" width="19" height="19" rx="2.5" stroke="#B8B8B8"/>
<text x="12" y="16.5" font-family="sans-serif" font-size="13" font-weight="600" fill="#AFCEA1" text-anchor="middle">C</text>
</svg>
`,
    keyV: `
    <svg width="24" height="24" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
<rect x="2.5" y="2.5" width="19" height="19" rx="2.5" stroke="#B8B8B8"/>
<text x="12" y="16.5" font-family="sans-serif" font-size="13" font-weight="600" fill="#AFCEA1" text-anchor="middle">V</text>
</svg>
`,
    keyH: `
    <svg width="24" height="24" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
<rect x="2.5" y="2.5" width="19" height="19" rx="2.5" stroke="#B8B8B8"/>
<text x="12" y="16.5" font-family="sans-serif" font-size="13" font-weight="600" fill="#AFCEA1" text-anchor="middle">H</text>
</svg>
`,
  };

  const shortcutsFolder = pane.addFolder({
    title: "Shortcuts",
    expanded: false,
  });

  const legend = document.createElement("div");
  legend.className = "tp-legend";
  legend.innerHTML = `
    <div class="tp-legend-row"><span class="tp-legend-ico">${svgs.wheel}</span><span>Wheel: zoom</span></div>
    <div class="tp-legend-row"><span class="tp-legend-ico">${svgs.drag}</span><span>Drag: pan</span></div>
    <div class="tp-legend-row"><span class="tp-legend-ico">${svgs.dbl}</span><span>Double click: reset position</span></div>
    <div class="tp-legend-row"><span class="tp-legend-ico">${svgs.keyR}</span><span>R: randomize all</span></div>
    <div class="tp-legend-row"><span class="tp-legend-ico">${svgs.keyS}</span><span>S: save image</span></div>
    <div class="tp-legend-row"><span class="tp-legend-ico">${svgs.keyC}</span><span>C: save current view</span></div>
    <div class="tp-legend-row"><span class="tp-legend-ico">${svgs.keyV}</span><span>V: save GIF loop</span></div>
    <div class="tp-legend-row"><span class="tp-legend-ico">${svgs.keyH}</span><span>H: toggle controls</span></div>
  `;
  const container = shortcutsFolder.element.querySelector(".tp-fldv_c");
  if (container) {
    container.appendChild(legend);
  }

  const waveTypeOptions = WAVE_TYPE_KEYS.reduce((obj, key) => {
    obj[key] = key;
    return obj;
  }, {});

  const guideFolder = pane.addFolder({
    title: "Guide",
    expanded: false,
  });
  guideFolder.addBlade({
    view: "infodump",
    markdown: true,
    content: [
      "**Waves + Modes:** Pick `Wave A`, `Wave B`, `Mix`, then adjust `modeCount`/ranges. Use `Randomize All` (or `R`) to quickly explore.",
      "**Logo/Text mask:** Upload a logo or add text. `Mask` supports `Auto`, `Alpha`, `Luminance`. If targets vanish, lower `Threshold` or switch mask mode.",
      "**Capture:** `Image resolution` controls PNG width. Use `Save GIF` (`Duration`, `GIF fps`, `GIF scale`, `GIF compression`) or `Save WebM` for video.",
    ].join("\n\n"),
  });

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
      .addBinding(CONFIG, "particleOpacity", {
        min: 0.1,
        max: 1.0,
        step: 0.01,
        label: "Opacity",
      })
      .on("change", (ev) => {
        points.material.opacity = ev.value;
      }),
  );
  inputs.push(
    particlesFolder
      .addBinding(CONFIG, "color", { label: "Color" })
      .on("change", (ev) => applyParticleColor(ev.value)),
  );
  const motionFolder = pane.addFolder({ title: "Motion" });
  inputs.push(
    motionFolder
      .addBinding(CONFIG, "settleStrength", {
        min: 0.5,
        max: 10.0,
        step: 0.1,
      })
      .on("change", () => rebuildParticles()),
  );
  inputs.push(
    motionFolder
      .addBinding(CONFIG, "jitter", {
        min: 0.0,
        max: 0.3,
        step: 0.01,
      })
      .on("change", () => rebuildParticles()),
  );
  inputs.push(
    motionFolder
      .addBinding(CONFIG, "drag", { min: 0.7, max: 0.9, step: 0.01 })
      .on("change", () => rebuildParticles()),
  );
  inputs.push(
    motionFolder
      .addBinding(CONFIG, "speedLimit", {
        min: 0.2,
        max: 3.0,
        step: 0.1,
      })
      .on("change", () => rebuildParticles()),
  );

  const wavesFolder = pane.addFolder({ title: "Waves" });
  inputs.push(
    wavesFolder
      .addBinding(CONFIG, "waveTypeA", {
        options: waveTypeOptions,
        label: "Wave A",
      })
      .on("change", () => {
        rebuildField();
      }),
  );
  inputs.push(
    wavesFolder
      .addBinding(CONFIG, "waveTypeB", {
        options: waveTypeOptions,
        label: "Wave B",
      })
      .on("change", () => {
        rebuildField();
      }),
  );
  inputs.push(
    wavesFolder
      .addBinding(CONFIG, "waveMix", {
        min: 0.0,
        max: 1.0,
        step: 0.01,
        label: "Mix",
      })
      .on("change", () => {
        rebuildField();
      }),
  );
  inputs.push(
    wavesFolder
      .addBinding(CONFIG, "fieldScale", {
        min: 0.3,
        max: 2.0,
        step: 0.1,
        label: "Field scale",
      })
      .on("change", () => {
        rebuildField();
      }),
  );
  const modesFolder = pane.addFolder({ title: "Modes" });
  inputs.push(
    modesFolder
      .addBinding(CONFIG, "modeCount", { min: 1, max: 10, step: 1 })
      .on("change", () => {
        if (!suppressUIEvents) {
          rebuildModesFromConfig();
          rebuildParticles();
        }
      }),
  );
  inputs.push(
    modesFolder
      .addBinding(CONFIG, "mRange", { min: 1, max: 10, step: 0.1 })
      .on("change", () => {
        if (!suppressUIEvents) {
          rebuildModesFromConfig();
          rebuildParticles();
        }
      }),
  );
  inputs.push(
    modesFolder
      .addBinding(CONFIG, "nRange", { min: 1, max: 10, step: 0.1 })
      .on("change", () => {
        if (!suppressUIEvents) {
          rebuildModesFromConfig();
          rebuildParticles();
        }
      }),
  );
  modesFolder.addButton({ title: "Randomize All" }).on("click", () => {
    randomizeAll();
  });

  const imageFolder = paneSecondary.addFolder({ title: "Logo" });
  imageFolder.addButton({ title: "Upload logo" }).on("click", () => {
    promptImageUpload();
  });
  imageFolder.addButton({ title: "Clear logo" }).on("click", () => {
    clearImageData();
  });
  inputs.push(
    imageFolder
      .addBinding(CONFIG, "imageScale", {
        min: 0.1,
        max: 1.5,
        step: 0.01,
        label: "Scale",
      })
      .on("change", () => {
        if (imagePointCloud) {
          rebuildImageTargets();
          if (CONFIG.imageMode) snapParticlesToImage();
        }
        if (textPointCloud) {
          rebuildTextTargets();
          if (CONFIG.imageMode) snapParticlesToText();
        }
        refreshImageModeFromTargets();
      }),
  );
  inputs.push(
    imageFolder
      .addBinding(CONFIG, "iconOffsetX", {
        min: -1.5,
        max: 1.5,
        step: 0.01,
        label: "X offset",
      })
      .on("change", () => {
        if (!logoImageSource) return;
        if (rebuildImagePointCloud()) {
          rebuildImageTargets();
          if (CONFIG.imageMode) snapParticlesToImage();
        }
        refreshImageModeFromTargets();
      }),
  );
  inputs.push(
    imageFolder
      .addBinding(CONFIG, "logoColor", {
        label: "Logo color",
      })
      .on("change", (ev) => applyLogoColor(ev.value)),
  );
  inputs.push(
    imageFolder
      .addBinding(CONFIG, "logoParticleCount", {
        min: 1000,
        max: 200000,
        step: 1000,
        label: "Logo count",
      })
      .on("change", () => {
        rebuildLogoParticles();
      }),
  );
  inputs.push(
    imageFolder
      .addBinding(CONFIG, "logoParticleSize", {
        min: 1.0,
        max: 6.0,
        step: 0.1,
        label: "Particles size",
      })
      .on("change", (ev) => {
        if (logoPoints?.material) logoPoints.material.size = ev.value;
      }),
  );
  inputs.push(
    imageFolder
      .addBinding(CONFIG, "logoParticleOpacity", {
        min: 0.0,
        max: 1.0,
        step: 0.01,
        label: "Opacity",
      })
      .on("change", (ev) => {
        if (logoPoints?.material) logoPoints.material.opacity = ev.value;
        syncLogoVisibility();
      }),
  );
  inputs.push(
    imageFolder
      .addBinding(CONFIG, "imageThreshold", {
        min: 0.0,
        max: 1.0,
        step: 0.01,
        label: "Threshold",
      })
      .on("change", () => {
        if (!hasImageSources()) return;
        if (rebuildImagePointCloud()) {
          rebuildImageTargets();
          if (CONFIG.imageMode) snapParticlesToImage();
        }
        if (rebuildTextPointCloud()) {
          rebuildTextTargets();
          if (CONFIG.imageMode) snapParticlesToText();
        }
        refreshImageModeFromTargets();
      }),
  );
  inputs.push(
    imageFolder
      .addBinding(CONFIG, "imageMaskMode", {
        options: {
          Auto: "Auto",
          Alpha: "Alpha",
          Luminance: "Luminance",
        },
        label: "Mask",
      })
      .on("change", () => {
        if (!hasImageSources()) return;
        if (rebuildImagePointCloud()) {
          rebuildImageTargets();
          if (CONFIG.imageMode) snapParticlesToImage();
        }
        if (rebuildTextPointCloud()) {
          rebuildTextTargets();
          if (CONFIG.imageMode) snapParticlesToText();
        }
        refreshImageModeFromTargets();
      }),
  );
  inputs.push(
    imageFolder
      .addBinding(CONFIG, "imageStyle", {
        options: {
          Fill: "Fill",
          Outline: "Outline",
        },
        label: "Style",
      })
      .on("change", () => {
        if (!hasImageSources()) return;
        if (rebuildImagePointCloud()) {
          rebuildImageTargets();
          if (CONFIG.imageMode) snapParticlesToImage();
        }
        if (rebuildTextPointCloud()) {
          rebuildTextTargets();
          if (CONFIG.imageMode) snapParticlesToText();
        }
        refreshImageModeFromTargets();
      }),
  );
  inputs.push(
    imageFolder
      .addBinding(CONFIG, "imageInvert", {
        label: "Invert",
      })
      .on("change", () => {
        if (!hasImageSources()) return;
        if (rebuildImagePointCloud()) {
          rebuildImageTargets();
          if (CONFIG.imageMode) snapParticlesToImage();
        }
        if (rebuildTextPointCloud()) {
          rebuildTextTargets();
          if (CONFIG.imageMode) snapParticlesToText();
        }
        refreshImageModeFromTargets();
      }),
  );
  inputs.push(
    imageFolder.addBinding(CONFIG, "imageBlend", {
      min: 0.0,
      max: 1.0,
      step: 0.01,
      label: "Image mix",
    }),
  );

  const textFolder = paneSecondary.addFolder({ title: "Text" });
  inputs.push(
    textFolder
      .addBinding(CONFIG, "textValue", {
        label: "Text",
      })
      .on("change", () => {
        applyTextSource();
      }),
  );
  inputs.push(
    textFolder
      .addBinding(CONFIG, "textAlign", {
        options: {
          Left: "Left",
          Center: "Center",
          Right: "Right",
        },
        label: "Align",
      })
      .on("change", () => {
        applyTextSource();
      }),
  );
  inputs.push(
    textFolder
      .addBinding(CONFIG, "textOffsetX", {
        min: -1.5,
        max: 1.5,
        step: 0.01,
        label: "X offset",
      })
      .on("change", () => {
        if (!textImageSource) return;
        if (rebuildTextPointCloud()) {
          rebuildTextTargets();
          if (CONFIG.imageMode) snapParticlesToText();
        }
      }),
  );
  inputs.push(
    textFolder
      .addBinding(CONFIG, "textColor", {
        label: "Text color",
      })
      .on("change", (ev) => applyTextColor(ev.value)),
  );
  inputs.push(
    textFolder
      .addBinding(CONFIG, "textParticleCount", {
        min: 10000,
        max: 200000,
        step: 1000,
        label: "Particle count",
      })
      .on("change", () => {
        rebuildTextParticles();
      }),
  );
  inputs.push(
    textFolder
      .addBinding(CONFIG, "textParticleSize", {
        min: 1.0,
        max: 4.0,
        step: 0.1,
        label: "Particle size",
      })
      .on("change", (ev) => {
        if (textPoints?.material) textPoints.material.size = ev.value;
      }),
  );
  inputs.push(
    textFolder
      .addBinding(CONFIG, "textParticleOpacity", {
        min: 0.0,
        max: 1.0,
        step: 0.01,
        label: "Text opacity",
      })
      .on("change", (ev) => {
        if (textPoints?.material) textPoints.material.opacity = ev.value;
        syncLogoVisibility();
      }),
  );
  inputs.push(
    textFolder
      .addBinding(CONFIG, "textFontSize", {
        min: 16,
        max: 360,
        step: 1,
        label: "Font size",
      })
      .on("change", () => {
        applyTextSource();
      }),
  );
  textFolder.addButton({ title: "Apply text" }).on("click", () => {
    applyTextSource();
  });
  textFolder.addButton({ title: "Remove text" }).on("click", () => {
    clearTextData();
  });

  const captureFolder = pane.addFolder({ title: "Capture" });
  inputs.push(
    captureFolder
      .addBinding(CONFIG, "exportBaseSize", {
        min: 1024,
        max: 7168,
        step: 512,
        label: "Image resolution",
      })
      .on("change", () => rebuildParticles()),
  );
  inputs.push(
    captureFolder
      .addBinding(CONFIG, "cameraControl", {
        x: { min: -CONFIG.viewScale, max: CONFIG.viewScale },
        y: { min: -CONFIG.viewScale, max: CONFIG.viewScale },
        z: { min: 0.5, max: 10, step: 0.1 },
      })
      .on("change", (ev) => applyCameraFromControl(ev.value)),
  );

  captureFolder
    .addButton({ title: "Save image" })
    .on("click", () => saveImage({ useViewport: false }));
  captureFolder
    .addButton({ title: "Save current view" })
    .on("click", () => saveImage({ useViewport: true }));
  inputs.push(
    captureFolder.addBinding(CONFIG, "videoLoopSeconds", {
      min: 1,
      max: 12,
      step: 0.5,
      label: "Duration (sec)",
    }),
  );
  inputs.push(
    captureFolder.addBinding(CONFIG, "videoFps", {
      min: 6,
      max: 24,
      step: 1,
      label: "GIF fps",
    }),
  );
  inputs.push(
    captureFolder.addBinding(CONFIG, "gifScale", {
      min: 0.25,
      max: 1.0,
      step: 0.05,
      label: "GIF scale",
    }),
  );
  inputs.push(
    captureFolder.addBinding(CONFIG, "gifQuality", {
      min: 10,
      max: 60,
      step: 1,
      label: "GIF compression",
    }),
  );
  inputs.push(
    captureFolder.addBinding(CONFIG, "formAnimate", {
      label: "Animate formation",
    }),
  );
  inputs.push(
    captureFolder.addBinding(CONFIG, "formAnimSpeed", {
      min: 0.25,
      max: 4.0,
      step: 0.25,
      label: "Formation speed",
    }),
  );
  captureFolder
    .addButton({ title: "Save GIF" })
    .on("click", () => recordVideoLoop());
  captureFolder
    .addButton({ title: "Save WebM" })
    .on("click", () => recordWebmLoop());

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
  const oldLogoPointSize = logoPoints ? logoPoints.material.size : 0;
  const oldTextPointSize = textPoints ? textPoints.material.size : 0;
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
  renderer.setClearColor(CONFIG.backgroundColor, 1);

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
  if (logoPoints)
    logoPoints.material.size = oldLogoPointSize * (exportWidth / oldSize.x);
  if (textPoints)
    textPoints.material.size = oldTextPointSize * (exportWidth / oldSize.x);
  renderer.render(scene, camera);

  const link = document.createElement("a");
  link.download = `chladni-${Date.now()}-${exportWidth}x${exportHeight}.png`;
  link.href = renderer.domElement.toDataURL("image/png");
  link.click();

  points.material.size = oldPointSize;
  if (logoPoints) logoPoints.material.size = oldLogoPointSize;
  if (textPoints) textPoints.material.size = oldTextPointSize;
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

function loadGifLibrary() {
  if (typeof window === "undefined") {
    return Promise.reject(new Error("GIF encoding unavailable"));
  }
  if (gifLibraryPromise) return gifLibraryPromise;

  const candidates = [
    "https://cdn.jsdelivr.net/npm/gifenc@1.0.3/+esm",
    "https://esm.sh/gifenc@1.0.3",
  ];
  gifLibraryPromise = (async () => {
    let lastError = null;
    for (const url of candidates) {
      try {
        const mod = await import(url);
        const source =
          mod?.default && typeof mod.default === "object" ? mod.default : mod;
        const GIFEncoder = source?.GIFEncoder || mod?.GIFEncoder;
        const quantize = source?.quantize || mod?.quantize;
        const applyPalette = source?.applyPalette || mod?.applyPalette;
        if (GIFEncoder && quantize && applyPalette) {
          return { GIFEncoder, quantize, applyPalette };
        }
        lastError = new Error("GIF library missing required exports");
      } catch (error) {
        lastError = error;
      }
    }
    throw lastError || new Error("Failed to load GIF library");
  })().catch((error) => {
    gifLibraryPromise = null;
    throw error;
  });

  return gifLibraryPromise;
}

function ensureCaptureStatusEl() {
  if (captureStatusEl) return captureStatusEl;
  const el = document.createElement("div");
  el.style.position = "fixed";
  el.style.left = "12px";
  el.style.bottom = "12px";
  el.style.padding = "8px 10px";
  el.style.borderRadius = "6px";
  el.style.background = "rgba(0,0,0,0.75)";
  el.style.color = "#d6e8c5";
  el.style.font = "12px/1.2 ui-sans-serif, system-ui, -apple-system, Segoe UI";
  el.style.zIndex = "9999";
  el.style.pointerEvents = "none";
  el.style.opacity = "0";
  el.style.transition = "opacity 0.15s ease";
  document.body.appendChild(el);
  captureStatusEl = el;
  return el;
}

function setCaptureStatus(message, timeoutMs = 1600) {
  const el = ensureCaptureStatusEl();
  el.textContent = message;
  el.style.opacity = "1";
  if (captureStatusTimer) {
    clearTimeout(captureStatusTimer);
    captureStatusTimer = null;
  }
  if (timeoutMs > 0) {
    captureStatusTimer = window.setTimeout(() => {
      if (!captureStatusEl) return;
      captureStatusEl.style.opacity = "0";
      captureStatusTimer = null;
    }, timeoutMs);
  }
}

function downloadBlob(blob, filename) {
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  link.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function getSupportedWebmMimeType() {
  if (typeof MediaRecorder === "undefined") return "";
  const options = [
    "video/webm;codecs=vp9",
    "video/webm;codecs=vp8",
    "video/webm",
  ];
  for (const option of options) {
    if (MediaRecorder.isTypeSupported(option)) return option;
  }
  return "";
}

function stopVideoLoopRecording() {
  if (gifRecordingTimer) {
    clearTimeout(gifRecordingTimer);
    gifRecordingTimer = null;
  }
  if (gifFrameTimer) {
    clearTimeout(gifFrameTimer);
    gifFrameTimer = null;
  }
  if (gifRecorder) {
    gifRecorder.cancelled = true;
    gifRecorder = null;
    gifIsRendering = false;
    setCaptureStatus("GIF capture canceled");
  }
}

async function recordVideoLoop() {
  if (!renderer || !renderer.domElement) return;
  if (gifRecorder || gifIsRendering || webmRecorder) {
    setCaptureStatus("Recording already in progress");
    return;
  }
  resetAllParticles();
  if (CONFIG.formAnimate) {
    scatterParticles();
    formAnimTime = 0;
  }

  let gifenc;
  try {
    gifenc = await loadGifLibrary();
  } catch (error) {
    console.error(error);
    setCaptureStatus("GIF encoder unavailable");
    window.alert("Could not load GIF encoder in this browser.");
    return;
  }
  const { GIFEncoder, quantize, applyPalette } = gifenc;

  const durationSeconds = clamp(Number(CONFIG.videoLoopSeconds) || 6, 1, 12);
  const fps = clamp(Math.round(Number(CONFIG.videoFps) || 12), 6, 24);
  const gifScale = clamp(Number(CONFIG.gifScale) || 1, 0.25, 1.0);
  const gifQuality = clamp(Math.round(Number(CONFIG.gifQuality) || 20), 10, 60);
  const frameDelayMs = Math.max(8, Math.round(1000 / fps));
  const totalFrames = Math.max(1, Math.round(durationSeconds * fps));
  const baseWidth = Math.max(16, Math.round(Number(CONFIG.videoWidth) || 1280));
  const baseHeight = Math.max(
    16,
    Math.round(Number(CONFIG.videoHeight) || 720),
  );
  const targetWidth = Math.max(16, Math.round(baseWidth * gifScale));
  const targetHeight = Math.max(16, Math.round(baseHeight * gifScale));
  const oldSize = new THREE.Vector2();
  renderer.getSize(oldSize);
  const oldPixelRatio = renderer.getPixelRatio();
  const oldPointSize = points.material.size;
  const oldLogoPointSize = logoPoints ? logoPoints.material.size : 0;
  const oldTextPointSize = textPoints ? textPoints.material.size : 0;
  let captureStateApplied = false;

  const restoreAfterCapture = () => {
    if (!captureStateApplied) return;
    captureStateApplied = false;
    points.material.size = oldPointSize;
    if (logoPoints) logoPoints.material.size = oldLogoPointSize;
    if (textPoints) textPoints.material.size = oldTextPointSize;
    renderer.setPixelRatio(oldPixelRatio);
    renderer.setSize(oldSize.x, oldSize.y, false);
    onWindowResize();
  };

  renderer.setPixelRatio(1);
  renderer.setSize(targetWidth, targetHeight, false);
  {
    const aspect = targetWidth / targetHeight;
    const frustumSize = CONFIG.viewScale * 2;
    camera.left = (-frustumSize * aspect) / 2;
    camera.right = (frustumSize * aspect) / 2;
    camera.top = frustumSize / 2;
    camera.bottom = -frustumSize / 2;
    camera.updateProjectionMatrix();
  }
  const pointScale = targetWidth / Math.max(1, oldSize.x);
  points.material.size = oldPointSize * pointScale;
  if (logoPoints) logoPoints.material.size = oldLogoPointSize * pointScale;
  if (textPoints) textPoints.material.size = oldTextPointSize * pointScale;
  captureStateApplied = true;

  const startTimestamp = Date.now();
  const state = { cancelled: false };
  gifRecorder = state;
  gifIsRendering = true;

  const tempCanvas = document.createElement("canvas");
  tempCanvas.width = targetWidth;
  tempCanvas.height = targetHeight;
  const tempCtx = tempCanvas.getContext("2d", { willReadFrequently: true });
  if (!tempCtx) {
    restoreAfterCapture();
    gifRecorder = null;
    gifIsRendering = false;
    setCaptureStatus("GIF context unavailable");
    return;
  }

  let encoder;
  try {
    encoder = GIFEncoder();
  } catch (error) {
    console.error("Could not start GIF encoder", error);
    restoreAfterCapture();
    gifRecorder = null;
    gifIsRendering = false;
    setCaptureStatus("GIF encoder failed");
    window.alert("Could not start GIF export in this browser.");
    return;
  }

  let capturedFrames = 0;
  let palette = null;
  const maxColors = clamp(
    Math.round(256 / Math.pow(2, (gifQuality - 10) / 10)),
    8,
    256,
  );
  const paletteFormat = gifQuality >= 20 ? "rgb444" : "rgb565";

  const finishCapture = () => {
    if (gifRecordingTimer) {
      clearTimeout(gifRecordingTimer);
      gifRecordingTimer = null;
    }
    if (gifFrameTimer) {
      clearTimeout(gifFrameTimer);
      gifFrameTimer = null;
    }
    if (!state.cancelled) {
      try {
        encoder.finish();
        const bytes =
          typeof encoder.bytesView === "function"
            ? encoder.bytesView()
            : encoder.bytes();
        const blob = new Blob([bytes], { type: "image/gif" });
        downloadBlob(
          blob,
          `chladni-loop-${startTimestamp}-${targetWidth}x${targetHeight}.gif`,
        );
        setCaptureStatus("GIF saved");
      } catch (error) {
        console.error("GIF finalization failed", error);
        setCaptureStatus("GIF export failed");
      }
    }
    if (gifRecorder === state) {
      gifRecorder = null;
    }
    gifIsRendering = false;
    restoreAfterCapture();
  };

  const captureFrame = () => {
    if (state.cancelled) {
      finishCapture();
      return;
    }
    try {
      renderer.render(scene, camera);
      tempCtx.clearRect(0, 0, targetWidth, targetHeight);
      tempCtx.drawImage(renderer.domElement, 0, 0, targetWidth, targetHeight);
      const rgba = tempCtx.getImageData(0, 0, targetWidth, targetHeight).data;

      if (!palette) {
        palette = quantize(rgba, maxColors, { format: paletteFormat });
      }
      const index = applyPalette(rgba, palette, paletteFormat);
      encoder.writeFrame(index, targetWidth, targetHeight, {
        palette,
        delay: frameDelayMs,
        repeat: 0,
      });
    } catch (error) {
      console.error("GIF frame encoding failed", error);
      state.cancelled = true;
      setCaptureStatus("GIF recording failed");
      finishCapture();
      return;
    }

    capturedFrames++;
    if (capturedFrames >= totalFrames) {
      finishCapture();
      return;
    }
    setCaptureStatus(`Recording GIF ${capturedFrames}/${totalFrames}...`, 0);
    gifFrameTimer = window.setTimeout(captureFrame, frameDelayMs);
  };

  gifRecordingTimer = window.setTimeout(
    () => {
      if (!state.cancelled) {
        state.cancelled = true;
        setCaptureStatus("GIF capture timeout");
      }
      finishCapture();
    },
    Math.round(durationSeconds * 1000 * 2.5),
  );
  setCaptureStatus(`Recording GIF 0/${totalFrames}...`, 0);
  captureFrame();
}

function stopWebmLoopRecording() {
  if (webmRecordingTimer) {
    clearTimeout(webmRecordingTimer);
    webmRecordingTimer = null;
  }
  if (webmRecorder && webmRecorder.state !== "inactive") {
    webmRecorder.stop();
  }
}

function recordWebmLoop() {
  if (!renderer || !renderer.domElement) return;
  if (gifRecorder || gifIsRendering || webmRecorder) {
    setCaptureStatus("Recording already in progress");
    return;
  }
  resetAllParticles();
  if (CONFIG.formAnimate) {
    scatterParticles();
    formAnimTime = 0;
  }
  if (typeof MediaRecorder === "undefined") {
    setCaptureStatus("WebM recording unsupported");
    window.alert("WebM recording is not supported in this browser.");
    return;
  }
  if (typeof renderer.domElement.captureStream !== "function") {
    setCaptureStatus("Canvas capture unsupported");
    window.alert("Canvas stream capture is not available in this browser.");
    return;
  }

  const durationSeconds = clamp(Number(CONFIG.videoLoopSeconds) || 6, 1, 12);
  const fps = clamp(Math.round(Number(CONFIG.videoFps) || 12), 6, 60);
  const mimeType = getSupportedWebmMimeType();
  const targetWidth = Math.max(
    16,
    Math.round(Number(CONFIG.videoWidth) || 1280),
  );
  const targetHeight = Math.max(
    16,
    Math.round(Number(CONFIG.videoHeight) || 720),
  );
  const oldSize = new THREE.Vector2();
  renderer.getSize(oldSize);
  const oldPixelRatio = renderer.getPixelRatio();
  const oldPointSize = points.material.size;
  const oldLogoPointSize = logoPoints ? logoPoints.material.size : 0;
  const oldTextPointSize = textPoints ? textPoints.material.size : 0;
  let captureStateApplied = false;

  const restoreAfterCapture = () => {
    if (!captureStateApplied) return;
    captureStateApplied = false;
    points.material.size = oldPointSize;
    if (logoPoints) logoPoints.material.size = oldLogoPointSize;
    if (textPoints) textPoints.material.size = oldTextPointSize;
    renderer.setPixelRatio(oldPixelRatio);
    renderer.setSize(oldSize.x, oldSize.y, false);
    onWindowResize();
  };

  renderer.setPixelRatio(1);
  renderer.setSize(targetWidth, targetHeight, false);
  {
    const aspect = targetWidth / targetHeight;
    const frustumSize = CONFIG.viewScale * 2;
    camera.left = (-frustumSize * aspect) / 2;
    camera.right = (frustumSize * aspect) / 2;
    camera.top = frustumSize / 2;
    camera.bottom = -frustumSize / 2;
    camera.updateProjectionMatrix();
  }
  const pointScale = targetWidth / Math.max(1, oldSize.x);
  points.material.size = oldPointSize * pointScale;
  if (logoPoints) logoPoints.material.size = oldLogoPointSize * pointScale;
  if (textPoints) textPoints.material.size = oldTextPointSize * pointScale;
  captureStateApplied = true;

  const stream = renderer.domElement.captureStream(fps);
  const startTimestamp = Date.now();
  let recorder;

  webmRecordingChunks = [];
  try {
    recorder = mimeType
      ? new MediaRecorder(stream, {
          mimeType,
          videoBitsPerSecond: 12000000,
        })
      : new MediaRecorder(stream, {
          videoBitsPerSecond: 12000000,
        });
  } catch (error) {
    console.error("Could not start WebM recorder", error);
    stream.getTracks().forEach((track) => track.stop());
    restoreAfterCapture();
    webmRecorder = null;
    setCaptureStatus("WebM recorder failed");
    window.alert("Could not start WebM recording in this browser.");
    return;
  }

  webmRecorder = recorder;
  recorder.addEventListener("dataavailable", (event) => {
    if (event.data && event.data.size > 0) {
      webmRecordingChunks.push(event.data);
    }
  });

  recorder.addEventListener("stop", () => {
    stream.getTracks().forEach((track) => track.stop());
    if (webmRecordingTimer) {
      clearTimeout(webmRecordingTimer);
      webmRecordingTimer = null;
    }
    restoreAfterCapture();
    if (!webmRecordingChunks.length) {
      webmRecorder = null;
      return;
    }
    const actualMime = recorder.mimeType || mimeType || "video/webm";
    const blob = new Blob(webmRecordingChunks, { type: actualMime });
    downloadBlob(
      blob,
      `chladni-loop-${startTimestamp}-${targetWidth}x${targetHeight}.webm`,
    );
    webmRecordingChunks = [];
    webmRecorder = null;
    setCaptureStatus("WebM saved");
  });

  recorder.addEventListener("error", (event) => {
    console.error("WebM recorder error", event);
    stream.getTracks().forEach((track) => track.stop());
    if (webmRecordingTimer) {
      clearTimeout(webmRecordingTimer);
      webmRecordingTimer = null;
    }
    restoreAfterCapture();
    webmRecordingChunks = [];
    webmRecorder = null;
    setCaptureStatus("WebM recording failed");
  });

  recorder.start(250);
  setCaptureStatus(`Recording WebM ${targetWidth}x${targetHeight}...`, 0);
  webmRecordingTimer = window.setTimeout(
    () => stopWebmLoopRecording(),
    Math.round(durationSeconds * 1000),
  );
}

function randomizeWaves(rebuild = true) {
  const keys = WAVE_TYPE_KEYS;
  CONFIG.waveTypeA = keys[Math.floor(Math.random() * keys.length)];
  CONFIG.waveTypeB = keys[Math.floor(Math.random() * keys.length)];

  if (rebuild) {
    rebuildField();
    resetParticles();
    if (refreshUI) refreshUI();
  }
}

function randomizeAll() {
  randomizeModes(false);
  randomizeWaves(false);

  CONFIG.waveMix = Math.round(Math.random() * 100) / 100;
  CONFIG.fieldScale = Math.round((0.2 + Math.random() * 1.8) * 10) / 10;

  // Randomize particles
  const minParticleCount = 50000;
  const maxParticleCount = 500000;
  const particleCountStep = 10000;
  const particleCountSteps = Math.floor(
    (maxParticleCount - minParticleCount) / particleCountStep,
  );
  CONFIG.particleCount =
    minParticleCount +
    Math.floor(Math.random() * (particleCountSteps + 1)) * particleCountStep;

  CONFIG.particleSize = Math.round((1.0 + Math.random() * 4.0) * 10) / 10;
  if (points) points.material.size = CONFIG.particleSize;

  CONFIG.particleOpacity = Math.round((0.1 + Math.random() * 0.9) * 100) / 100;
  if (points) points.material.opacity = CONFIG.particleOpacity;

  CONFIG.color =
    "#" +
    new THREE.Color(Math.random(), Math.random(), Math.random()).getHexString();

  rebuildField();
  rebuildParticles();
  if (refreshUI) refreshUI();
}
