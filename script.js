import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { GUI } from "lil-gui";

// --- Configuration ---
const CONFIG = {
  particleCount: 200000,
  gridSize: 256, // Resolution of the underlying energy field
  settleStrength: 1.5, // How fast particles move to nodal lines
  jitter: 0.2, // Brownian motion to keep them "alive"
  drag: 0.88, // Friction
  speedLimit: 2.0,
  viewScale: 700, // Size of the world
  color: "#9fd3ff",
  particleSize: 2.5,
  fogDensity: 0.0005,
  cameraZoom: 1,
  captureScale: 2,
  exportBaseSize: 2048,
  recordFps: 30,
  recordWidth: 1280,
  recordHeight: 720,
  modeCount: 4,
  modeMin: 2,
  modeMax: 12,
  shape: "Rectangle",
  rectAspect: 1,
  waveTypeA: "Cartesian",
  waveTypeB: "Radial",
  waveMix: 0.5,
  randomPhase: true,
  randomRotation: true,
  rotationMaxDeg: 30,
  modeScaleJitter: 0.15,
  modeStyle: "Balanced",
  // Wave types: "Cartesian", "Radial", "Spiral", "Hexagonal", "Quasicrystal", "Lissajous", "Moire", "Chevron", "Ring", "Flower", "GridWarp", "Diamond", "Square", "Lattice", "Kaleidoscope", "Voronoi", "Perlin"
};

// --- Global Variables ---
let scene, camera, renderer, geometry, points;
let positions, velocities; // Float32Arrays for high performance
let frameId;
let recorder = null;
let recordingChunks = [];
let isRecording = false;
let recordCanvas = null;
let recordCtx = null;

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
  scene.fog = null;

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
  camera.zoom = CONFIG.cameraZoom;
  camera.lookAt(0, 0, 0);
  camera.updateProjectionMatrix();

  // 2. Renderer
  renderer = new THREE.WebGLRenderer({
    antialias: true,
    powerPreference: "high-performance",
    preserveDrawingBuffer: true,
  });
  renderer.setSize(window.innerWidth, window.innerHeight);
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2)); // Limit pixel ratio for performance
  renderer.setClearColor(0x000000, 1);
  document.body.appendChild(renderer.domElement);

  // 3. Controls
  const controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  controls.dampingFactor = 0.05;
  controls.maxDistance = 1500;
  controls.minDistance = 100;
  controls.target.set(0, 0, 0);
  controls.enableRotate = false; // Keep top-down view locked
  controls.minPolarAngle = Math.PI / 2;
  controls.maxPolarAngle = Math.PI / 2;

  // 4. Particle System Setup
  buildParticles();

  // 5. Material
  // We use a simple circular sprite texture for efficiency
  const sprite = new THREE.TextureLoader().load(
    "https://threejs.org/examples/textures/sprites/disc.png",
  );

  const material = new THREE.PointsMaterial({
    color: CONFIG.color,
    vertexColors: false,
    size: CONFIG.particleSize,
    map: sprite,
    alphaTest: 0.2, // Performance optimization
    transparent: false,
    opacity: 1,
    blending: THREE.AdditiveBlending, // Makes overlapping particles glow
    depthWrite: false, // Crucial for transparency performance
    sizeAttenuation: false, // Keep points visible when zooming out
  });

  points = new THREE.Points(geometry, material);
  scene.add(points);

  // 6. Event Listeners
  window.addEventListener("resize", onWindowResize);
  // Randomization is now only via GUI

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

  // Allocate memory
  positions = new Float32Array(CONFIG.particleCount * 3); // x, y, z
  velocities = new Float32Array(CONFIG.particleCount * 2); // vx, vy (we don't need vz physics)

  // Initialize particles randomly within the chosen shape
  for (let i = 0; i < CONFIG.particleCount; i++) {
    const [x, y] = randomPointInShape();
    positions[i * 3 + 0] = x;
    positions[i * 3 + 1] = y;
    positions[i * 3 + 2] = 0; // z

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
  // Calculates the standing wave pattern on a grid
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
  // Helpers below used by select wave types.

  const baseCartesian = (cx, cy, mode) => {
    const cosR = mode.cos;
    const sinR = mode.sin;
    const rx = (cx * cosR - cy * sinR) * mode.sx;
    const ry = (cx * sinR + cy * cosR) * mode.sy;
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
    if (type === "Hexagonal") {
      const cosR = mode.cos;
      const sinR = mode.sin;
      const rx = (cx * cosR - cy * sinR) * mode.sx;
      const ry = (cx * sinR + cy * cosR) * mode.sy;
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
      const rx = (cx * cosR - cy * sinR) * mode.sx;
      const ry = (cx * sinR + cy * cosR) * mode.sy;
      const kFreq = mode.m * Math.PI;
      let sum = 0;
      for (let j = 0; j < 5; j++) {
        const theta = (Math.PI * 2 * j) / 5;
        const v = rx * Math.cos(theta) + ry * Math.sin(theta);
        sum += Math.cos(kFreq * v + mode.px);
      }
      return sum;
    }
    if (type === "Lissajous") {
      const cosR = mode.cos;
      const sinR = mode.sin;
      const rx = (cx * cosR - cy * sinR) * mode.sx;
      const ry = (cx * sinR + cy * cosR) * mode.sy;
      return (
        Math.sin(mode.m * Math.PI * rx + mode.px) *
        Math.sin(mode.n * Math.PI * ry + mode.py)
      );
    }
    if (type === "Moire") {
      const cosR = mode.cos;
      const sinR = mode.sin;
      const rx = (cx * cosR - cy * sinR) * mode.sx;
      const ry = (cx * sinR + cy * cosR) * mode.sy;
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
      const rx = (cx * cosR - cy * sinR) * mode.sx;
      const ry = (cx * sinR + cy * cosR) * mode.sy;
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
    if (type === "Flower") {
      const r = Math.sqrt(cx * cx + cy * cy) * 2.0;
      const theta = Math.atan2(cy, cx);
      const nInt = Math.round(mode.n);
      return (
        Math.sin(mode.m * Math.PI * r + mode.px) *
        Math.cos(nInt * theta + mode.py)
      );
    }
    if (type === "GridWarp") {
      const cosR = mode.cos;
      const sinR = mode.sin;
      const rx = (cx * cosR - cy * sinR) * mode.sx;
      const ry = (cx * sinR + cy * cosR) * mode.sy;
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
      const rx = (cx * cosR - cy * sinR) * mode.sx;
      const ry = (cx * sinR + cy * cosR) * mode.sy;
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
      const rx = (cx * cosR - cy * sinR) * mode.sx;
      const ry = (cx * sinR + cy * cosR) * mode.sy;
      const scale = Math.max(1.5, mode.m);
      const n = valueNoise(rx * scale + mode.px, ry * scale + mode.py);
      return n;
    }

    // Cartesian
    return baseCartesian(cx, cy, mode);
  };

  // 1. Compute Energy Field (Amplitude^2)
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
        const wave = waveA * (1 - mix) + waveB * mix;
        phi += mode.a * wave;
      }
      energy[idx] = phi * phi;
    }
  }

  // 2. Compute Gradients (dE/dx, dE/dy)
  // This tells particles which way is "downhill" towards the nodal lines
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

  // Optimization: Pre-calculate constants to avoid division and recalculation inside the loop
  const gridScale = (G - 1) / fullRange;
  const shape = CONFIG.shape;
  const rectRy = Math.max(1, Math.round(range / CONFIG.rectAspect));
  const ovalRy = Math.max(1, Math.round(range * 0.65));
  const rangeSq = range * range;
  const ovalRySq = ovalRy * ovalRy;

  for (let i = 0; i < count; i++) {
    const i3 = i * 3;
    const i2 = i * 2;

    let x = pos[i3];
    let y = pos[i3 + 1];
    let vx = vel[i2];
    let vy = vel[i2 + 1];

    // 1. Boundary Check (Inlined for performance)
    let inside = true;
    if (shape === "Rectangle") {
      if (x < -range || x > range || y < -rectRy || y > rectRy) inside = false;
    } else if (shape === "Circle") {
      if (x * x + y * y > rangeSq) inside = false;
    } else if (shape === "Oval") {
      const nx = x / range;
      const ny = y / ovalRy;
      if (nx * nx + ny * ny > 1) inside = false;
    }

    // Teleport if out of bounds
    if (!inside) {
      const respawn = randomPointInShape();
      x = respawn[0];
      y = respawn[1];
      vx = 0;
      vy = 0;
    }

    // 2. Bilinear Interpolation of Gradients
    // Convert world position to grid coordinates
    let gx_pos = (x + range) * gridScale;
    let gy_pos = (y + range) * gridScale;

    // Clamp to grid edges
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

    // Interpolate Gradient X
    const gxVal =
      (gradX[idx00] * (1 - tx) + gradX[idx10] * tx) * (1 - ty) +
      (gradX[idx01] * (1 - tx) + gradX[idx11] * tx) * ty;

    // Interpolate Gradient Y
    const gyVal =
      (gradY[idx00] * (1 - tx) + gradY[idx10] * tx) * (1 - ty) +
      (gradY[idx01] * (1 - tx) + gradY[idx11] * tx) * ty;

    // Interpolate Energy (for Z height)
    const eVal =
      (energy[idx00] * (1 - tx) + energy[idx10] * tx) * (1 - ty) +
      (energy[idx01] * (1 - tx) + energy[idx11] * tx) * ty;

    // 3. Update Velocity (Gradient Descent)
    // Move opposite to gradient (downhill)
    vx -= gxVal * CONFIG.settleStrength;
    vy -= gyVal * CONFIG.settleStrength;

    // Jitter (Brownian motion)
    vx += (Math.random() - 0.5) * CONFIG.jitter;
    vy += (Math.random() - 0.5) * CONFIG.jitter;

    // Apply Drag
    vx *= CONFIG.drag;
    vy *= CONFIG.drag;

    // Speed Limit
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
    pos[i3 + 2] = eVal * 2; // Minimal Z height for 2D plate look

    vel[i2] = vx;
    vel[i2 + 1] = vy;
  }

  // Mark geometry as needing an update on the GPU
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

function isInsideShape(x, y) {
  const range = CONFIG.viewScale;
  const rx = range;
  const rectRy = Math.max(1, Math.round(range / CONFIG.rectAspect));
  const ry =
    CONFIG.shape === "Oval" ? Math.max(1, Math.round(range * 0.65)) : range;

  if (CONFIG.shape === "Circle") {
    return x * x + y * y <= range * range;
  }

  if (CONFIG.shape === "Oval") {
    const nx = x / rx;
    const ny = y / ry;
    return nx * nx + ny * ny <= 1;
  }

  if (CONFIG.shape === "Rectangle") {
    return x >= -range && x <= range && y >= -rectRy && y <= rectRy;
  }

  // Rectangle (default)
  return x >= -range && x <= range && y >= -rectRy && y <= rectRy;
}

function randomPointInShape() {
  const range = CONFIG.viewScale;
  const rx = range;
  const rectRy = Math.max(1, Math.round(range / CONFIG.rectAspect));
  const ry =
    CONFIG.shape === "Oval" ? Math.max(1, Math.round(range * 0.65)) : range;

  // Rejection sampling for circle/oval
  if (CONFIG.shape === "Circle" || CONFIG.shape === "Oval") {
    while (true) {
      const x = (Math.random() - 0.5) * 2 * rx;
      const y = (Math.random() - 0.5) * 2 * ry;
      if (isInsideShape(x, y)) {
        return [x, y];
      }
    }
  }

  // Rectangle (default)
  return [
    (Math.random() - 0.5) * 2 * range,
    (Math.random() - 0.5) * 2 * rectRy,
  ];
}

function randomizeModes() {
  // Generate 3 random wave modes
  modes = [];
  for (let i = 0; i < CONFIG.modeCount; i++) {
    const rot =
      CONFIG.randomRotation && CONFIG.rotationMaxDeg > 0
        ? ((Math.random() * 2 - 1) * CONFIG.rotationMaxDeg * Math.PI) / 180
        : 0;
    const sx = 1 + (Math.random() * 2 - 1) * CONFIG.modeScaleJitter;
    const sy = 1 + (Math.random() * 2 - 1) * CONFIG.modeScaleJitter;
    const mRaw =
      Math.random() ** 0.6 * (CONFIG.modeMax - CONFIG.modeMin) + CONFIG.modeMin;
    const nRaw =
      Math.random() ** 0.6 * (CONFIG.modeMax - CONFIG.modeMin) + CONFIG.modeMin;
    const m = mRaw;
    const n = nRaw;
    modes.push({
      m,
      n,
      a: Math.random() * (1 - i * 0.2), // Dampen amplitude for secondary modes
      px: CONFIG.randomPhase ? Math.random() * Math.PI * 2 : 0,
      py: CONFIG.randomPhase ? Math.random() * Math.PI * 2 : 0,
      cos: Math.cos(rot),
      sin: Math.sin(rot),
      sx,
      sy,
    });
  }

  normalizeModeAmplitudes();
  rebuildField();

  // "Kick" the particles so they find the new pattern
  for (let i = 0; i < CONFIG.particleCount; i++) {
    velocities[i * 2] += (Math.random() - 0.5) * 5;
    velocities[i * 2 + 1] += (Math.random() - 0.5) * 5;
  }
}

function initModes() {
  modes = [];
  for (let i = 0; i < CONFIG.modeCount; i++) {
    const rot =
      CONFIG.randomRotation && CONFIG.rotationMaxDeg > 0
        ? ((Math.random() * 2 - 1) * CONFIG.rotationMaxDeg * Math.PI) / 180
        : 0;
    const sx = 1 + (Math.random() * 2 - 1) * CONFIG.modeScaleJitter;
    const sy = 1 + (Math.random() * 2 - 1) * CONFIG.modeScaleJitter;
    const mRaw =
      Math.random() ** 0.6 * (CONFIG.modeMax - CONFIG.modeMin) + CONFIG.modeMin;
    const nRaw =
      Math.random() ** 0.6 * (CONFIG.modeMax - CONFIG.modeMin) + CONFIG.modeMin;
    const m = mRaw;
    const n = nRaw;
    modes.push({
      m,
      n,
      a: Math.random() * (1 - i * 0.2),
      px: CONFIG.randomPhase ? Math.random() * Math.PI * 2 : 0,
      py: CONFIG.randomPhase ? Math.random() * Math.PI * 2 : 0,
      cos: Math.cos(rot),
      sin: Math.sin(rot),
      sx,
      sy,
    });
  }
  normalizeModeAmplitudes();
}

function normalizeModeAmplitudes() {
  let sum = 0;
  for (let i = 0; i < modes.length; i++) {
    sum += Math.abs(modes[i].a);
  }
  if (sum <= 0) return;
  const inv = 1 / sum;
  for (let i = 0; i < modes.length; i++) {
    modes[i].a *= inv;
  }
}

function animate() {
  frameId = requestAnimationFrame(animate);
  updateParticles();
  renderer.render(scene, camera);
  if (isRecording && recordCtx && recordCanvas) {
    recordCtx.drawImage(
      renderer.domElement,
      0,
      0,
      recordCanvas.width,
      recordCanvas.height,
    );
  }
}

function setupGUI() {
  const gui = new GUI({ width: 320 });
  const controllers = [];
  let modeMinController;
  let modeMaxController;
  let modeCountController;

  const sim = gui.addFolder("Simulation");
  controllers.push(
    sim.add(CONFIG, "particleCount", 10000, 500000, 1000).onFinishChange(() => {
      rebuildParticles();
    }),
  );
  controllers.push(sim.add(CONFIG, "settleStrength", 0.1, 4.0, 0.01));
  controllers.push(sim.add(CONFIG, "jitter", 0.0, 0.25, 0.01));
  controllers.push(sim.add(CONFIG, "drag", 0.7, 0.9, 0.01));
  controllers.push(sim.add(CONFIG, "speedLimit", 0.5, 10.0, 0.1));
  // viewScale is fixed via CONFIG default
  controllers.push(
    sim.add(CONFIG, "viewScale", 300, 800, 1).onFinishChange(() => {
      rebuildParticles();
      rebuildField();
      camera.zoom = CONFIG.cameraZoom;
      camera.updateProjectionMatrix();
      onWindowResize();
    }),
  );
  controllers.push(
    sim
      .add(CONFIG, "shape", ["Rectangle", "Circle", "Oval"])
      .onFinishChange(() => {
        rebuildParticles();
      }),
  );
  controllers.push(
    sim
      .add(CONFIG, "waveTypeA", [
        "Cartesian",
        "Radial",
        "Spiral",
        "Hexagonal",
        "Quasicrystal",
        "Lissajous",
        "Moire",
        "Chevron",
        "Ring",
        "Flower",
        "GridWarp",
        "Diamond",
        "Square",
        "Lattice",
        "Kaleidoscope",
        "Voronoi",
        "Perlin",
      ])
      .onFinishChange(() => {
        rebuildField();
      }),
  );
  controllers.push(
    sim
      .add(CONFIG, "waveTypeB", [
        "Cartesian",
        "Radial",
        "Spiral",
        "Hexagonal",
        "Quasicrystal",
        "Lissajous",
        "Moire",
        "Chevron",
        "Ring",
        "Flower",
        "GridWarp",
        "Diamond",
        "Square",
        "Lattice",
        "Kaleidoscope",
        "Voronoi",
        "Perlin",
      ])
      .onFinishChange(() => {
        rebuildField();
      }),
  );
  controllers.push(
    sim.add(CONFIG, "waveMix", 0, 1, 0.01).onFinishChange(() => {
      rebuildField();
    }),
  );
  controllers.push(
    sim.add(CONFIG, "rectAspect", 1, 3, 0.01).onFinishChange(() => {
      rebuildParticles();
    }),
  );

  const look = gui.addFolder("Look");
  controllers.push(
    look.addColor(CONFIG, "color").onChange((value) => {
      points.material.color.set(value);
    }),
  );
  controllers.push(
    look.add(CONFIG, "particleSize", 1, 3, 0.1).onChange((value) => {
      points.material.size = value;
    }),
  );

  const cameraFolder = gui.addFolder("Camera");
  controllers.push(
    cameraFolder.add(CONFIG, "cameraZoom", 0.5, 5, 0.1).onChange((value) => {
      camera.zoom = value;
      camera.updateProjectionMatrix();
    }),
  );

  const modesFolder = gui.addFolder("Modes");
  const modeStyles = {
    Clean: "Clean",
    Balanced: "Balanced",
    Complex: "Complex",
  };
  const applyModeStyle = (style) => {
    if (style === "Clean") {
      CONFIG.modeCount = 2;
      CONFIG.modeMin = 1;
      CONFIG.modeMax = 8;
    } else if (style === "Complex") {
      CONFIG.modeCount = 8;
      CONFIG.modeMin = 3;
      CONFIG.modeMax = 20;
    } else {
      CONFIG.modeCount = 4;
      CONFIG.modeMin = 2;
      CONFIG.modeMax = 12;
    }
  };
  controllers.push(
    modesFolder.add(CONFIG, "modeStyle", modeStyles).onChange((value) => {
      applyModeStyle(value);
      if (modeCountController) modeCountController.updateDisplay();
      if (modeMinController) modeMinController.updateDisplay();
      if (modeMaxController) modeMaxController.updateDisplay();
      randomizeModes();
    }),
  );
  controllers.push(
    (modeCountController = modesFolder
      .add(CONFIG, "modeCount", 1, 20, 1)
      .onChange(() => {
        randomizeModes();
      })),
  );
  modeMinController = modesFolder
    .add(CONFIG, "modeMin", 0, 20, 0.1)
    .onChange(() => {
      if (CONFIG.modeMin > CONFIG.modeMax) {
        CONFIG.modeMax = CONFIG.modeMin;
        if (modeMaxController) {
          modeMaxController.updateDisplay();
        }
      }
      randomizeModes();
    });
  controllers.push(modeMinController);
  modeMaxController = modesFolder
    .add(CONFIG, "modeMax", 0, 50, 0.1)
    .onChange(() => {
      if (CONFIG.modeMax < CONFIG.modeMin) {
        CONFIG.modeMin = CONFIG.modeMax;
        if (modeMinController) {
          modeMinController.updateDisplay();
        }
      }
      randomizeModes();
    });
  controllers.push(modeMaxController);
  controllers.push(
    modesFolder.add(CONFIG, "randomPhase").onChange(() => {
      randomizeModes();
    }),
  );
  controllers.push(
    modesFolder.add(CONFIG, "randomRotation").onChange(() => {
      randomizeModes();
    }),
  );
  controllers.push(
    modesFolder
      .add(CONFIG, "rotationMaxDeg", 0, 90, 1)
      .name("rotationMaxDeg (deg)")
      .onChange(() => {
        randomizeModes();
      }),
  );
  controllers.push(
    modesFolder
      .add(CONFIG, "modeScaleJitter", 0, 0.5, 0.01)
      .name("modeScaleJitter (scale)")
      .onChange(() => {
        randomizeModes();
      }),
  );

  const captureFolder = gui.addFolder("Capture");
  captureFolder.add({ saveImage }, "saveImage");
  captureFolder.add(CONFIG, "captureScale", 1, 6, 1).name("captureScale (x)");
  captureFolder
    .add(CONFIG, "exportBaseSize", 512, 8192, 128)
    .name("exportBaseSize (px)");
  captureFolder.add({ saveImageHQ }, "saveImageHQ");
  captureFolder.add(CONFIG, "recordFps", 15, 60, 1).name("recordFps");
  captureFolder
    .add(CONFIG, "recordWidth", 320, 3840, 10)
    .name("recordWidth (px)");
  captureFolder
    .add(CONFIG, "recordHeight", 240, 2160, 10)
    .name("recordHeight (px)");
  captureFolder.add({ startRecording }, "startRecording");
  captureFolder.add({ stopRecording }, "stopRecording");

  const actionsFolder = gui.addFolder("Actions");
  actionsFolder.add({ randomizeAll }, "randomizeAll");

  function syncGUI() {
    controllers.forEach((controller) => controller.updateDisplay());
  }

  randomizeAll.syncGUI = syncGUI;

  sim.open();
  look.open();
}

function saveImage() {
  const link = document.createElement("a");
  link.download = `chladni-${Date.now()}.png`;
  link.href = renderer.domElement.toDataURL("image/png");
  link.click();
}

function saveImageHQ() {
  const scale = Math.max(1, Math.round(CONFIG.captureScale));
  const oldSize = new THREE.Vector2();
  renderer.getSize(oldSize);
  const oldPixelRatio = renderer.getPixelRatio();
  const oldClearAlpha = renderer.getClearAlpha();
  const oldClearColor = renderer.getClearColor(new THREE.Color());
  const oldPointSize = points.material.size;

  const baseSize = Math.max(256, Math.round(CONFIG.exportBaseSize));

  const range = CONFIG.viewScale;
  const rectRy = Math.max(1, Math.round(range / CONFIG.rectAspect));
  const ovalRy = Math.max(1, Math.round(range * 0.65));
  const exportWorldWidth = range * 2;
  const exportWorldHeight =
    CONFIG.shape === "Circle"
      ? range * 2
      : CONFIG.shape === "Oval"
        ? ovalRy * 2
        : rectRy * 2;

  const aspect = exportWorldWidth / exportWorldHeight;
  const exportWidth = Math.max(1, Math.round(baseSize * scale));
  const exportHeight = Math.max(1, Math.round(exportWidth / aspect));

  renderer.setPixelRatio(1);
  renderer.setSize(exportWidth, exportHeight, false);
  renderer.setClearColor(oldClearColor, 0);

  const frustumWidth = exportWorldWidth;
  const frustumHeight = exportWorldHeight;
  camera.left = -frustumWidth / 2;
  camera.right = frustumWidth / 2;
  camera.top = frustumHeight / 2;
  camera.bottom = -frustumHeight / 2;
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

function startRecording() {
  if (isRecording) return;
  const width = Math.max(1, Math.round(CONFIG.recordWidth));
  const height = Math.max(1, Math.round(CONFIG.recordHeight));
  recordCanvas = document.createElement("canvas");
  recordCanvas.width = width;
  recordCanvas.height = height;
  recordCtx = recordCanvas.getContext("2d");

  const stream = recordCanvas.captureStream(CONFIG.recordFps);
  const mimeType = MediaRecorder.isTypeSupported("video/webm;codecs=vp9")
    ? "video/webm;codecs=vp9"
    : "video/webm";
  recorder = new MediaRecorder(stream, { mimeType });
  recordingChunks = [];
  recorder.ondataavailable = (event) => {
    if (event.data && event.data.size > 0) {
      recordingChunks.push(event.data);
    }
  };
  recorder.onstop = () => {
    const blob = new Blob(recordingChunks, { type: mimeType });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.download = `chladni-${Date.now()}.webm`;
    link.href = url;
    link.click();
    URL.revokeObjectURL(url);
    recordingChunks = [];
    recorder = null;
    recordCanvas = null;
    recordCtx = null;
  };
  recorder.start();
  isRecording = true;
}

function stopRecording() {
  if (!recorder || !isRecording) return;
  recorder.stop();
  isRecording = false;
}

function randomizeAll() {
  const waveTypes = [
    "Cartesian",
    "Radial",
    "Spiral",
    "Hexagonal",
    "Quasicrystal",
    "Lissajous",
    "Moire",
    "Chevron",
    "Ring",
    "Flower",
    "GridWarp",
    "Diamond",
    "Square",
    "Lattice",
    "Kaleidoscope",
    "Voronoi",
    "Perlin",
  ];

  CONFIG.particleCount = Math.floor(10000 + Math.random() * 490000);
  CONFIG.gridSize = 512;
  CONFIG.settleStrength = 0.1 + Math.random() * 3.9;
  CONFIG.jitter = Math.random() * 0.25;
  CONFIG.drag = 0.7 + Math.random() * 0.2;
  CONFIG.speedLimit = 0.5 + Math.random() * 9.5;
  CONFIG.viewScale = 300 + Math.random() * 500;
  CONFIG.waveTypeA = waveTypes[Math.floor(Math.random() * waveTypes.length)];
  CONFIG.waveTypeB = waveTypes[Math.floor(Math.random() * waveTypes.length)];
  CONFIG.color = `#${Math.floor(Math.random() * 0xffffff)
    .toString(16)
    .padStart(6, "0")}`;
  CONFIG.particleSize = 1 + Math.random() * 2;
  CONFIG.particleOpacity = 0.2 + Math.random() * 0.8;
  CONFIG.cameraZoom = 0.5 + Math.random() * 4.5;
  CONFIG.modeCount = Math.floor(1 + Math.random() * 20);
  CONFIG.modeMin = Math.random() * 20;
  CONFIG.modeMax = CONFIG.modeMin + Math.random() * (50 - CONFIG.modeMin);
  CONFIG.randomPhase = Math.random() < 0.5;
  CONFIG.randomRotation = Math.random() < 0.5;
  CONFIG.rotationMaxDeg = Math.floor(Math.random() * 91);
  CONFIG.modeScaleJitter = Math.random() * 0.5;
  CONFIG.waveMix = Math.random();

  rebuildParticles();
  allocateField();
  rebuildField();
  randomizeModes();

  points.material.color.set(CONFIG.color);
  points.material.size = CONFIG.particleSize;

  camera.zoom = CONFIG.cameraZoom;
  camera.updateProjectionMatrix();

  if (randomizeAll.syncGUI) {
    randomizeAll.syncGUI();
  }
}
