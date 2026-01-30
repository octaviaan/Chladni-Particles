# Chladni Particle Flow
Generative art simulation built with **Three.js**. This project simulates Chladni figures (nodal patterns formed on vibrating plates) using a particle system that reacts to a computed energy field.
Instead of calculating simple 2D lines, it generates a "terrain" of standing waves. Millions of particles flow downhill (via gradient descent) to settle in the "valleys" (nodes) where vibration is lowest.

## ✨ Key Features
* **Particle System:** Renders up to ~1.3 million particles using `THREE.Points` and low-level `Float32Array`.
* **Dual-Wave Synthesis:** Blends two distinct wave algorithms (Type A and Type B) with adjustable mixing to create hybrid geometries.
* **25+ Wave Algorithms:** Includes classic Cartesian/Radial waves, plus exotic types like Gyroid, Quasicrystal, Voronoi, Perlin Noise, and Bessel functions.
* **Physics Simulation:** Particles exhibit Brownian motion (jitter), drag, and momentum, giving the visual an organic, fluid-like behavior.
* **High-Res Export:** Export high-resolution (up to 8k) PNG snapshots with optional transparency.
* **Tweakpane Integration:** Full GUI control over physics, wave modes, colors, and camera settings.

## 🚀 Installation & Usage
* [Three.js](https://threejs.org/)
* [Tweakpane](https://tweakpane.github.io/docs/)
### Physics & Performance

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `particleCount` | `int` | `400,000` | The total number of particles in the simulation. Higher numbers create denser lines but require more GPU power. |
| `gridSize` | `int` | `512` | The resolution of the underlying energy grid. Higher values result in smoother curves but slower recalculations when wave settings change. |
| `settleStrength` | `float` | `1.5` | How strongly particles are pulled toward nodal lines. Higher values make lines form instantly; lower values create a "drifting" effect. |
| `jitter` | `float` | `0.1` | Adds random Brownian motion to particles. Keeps them "alive" and prevents them from getting stuck in local minima. |
| `drag` | `float` | `0.85` | Friction applied to particle velocity. Determines how quickly particles lose momentum. |

### Wave & Math Generation

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `waveTypeA/B` | `string` | `Cartesian` | The mathematical algorithms used to generate the field (e.g., "Radial", "Hexagonal", "Voronoi"). |
| `waveMix` | `float` | `0.5` | Blends between Wave Type A (0.0) and Wave Type B (1.0). |
| `modeCount` | `int` | `4` | The number of overlapping wave layers (harmonics) used to generate the interference pattern. |
| `mMin` / `mMax` | `float` | `2` - `12` | The frequency range for the primary wave axis. Higher values create more complex, tighter patterns. |
| `nMin` / `nMax` | `float` | `2` - `12` | The frequency range for the secondary wave axis. |

### Visuals & Camera

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `color` | `hex` | `#72dcff` | The tint color of the particles. |
| `viewScale` | `int` | `800` | The coordinate bounds of the world. Increases the "zoom out" factor of the simulation space. |
| `exportBaseSize` | `int` | `2048` | The base width (in pixels) for image export. |

## 🕹️ Controls (GUI)

A Tweakpane UI is generated automatically on launch with the following folders:

1. **Particles:** Adjust count, size, and color.
2. **Motion:** Fine-tune the physics (speed, drag, jitter).
3. **Waves:** Select algorithms (Type A/B) and mix ratio.
4. **Modes:** Control the complexity (frequency ranges) and randomness of the wave layers.
5. **Capture:** Camera controls (Pan/Zoom) and **Save Image** button.
4. **Particle Update Loop:**
* Particles check their position against the grid.
* Bilinear interpolation retrieves the precise gradient at that location.
* Velocity is updated to move the particle *down* the gradient (towards 0 energy).
* Jitter and Drag are applied to simulate fluid dynamics.
