# Chladni Particle Flow

Generative art simulation built with Three.js. This project simulates Chladni figures (nodal patterns formed on vibrating plates) using a particle system that reacts to a computed energy field. Instead of drawing explicit 2D contour lines, it generates an interference “terrain” from blended wave functions. Particles flow downhill (via gradient descent) and settle into low-energy valleys (nodes).
## [> Test it out](https://octaviaan.github.io/Chladni-Particles/)

## Key Features

- **Particle System:** Renders a particle field using `THREE.Points` backed by low-level `Float32Array` buffers.
- **Dual-Wave Synthesis:** Blends two distinct wave algorithms (**Type A** and **Type B**) using a global mix (`waveMix`) plus a spatially varying mask for hybrid geometries.
- **25+ Wave Algorithms:** Includes classic Cartesian/Radial waves plus stylized variants like Gyroid, Quasicrystal, Voronoi, Hexagonal, Lattice, RingInterference, Wallpaper, and more.
- **Physics Simulation:** Particles include Brownian motion (jitter), drag, momentum, and a speed clamp for stable, organic motion.
- **High-Res Export:** Save PNG snapshots with optional transparent background.
- **Tweakpane Integration:** GUI control over particles, motion, waves, modes, camera, and capture.


## Examples

![fav](https://github.com/user-attachments/assets/e28b23bc-3a86-4206-b639-8119bf60b763)


## Installation & Usage

- Three.js
- Tweakpane
Run with any local server (recommended) and open `index.html`.

## Physics & Performance

| Parameter        | Type    | Default  | Description                                                                                                   |
| ---------------- | ------- | -------- | ------------------------------------------------------------------------------------------------------------- |
| `particleCount`  | `int`   | `300000` | Total number of particles. Higher values create denser patterns but require more CPU/GPU bandwidth.           |
| `gridSize`       | `int`   | `256`    | Resolution of the underlying energy grid. Higher values produce smoother gradients but slower field rebuilds. |
| `settleStrength` | `float` | `2.0`    | Strength of the gradient-descent pull (how quickly particles lock into valleys).                              |
| `jitter`         | `float` | `0.1`    | Random Brownian motion added each frame to keep motion alive and avoid local sticking.                        |
| `drag`           | `float` | `0.85`   | Velocity damping (friction). Higher = faster settling, lower = more drift.                                    |
| `speedLimit`     | `float` | `2.0`    | Caps particle speed for stability and cleaner lines.                                                          |

## Wave & Math Generation

| Parameter                 | Type      | Default                | Description                                                                       |
| ------------------------- | --------- | ---------------------- | --------------------------------------------------------------------------------- |
| `waveTypeA` / `waveTypeB` | `string`  | `Cartesian` / `Radial` | Wave algorithms used to generate the field.                                       |
| `waveMix`                 | `float`   | `0.5`                  | Global blend between Type A (0.0) and Type B (1.0).                               |
| `fieldScale`              | `float`   | `1.0`                  | Scales the coordinate domain used when evaluating wave functions.                 |
| `modeCount`               | `int`     | `4`                    | Number of overlapping mode layers used in field construction.                     |
| `mRange` / `nRange`       | `object`  | `{2→4}` / `{4→8}`      | Frequency ranges used to generate per-mode parameters (`m`, `n`).                 |
| `integerModes`            | `boolean` | `true`                 | Quantizes `m`/`n` to integers (classic mode stepping).                            |

## Visuals & Camera

| Parameter        | Type    | Default   | Description                                               |
| ---------------- | ------- | --------- | --------------------------------------------------------- |
| `color`          | `hex`   | `#b3a79b` | Particle color tint.                                      |
| `particleSize`   | `float` | `2.0`     | Point size (screen-space).                                |
| `exportBaseSize` | `int`   | `3000`    | Base export width in pixels (height derived from aspect). |
| `rectAspect`     | `float` | `1.78`    | Export aspect ratio (1 = square, 1.78 ≈ landscape).       |
| `exportOpaque`   | `bool`  | `false`   | When false, export uses a transparent background.         |

## Controls (GUI)

A Tweakpane UI is generated automatically on launch with the following folders:

1. **Particles:** Adjust count, size, and color (count change rebuilds buffers).
2. **Motion:** Fine-tune physics (settle strength, drag, jitter, speed limit).
3. **Waves:** Select algorithms (Type A/B), mix ratio, noise, and field scale.
4. **Modes:** Control complexity (`modeCount`, `mRange`, `nRange`, integer quantization).
5. **Capture:** Aspect, export size, background transparency, camera pan/zoom, and export buttons.

### Mouse & Keyboard

- **Mouse wheel:** Zoom
- **Mouse drag:** Pan
- **Double click:** Reset view
- **R key:** Randomize modes and waves

## Wave Functions

![wave-functions](https://github.com/user-attachments/assets/cf0a15f9-650d-4795-9d3c-2c4fca027ce9)


## Wave Mixing Examples

![Frame 119](https://github.com/user-attachments/assets/d5492f8f-a354-45d2-8fdc-333ddbb8056d)
