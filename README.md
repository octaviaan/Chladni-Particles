# Chladni Particle Flow

Generative art simulation built with Three.js. This project simulates Chladni figures (nodal patterns formed on vibrating plates) using a particle system that reacts to a computed energy field. Instead of drawing explicit 2D contour lines, it generates an interference “terrain” from blended wave functions. Particles flow downhill (via gradient descent) and settle into low-energy valleys (nodes).
## [> Test it out](https://octaviaan.github.io/Chladni-Particles/)

## Key Features

- **Particle System:** Renders a particle field using `THREE.Points` backed by low-level `Float32Array` buffers.
- **Dual-Wave Synthesis:** Blends two distinct wave algorithms (**Type A** and **Type B**) using a global mix (`waveMix`) for hybrid geometries.
- **25+ Wave Algorithms:** Includes classic Cartesian/Radial waves plus stylized variants like Gyroid, Quasicrystal, Voronoi, Hexagonal, Lattice, RingInterference, Wallpaper, and more.
- **Physics Simulation:** Particles include Brownian motion (jitter), drag, momentum, and a speed clamp for stable, organic motion.
- **Logo Integration Mode:** Upload a logo (transparent PNG recommended) and map particles to fill or outline targets with alpha/luminance extraction.
- **Text Layer Mode:** Render text as its own particle target layer using **BDO Grotesk** (bold), with independent text particle controls.
- **Export Tools:** Save PNG snapshots (full frame or current view), plus loop exports as GIF or WebM.
- **Tweakpane Integration:** GUI control over particles, motion, waves, modes, camera, and capture.

## Recent Updates

- Added two-column controls layout (`Controls` + `Logo & Text`).
- Added logo mask mode switch (`Auto` / `Alpha` / `Luminance`) and threshold workflow improvements.
- Added loop export controls for GIF (`Duration`, `GIF fps`, `GIF scale`, `GIF compression`) and `Save WebM`.
- Updated capture naming from `exportBaseSize` to **Image resolution** in the GUI.
- Simplified controls by removing wave-only randomization and respawn behavior.


## Examples

![fav](https://github.com/user-attachments/assets/e28b23bc-3a86-4206-b639-8119bf60b763)


## Installation & Usage

- Three.js
- Tweakpane
Run with any local server (recommended) and open `index.html`.

## Physics & Performance

| Parameter        | Type    | Default  | Description                                                                                                   |
| ---------------- | ------- | -------- | ------------------------------------------------------------------------------------------------------------- |
| `particleCount`  | `int`   | `500000` | Total number of particles. Higher values create denser patterns but require more CPU/GPU bandwidth.           |
| `gridSize`       | `int`   | `256`    | Resolution of the underlying energy grid. Higher values produce smoother gradients but slower field rebuilds. |
| `settleStrength` | `float` | `5.0`    | Strength of the gradient-descent pull (how quickly particles lock into valleys).                              |
| `jitter`         | `float` | `0.07`   | Random Brownian motion added each frame to keep motion alive and avoid local sticking.                        |
| `drag`           | `float` | `0.85`   | Velocity damping (friction). Higher = faster settling, lower = more drift.                                    |
| `speedLimit`     | `float` | `2.0`    | Caps particle speed for stability and cleaner lines.                                                          |

## Wave & Math Generation

| Parameter                 | Type      | Default                | Description                                                                       |
| ------------------------- | --------- | ---------------------- | --------------------------------------------------------------------------------- |
| `waveTypeA` / `waveTypeB` | `string`  | `Cartesian` / `Gyroid` | Wave algorithms used to generate the field.                                       |
| `waveMix`                 | `float`   | `0.5`                  | Global blend between Type A (0.0) and Type B (1.0).                               |
| `fieldScale`              | `float`   | `0.9`                  | Scales the coordinate domain used when evaluating wave functions.                 |
| `modeCount`               | `int`     | `6`                    | Number of overlapping mode layers used in field construction.                     |
| `mRange` / `nRange`       | `object`  | `{2→6}` / `{4→10}`     | Frequency ranges used to generate per-mode parameters (`m`, `n`).                 |
| `integerModes`            | `boolean` | `true`                 | Quantizes `m`/`n` to integers (classic mode stepping).                            |

## Visuals & Camera

| Parameter        | Type    | Default   | Description                                                                             |
| ---------------- | ------- | --------- | --------------------------------------------------------------------------------------- |
| `color`          | `hex`   | `#4f4030` | Particle color tint.                                                                    |
| `backgroundColor` | `hex`   | `#1C1A17` | Render background color used in viewport and opaque exports.                            |
| `particleSize`   | `float` | `2.0`     | Point size (screen-space).                                                              |
| `exportBaseSize` | `int`   | `2048`    | Base width for **Save image** export (GUI label: `Image resolution`).                  |
| `rectAspect`     | `float` | `1.78`    | Internal landscape aspect used for full-frame export framing (not exposed in the GUI). |

## Controls (GUI)

A Tweakpane UI is generated automatically on launch with the following folders:

1. **Particles:** Adjust count, size, and color (count change rebuilds buffers).
2. **Motion:** Fine-tune physics (settle strength, drag, jitter, speed limit).
3. **Waves:** Select algorithms (`Wave A`/`Wave B`), blend (`Mix`), and `Field scale`.
4. **Modes:** Control complexity (`modeCount`, `mRange`, `nRange`) and use `Randomize All`.
5. **Logo:** Upload/clear logo, tune `Mask`, `Threshold`, `Style`, `Invert`, `Image mix`, plus logo particle color/count/size/opacity and offsets/scale.
6. **Text:** Enter text and control `Align`, `X offset`, `Particle count`, `Particle size`, `Text opacity`, `Text color`, and `Font size`.
7. **Capture:** `Image resolution`, camera controls, PNG save buttons, GIF settings, and `Save GIF`/`Save WebM`.

### Mouse & Keyboard

- **Mouse wheel:** Zoom
- **Mouse drag:** Pan
- **Double click:** Reset view
- **R key:** Randomize all
- **S key:** Save image
- **C key:** Save current view
- **V key:** Save GIF loop
- **H key:** Toggle controls

## Wave Functions

![wave-functions](https://github.com/user-attachments/assets/cf0a15f9-650d-4795-9d3c-2c4fca027ce9)


## Wave Mixing Examples

![Frame 119](https://github.com/user-attachments/assets/d5492f8f-a354-45d2-8fdc-333ddbb8056d)
