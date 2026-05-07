<p align="center">
  <video src="turtle_video.MP4" width="640" autoplay loop muted playsinline controls>
    Your browser does not support embedded video. <a href="turtle_video.MP4">Download turtle_video.MP4</a>.
  </video>
</p>

<h1 align="center">Going with the Flow</h1>

<p align="center">
  <strong>A C++ implementation of "Going with the Flow" fluid-structure interaction for swimming creatures.</strong>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/Language-C%2B%2B20-blue" alt="C++20"/>
  <img src="https://img.shields.io/badge/Build-CMake-064F8C" alt="CMake"/>
  <img src="https://img.shields.io/badge/Math-Eigen-brightgreen" alt="Eigen"/>
  <img src="https://img.shields.io/badge/UI-Qt6%20%2B%20OpenGL-41cd52" alt="Qt6 + OpenGL"/>
  <img src="https://img.shields.io/badge/Paper-SIGGRAPH%202024-orange" alt="SIGGRAPH 2024"/>
</p>

---

A from-scratch C++ implementation of Soliman, Padilla, Gross, Knöppel, Pinkall, and Schröder's *Going with the Flow*, which appeared at SIGGRAPH 2024. The simulation drives a deforming surface mesh through a surrounding fluid by exchanging momentum directly with the shape change. There is no fluid grid, no particle solver, and no pressure projection. A swimming turtle or worm is animated by replaying captured pose sequences, and the body's rigid trajectory through the world emerges from the variational integrator alone.

Beyond the paper, the project adds a **piecewise added-mass** extension. The body is segmented into anatomically meaningful regions for the head, the four flippers, and the shell, and each region carries its own fluid-coupling coefficient. This produces slightly different swimming dynamics than the global formulation.

---

## Table of Contents

- [Authors](#authors)
- [Features](#features)
- [Architecture](#architecture)
- [Tech Stack](#tech-stack)
- [Project Structure](#project-structure)
- [Getting Started](#getting-started)
- [Run Modes](#run-modes)
- [Algorithm Reference](#algorithm-reference)
- [Paper](#paper)

---

## Authors

The Turtle Club: Gavin Dhanda, Ryan Aubrey, and Briana Fedkiw.

Slides: [Going with the Flow presentation](https://docs.google.com/presentation/d/1YBkFX76Jy7dQEbvgLU_Urhl9qnAMPBsrKDC2PiWQLVo/edit?usp=sharing).

---

## Features

### Simulation Core

| Feature | Description |
|:--|:--|
| **Variational Integrator** | Full implementation of the paper's Algorithm 1 — Kirchhoff tensor assembly, body and fluid momentum, Newton solve with finite-difference Jacobian, Cayley-map state reconstruction on SE(3). |
| **Body Inertia** | Time-varying 6×6 inertia matrix `M_b` reassembled every substep from the deforming vertex cloud. |
| **Added Mass** | Surface-integral added-mass tensor `M_f`, parameterized by a single shape coefficient δ derived from total area and mean bending. |
| **Lift and Drag** | Per-face flat-plate force (Proposition 1) plus gravity, computed against the relative body-fluid velocity. |
| **Pose Replay** | Loads a sequence of OBJ frames and interpolates linearly between consecutive poses across configurable substeps, for any number of strokes. |

### Piecewise Added-Mass Extension

| Feature | Description |
|:--|:--|
| **Mesh Segmentation** | Automatic partitioning of the body into six regions (head, left/right front flipper, left/right rear flipper, shell) using face-adjacency Dijkstra distances, quantile-based tip detection, and PCA-derived front/side axes. |
| **Per-Region δ** | Each region's δ is computed independently from its own area and bending edges, then applied face-by-face in both the added-mass tensor and the fluid momentum. |
| **Cached Labels** | Partition labels are persisted to `turtle_poses/frame_20_partition_labels.txt` so the segmentation cost is paid once. |
| **Side-by-Side Compare** | A built-in `compare` mode runs both global and piecewise simulations and prints a centroid-trajectory diff. |

### Visualization

| Feature | Description |
|:--|:--|
| **Qt6 + OpenGL Viewer** | Real-time playback of the simulated mesh with texture mapping, free-fly camera, ground plane, and play/pause controls. |
| **Partition Editor** | Interactive tool for inspecting and editing the mesh partition on the reference frame. |
| **Debug Widgets** | Standalone widgets for visualizing body inertia axes, fluid momentum vectors, and partition coloring. |

---

## Architecture

```
+----------------------------------------------------------+
|                         main.cpp                         |
|        argv => { "" | "piecewise" | "compare" |          |
|                  "partition" }                           |
+--------------------------+-------------------------------+
                           |
              +------------+-------------+
              |                          |
              v                          v
    +-------------------+      +-----------------------+
    |    Integrator     |      |  PartitionDebugWidget |
    |                   |      |     (label editor)    |
    |  LoadAllPoses     |      +-----------------------+
    |  Simulate(mode)   |
    +---------+---------+
              |
              | per substep
              v
    +-------------------+      +---------------------+
    |   Algorithm 1     | <--- |     se3.h           |
    |   variational     |      |  cayley_map         |
    |   integrator      |      |  dtau_inv_star      |
    +---------+---------+      +---------------------+
              |
   +----------+-----------+----------------+----------------+
   |          |           |                |                |
   v          v           v                v                v
+--------+ +--------+ +-------------+ +------------+ +-----------+
|inertia | |momentum| |lift_and_drag| | partition  | |mesh_loader|
| M_b M_f| |  l p   | | Prop. 1     | | regions    | |  OBJ+MTL  |
+--------+ +--------+ +-------------+ +------------+ +-----------+
              |
              v
    +-------------------+
    |  output_frames    |
    +---------+---------+
              |
              v
    +-------------------+
    |  SimViewerWidget  |  Qt6 + OpenGL playback
    +-------------------+
```

### Design Decisions

- **No fluid grid.** The fluid never appears as a discretized field. All coupling is integrated against the surface in closed form via Stokes' theorem, following the paper. This is what makes the simulation stable at large timesteps and inexpensive to run.
- **Body-frame state.** Pose is tracked as an `SE(3)` matrix and velocity as a body-frame `se(3)` twist. World-frame meshes are reconstructed only for visualization.
- **Cayley map over `exp`.** SE(3) updates use the Cayley map from Appendix B of the paper. It is rational and avoids the trig of the matrix exponential.
- **Header-only OBJ loader.** [`mesh_loader.h`](src/mesh_loader.h) parses OBJ + MTL with texcoords and fan-triangulates n-gons in a single header, keeping the build graph flat.
- **Mode-dispatched simulation.** Global vs. piecewise added-mass is a single enum at the call site; the integrator branches into either `calc_added_mass` / `calc_fluid_momentum` or their piecewise counterparts without any code duplication in the main loop.

---

## Tech Stack

| Layer | Technology |
|:--|:--|
| **Language** | C++20 |
| **Build** | CMake 3.16+ |
| **Linear Algebra** | [Eigen](https://eigen.tuxfamily.org) (vendored) |
| **Geometry I/O** | Custom OBJ + MTL parser (header-only) |
| **GUI** | Qt6 (Core, Widgets, OpenGL, OpenGLWidgets, Concurrent, Xml, Gui) |
| **Rendering** | OpenGL via GLEW (vendored as a static lib) |
| **Tests** | Pure-C++ test executable (`test_forces`), no Qt dependency |

---

## Project Structure

```
GoingWithTheFlow/
|-- CMakeLists.txt               # Two targets: GoingWithTheFlow + test_forces
|-- README.md
|-- Eigen/                       # Vendored linear algebra library
|-- glew/                        # Vendored OpenGL extension wrangler
|
|-- src/
|   |-- main.cpp                 # Entry point, argv mode dispatch
|   |
|   |-- integrator.h/.cpp        # Algorithm 1 driver: pose load, substep loop, Newton solve
|   |-- momentum.h/.cpp          # Body + fluid momentum (global + piecewise variants)
|   |-- inertia.h/.cpp           # Body inertia M_b and added mass M_f (global + piecewise)
|   |-- lift_and_drag.cpp        # Prop. 1 face force, gravity, total external wrench, mesh volume
|   |-- se3.h                    # Cayley map and dtau_inv_star on SE(3) / se(3)
|   |-- mesh_loader.h            # Header-only OBJ + MTL parser with texture coords
|   |-- partition.h/.cpp         # Mesh segmentation into 6 anatomical regions + label I/O
|   |
|   |-- simviewerwidget.h/.cpp   # Qt6 + OpenGL playback viewer (textured mesh, ground, camera)
|   |-- bodydebugwidget.*        # Debug viz for body momentum
|   |-- fluiddebugwidget.*       # Debug viz for fluid momentum
|   |-- inertiadebugwidget.*     # Debug viz for inertia tensor axes
|   |-- partitiondebugwidget.*   # Interactive partition label editor
|   |
|   |-- graphics/                # Camera, shader, shape, debug rendering helpers
|   |-- test_forces.cpp          # Pure-C++ unit tests for inertia + lift/drag + volume
|
|-- turtle_poses/                # 40 OBJ frames + texture + cached partition labels
|-- worm_poses/                  # 40 OBJ frames for the worm scene
|-- resources/                   # Qt resource bundle (shaders)
|-- util/                        # Misc utilities included on the compiler search path
|-- output/                      # Runtime-created directory for simulation artifacts
|-- test_images/                 # Reference images and a sample output video
|-- project-files/               # Paper PDF + project proposal
```

---

## Getting Started

### Prerequisites

- A C++20 compiler (Apple Clang 15+, GCC 11+, or MSVC 19.30+)
- CMake 3.16 or newer
- Qt6 (Core, Widgets, OpenGL, OpenGLWidgets, Concurrent, Xml, Gui) — required for the viewer
- An OpenGL-capable display

Eigen and GLEW are already vendored in the repository, so no external installs are needed for those.

### Build

From the project root:

```bash
cd /Users/ryana/Documents/Final/GoingWithTheFlow
rm -rf build
cmake -S . -B build
cmake --build build --target GoingWithTheFlow
```

`rm -rf build` is only needed when reconfiguring from scratch. For iterative rebuilds you can skip it and just rerun the `cmake --build` command.

### Run

```bash
./build/GoingWithTheFlow
./build/GoingWithTheFlow piecewise
./build/GoingWithTheFlow compare
```

The first command runs the standard simulation with the global added-mass formulation. The other two are described below.

---

## Run Modes

The executable takes a single optional argument that selects a simulation mode.

| Command | Mode | What it does |
|:--|:--|:--|
| `./build/GoingWithTheFlow` | **Global** (default) | Runs the full simulation with a single δ derived from the whole mesh, then opens a Qt viewer that plays back the resulting trajectory. |
| `./build/GoingWithTheFlow piecewise` | **Piecewise** | Same simulation pipeline, but uses per-region δ values from the partitioned mesh. Opens a single viewer for the piecewise result. |
| `./build/GoingWithTheFlow compare` | **Compare** | Runs both simulations back-to-back, prints centroid-trajectory metrics (mean / max / final difference plus total displacement vectors) to stdout, and opens two viewers side by side for direct visual comparison. |
| `./build/GoingWithTheFlow partition` | **Partition editor** | Skips the simulation and opens an interactive editor for inspecting and editing the mesh partition on the reference frame. |

### Switching the active creature

The integrator currently animates the worm scene. To switch to the turtle, open [`src/integrator.h`](src/integrator.h). The worm scene parameters block is active, and the turtle block is commented out directly above it. Swap which block is commented and rebuild.

---

## Algorithm Reference

| Algorithm | Implemented in | Notes |
|:--|:--|:--|
| **Alg. 1 — Variational integrator** | [`src/integrator.cpp`](src/integrator.cpp) | Substep loop assembles the Kirchhoff tensor and momentum, then solves a 6×6 nonlinear system per substep with a finite-difference Jacobian (10 iterations, residual tolerance `1e-8`). |
| **Alg. 2 — Lift and drag (Prop. 1)** | [`src/lift_and_drag.cpp`](src/lift_and_drag.cpp) | Per-face flat-plate force evaluated against `(rigid + shape-change) - background_flow`. |
| **Alg. 3 — Body inertia `M_b`** | [`src/inertia.cpp`](src/inertia.cpp) | 6×6 symmetric matrix with inertia tensor, first-moment cross-coupling, and total mass. |
| **Body momentum `[l_b; p_b]`** | [`src/momentum.cpp`](src/momentum.cpp) | Per-vertex contribution from `gamma_prime`. |
| **Added mass `M_f`** | [`src/inertia.cpp`](src/inertia.cpp) | Surface integral over faces. Both single-δ and per-face-δ variants. |
| **Fluid momentum `[l_f; p_f]`** | [`src/momentum.cpp`](src/momentum.cpp) | Surface integral against face-centroid velocity. Both single-δ and per-face-δ variants. |
| **SE(3) Cayley map** | [`src/se3.h`](src/se3.h) | `cayley_map` and the inverse-transpose differential `dtau_inv_star` used by the Newton solve. |
| **Mesh partition (extension)** | [`src/partition.cpp`](src/partition.cpp) | Six-region segmentation pipeline plus a bbox-based fallback. |

---

## Paper

Soliman, Y., Padilla, M., Gross, O., Knöppel, F., Pinkall, U., & Schröder, P. (2024). *Going with the Flow.* ACM Transactions on Graphics (SIGGRAPH 2024).

Local copy: [`project-files/GoingWithTheFlow.pdf`](project-files/GoingWithTheFlow.pdf). Hosted copy: [yousufsoliman.com/projects/GoingWithTheFlow.pdf](https://www.yousufsoliman.com/projects/download/GoingWithTheFlow.pdf).
