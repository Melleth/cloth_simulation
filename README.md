# Rusty cloth simulation.
Repository for the small project of the INFOMGP course, year '18-'19. All code by Gert Meijer (6026087) and Siemen Kraayenbrink (4032322)
![screenshot died, appearently](https://github.com/Melleth/cloth_simulation/blob/master/screenshot.png "SCREENSHOT :D")
## Features
1. Mass spring damper model simulation.
2. Wind simulation by calculating the wind direction versus the mesh triangles normal.
3. Self intersection resolvement through sphere colision along the vertices, optimized with a grid.
4. Multithreaded spring force calculation.

## Controls
Move camera with arrow keys, mouse click and drag to pan camera, scroll to zoom.
W - Toggle wireframe rendering of cloth.
Space - Toggle simulation play/pause.

## Build instructions.
1. This project was developed with the rust programming language. You can install all things rust very easily through rustup,
which can be found here: https://rustup.rs/
2. Clone the repo
3. cd into repo
4. Build and run the project by entering ``cargo run --release``

## Dependencies
We make use of a couple dependencies. They are listed in Cargo.toml, which is used by the build system to download and
compile the required source.

1. kiss3d: https://crates.io/crates/kiss3d
.. Amazingly easy graphics lib that allows very easy window creating, input management and provides a simple way to render 3D graphics.
2. Nalgebra: https://crates.io/crates/nalgebra
.. linear algebra lib for rust.
3. Rayon: https://crates.io/crates/rayon
.. Data parallelism lib that implements work stealing.
