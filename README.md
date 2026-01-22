# Eulerian fluid simulation

The project consists of an eulerian fluid simulation, using the Rust programming
language and the OpenGL graphics library, and also using the PIC (Particle in Cell)
and FLIP (Fluid Implicit Particle) methodologies. We explore the different phases
the simulation goes through and the data structures required to make a successful execution in a two-dimensional grid environment.

## Preview

https://github.com/user-attachments/assets/52550f8e-b802-42b8-b5bc-dac84823036c

## Technologies used

- Rust
- OpenGL
- GLSL
- Python

## Features

- `R`: Restarts the simulation
- `mousedrag`: Manipulates the fluid and reacts to the acceleration of the mouse drag.
- `P`: Toggles the particle visibility (shows dots as particles).
- `V`: Maps the velocity of each particle (in a blue-orange gradient).
- `S`: Toggles the blur and glow of particles in order to have a more "realistic" fluid,
    by giving the sensation of cohesion between particles.
- `G`: Toggles the grid visibility.

At the end of the simulation, the execution details of the simulation are
stored in the `benchmarks` folder, using a Python script, graphs can be generated
to compare performance under different parameters.

### Setup

The optimal parameters I've found are already set in the `config.rs` file:

- `WINDOW_SIZE = 1000.0`: Defines the window size in pixels.
- `RESOLUTION = 80.0`: Sets the simulation's resolution. In other words,
how many cells there are in a row or column.
- `FRAME_RATE = 60.0`: Sets the frame rate.
- `CELL_SIZE = WINDOW_SIZE / RESOLUTION`: Calculates the cell size, dividing
the window size between the resolution.
- `TOTAL_PARTICLES = RESOLUTION * RESOLUTION`: Defines the total number of
particles rendered in the simulation.
- `PARTICLE_RADIUS = CELL_SIZE * 0.4`: Sets the particle radius, in relation
to the cell size.
- `VISCOSITY = 0.0`: Defines the fluid's viscosity. In other words, the
quantity of friction or movement limitation of the fluid in the medium.
- `GLOW_MULTIPLIER = 6.0`: It's a multiplier that applies to the particle's glow
intensity. This allows us to have a more realistic effect for the fluid, so that
they don't appear round.
- `BLUR_MULTIPLIER = 10.0`: It's a multiplier that affects to the particle's blur,
in order to have a more realistic fluid, so that they don't appear round.
- `FLIP_RATIO = 0.9`: A parameter used in fluid simulations to balance the
FLIP technique, with PIC.
- `OVER_RELAXATION = 1.9`: A parameter used in the simulation's numerical relaxation, that adjusts the speed that the calculations get solved.
- `OBSTACLE_RADIUS = 50.0`: Defines the obstacles radius in the simulation. In our
case, it's the cursor, whenever a left click is triggered inside the window.
 
## The process (How I built it)

In order to build this project, I've taken as a starting point practial
implementations of fluid simulations developed by other autors, which I will link
below.

The methodology consists on implementing a simulation from scratch, structuring
the code and data according to Rust's best practices, and taking advantage of
its performance, concurrency capabilities and code and memory safety.

I've taken a special interest in maintaining a modular and clear architecture,
in order that every part of the system is easily modifiable or expandable in the
future.

For the visual representation, OpenGL has been used as the rendering API and
SDL2 as the window and event management layer. This allows fine-grained control
over the graphics pipeline, making it possible to directly send the vertices of
the polygons we want to render. In this way, the GPU handles rendering and
coordinate calculations, freeing the CPU for physical computations.

Regarding performance, we have employed CPU parallelization, taking advantage
of Rust's concurrency tools through the `Rayon` library. This has allowed us
to divide the processing load of the different simulation steps among the
processor's cores, enabling the handling of a large number of particles without
compromising the frame rate.

Similarly, GPU parallelization has been used for renderign, so we expect the 
simulation's performance to be significantly hihger than that of traditional
solutions.

For the simulation to function correctly, a series of steps defining each system
update cycle has been followed. These steps, executed in order, ensure that
the simulation is coherent and that the particles are represented correctly
throughout the fluid simulation.

Initially, we generate the SDL environment and add an OpenGL context in order to
access its API through our window. After that, we initialize all the data structures
required for the simulation, which will be discussed later in future sections,
with their initial state and limits according to the configurable parameters.

- We integrate the particles into the medium, applying external forces (in this case, gravity).
- We correct particle motion through collision detection between particles and
the environment
- We correct particle motion through particle-particle collision detection using
collision detection using collision detection algorithms.
- We transfer particle velocities to the grid.
- We compute the density in each cell according to the number of particles
contained in it.
- We adjust the grid velocity field to satisfy the zero-divergence condition (incompressible flow).
- We transfer grid velocities back to the particles (PIC/FLIP stage).
- We map the density of the cell in which each particle is located to a color.
- We map the density of each cell to a color.

## What I learned

The reasons as to why I decided to take on this project were both personal and
academic. I've been for a long time atracted to the idea of developing videogames
with the fewer amount of libraries possible, in order to learn the internal
mechanisms of whatever I was programming in that moment. I consider that the
majority of implementations and frameworks abstract a lot that logic and it feels
like magic to the novice, even though it's helpful for the experts in order to build
things fast.

- Developed a fully functional real-time fluid simulation using Rust, OpenGL,
and SDL2.
- Gained experience with real-time simulation concepts and performance constraints.
- Learned Rust from scratch and worked with its ownership, concurrency, and
tooling ecosystem.
- Used CPU parallelization with Rayon and GPU parallelization via OpenGL.
- Implemented multiple simulation features such as:
    - Density and velocity maps.
    - Parameterization of simulation behavior.
    - Visual effects like blur and glow.
    - Performance analysis through logging.
- Improved ability to research and learn from official documentation and
community resources.
- Experienced significant technical groth by tackling unfamiliar technologies.

## How it can be improved

- Overall performance could be improved, especially compared to other
optimized implementations.
- Data structures could be simplified (e.g., using arrays instead of many structs).
- The OpenGL context may introduce unnecesary CPU overhead.
- The initial code structure limited safe and effective parallelization.
- Designing with parallelism in mind from the beginning would improve scalability.
- Achievign full loop parallelization remains an open improvement.
- The simulation could be extended to:
    - A 3D environment.
    - Rotating containers with intertia effects.
- Rust was approached too much from an object-oriented mindset, limiting its
full potential.

## How to run the project

### Rust installation

Below, the steps required to install Rust on the most common operating systems
are described, along with the dependencies needed to compile and run the
simulation developed in this project.

This guide focuses on a development environment based on Rust, OpenGL and SDL2.

#### General requirements

- Access to a terminal.
- Connection to Internet.
- Admin rights to install software.
- A C/C++ compiler, required for some dependencies such as SDL2.

#### Linux/MacOS installation

To install on Linux/MacOS:

1. Open the terminal and execute the following command, to install Rust:

```bash
curl --proto ‘=https’ --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

2. We install the following system dependencies:

```bash
sudo apt install build-essential libgl1-mesa-dev libsdl2-dev
```

3. In order to verify the installation, we execute the following command:

```bash
rustc --version
cargo --version
```

#### Windows installation

1. We download the installer from [https://rustup.rs](https://rustup.rs).
2. We run the installer and select the default installation.
3. To verify the installation, we open a terminal and run the following commands:

```powershell
rustc --version
cargo --version
```

### Compilation of the simulation

This guide explains the steps required to compile and run the fluid dynamics
simulation developed in Rust.

The simulation runs in real time and requires the system dependencies (Rust,
OpenGL, and SDL2) to be correctly installed. Instruction on how to install
these dependencies are provided in previous sections.

Once these requirements are met, the compilations and execution process is
straightforward thanks to the integrated tools offered by Rust, such as Cargo.

Below, we describe the commands and steps required to build the project and
launch the simulation:

1. Navigate to the project directory:

```bash
cd path/to/project
```

2. Compile the project:

```bash
cargo build --release
```

3. Run the simulation:

```bash
cargo run --release
```

### Benchmarking

We're going to need Python in order to execute the `benchmarks.py` file. To install
Python, go to the [following page](https://www.python.org/downloads/).

To check if Python is installed, run the following command in the terminal:

```bash
python --version
```

Generally, `pip` comes installed with Python 3, but you can check it
with the following command:

```bash
pip --version
```

In case it's not installed, we execute the following command:

```bash
python -m ensurepip --upgrade
```

To install `matplotlib` (the library that we'll use to make the plot), we
execute the following command:

```bash
pip install matplotlib seaborn
```

To execute the script and generate the plot, we need to get in the directory 
where the project is located and execute the following command:

```bash
python benchmarks.py
```

This should execute the script and show the bar plot, coming from the log 
directory `benchmarks`.

## References

- Rust Programming Language. (n.d.). https://www.rust-lang.org/
- Introduction - Rust by example. (n.d.).
https://doc.rust-lang.org/stable/rust-by-example/index.html
- std - Rust. (n.d.). https://doc.rust-lang.org/std/
- sdl2 - Rust. (n.d.). https://docs.rs/sdl2/latest/sdl2/
- Rust-Sdl. (n.d.). GitHub - Rust-SDL2/rust-sdl2: SDL2 bindings for Rust.
- GitHub. https://github.com/Rust-SDL2/rust-sdl2
- Codeus. (n.d.). YouTube. https://www.youtube.com/@codedeus
- Coddeus. (n.d.). GitHub - Coddeus/Rust-OpenGL: An OpenGL guide in Rust. GitHub.
https://github.com/Coddeus/Rust-OpenGL
- Ten minute Physics. (n.d.). YouTube. https://www.youtube.com/@TenMinutePhysics
- Matthias-Research. (n.d.). pages/tenMinutePhysics/18-flip.html at master · matthiasresearch/pages. GitHub.
https://github.com/matthias-research/pages/blob/master/tenMinutePhysics/18-flip.html
- Matthias-Research. (n.d.-a). pages/tenMinutePhysics/17-fluidSim.html at master · matthiasresearch/pages. GitHub.
https://github.com/matthias-research/pages/blob/master/tenMinutePhysics/17-fluidSim.html
- Rustsim (2021). Rust for Scientific Computing: Rustsim Libraries. Github Repository.
https://github.com/rustsim
- Berg, M. et al. (2020). Parallel Computing for Fluid Simulations. Springer Handbook of
Computational Fluid Dynamics
- Versteeg, H. K., & Malalasekera, W. (2007). An Introduction to Computational Fluid
Dynamics
