// Screen.
pub const WINDOW_SIZE: f32 = 1000.0;
pub const RESOLUTION: usize = 100;
pub const FRAME_RATE: f32 = 60.0;
pub const CELL_SIZE: f32 = WINDOW_SIZE / RESOLUTION as f32;

// Particles.
pub const TOTAL_PARTICLES: usize = RESOLUTION * RESOLUTION;
pub const PARTICLE_RADIUS: f32 = CELL_SIZE * 0.4;
pub const PARTICLE_COLOR: [f32; 3] = [0.6, 0.6, 1.0];
pub const PARTICLE_VELOCITY_RANGE: [f32; 2] = [0.0, 10.0];

// Environment.
pub const TIME_STEP: f32 = 1.0 / 5.0;
pub const GRAVITY: f32 = -9.81 * 5.0;
pub const GLOW_MULTIPLIER: f32 = 6.0;
pub const BLUR_MULTIPLIER: f32 = 16.0;
pub const PARTICLE_ITERATIONS: usize = 2;
pub const PRESSURE_ITERATIONS: usize = 50;
pub const VISCOSITY: f32 = 0.0;
pub const DENSITY: f32 = 1000.0;
pub const FLIP_RATIO: f32 = 0.9;
pub const OVER_RELAXATION: f32 = 1.9;
pub const CELL_DENSITY_RANGE: [f32; 2] = [0.0, 3.0];

// Visibility.
pub const ARE_PARTICLES_BLURRED: bool = true;
pub const IS_VELOCITY_MAPPED: bool = false;
pub const IS_GRID_VISIBLE: bool = false;
pub const ARE_PARTICLES_VISIBLE: bool = true;

// Obstacle.
pub const OBSTACLE_RADIUS: f32 = 50.0;

pub const TOTAL_CELLS: usize = RESOLUTION * RESOLUTION as usize;
pub const PARTICLE_DIAMETER: f32 = PARTICLE_RADIUS * 2.0;