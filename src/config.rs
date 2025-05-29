pub static CONFIG: Config = Config {
    scene: Scene {
        window_size: 1300 as f32,
        resolution: 100 as f32,
        frame_rate: 60.0,
        time_step: 1.0 / 4.0,
        is_grid_visible: true,
        are_particles_visible: true,
        particle_iterations: 2,
        pressure_iterations: 50,
    },
    particle: Particle {
        total: 1000, // Aguanta las 10.000 bastante bien.
        radius: 5.0 , // Aguanta 2, pero puede que sea mejor 3.
    },
    environment: Environment {
        gravity: -9.81,
        viscosity: 0.0,
        density: 1000.0,
        flip_ratio: 0.9,
        over_relaxation: 1.9,
        obstacle_radius: 0.0,
        compensate_drift: true,
        color_diffusion_coeff: 0.001,
    },
};

pub struct Config {
    pub particle: Particle,
    pub environment: Environment,
    pub scene: Scene,
}

pub struct Particle {
    pub total: usize,
    pub radius: f32,
}

pub struct Environment {
    pub gravity: f32,
    pub viscosity: f32,
    pub density: f32,
    pub flip_ratio: f32,
    pub over_relaxation: f32,
    pub obstacle_radius: f32,
    pub compensate_drift: bool,
    pub color_diffusion_coeff: f32,
}

pub struct Scene {
    pub window_size: f32,
    pub resolution: f32,
    // pub height: f32,
    pub time_step: f32,
    pub are_particles_visible: bool,
    pub is_grid_visible: bool,
    pub frame_rate: f32,
    pub particle_iterations: u32,
    pub pressure_iterations: u32,
}
