use std::fs::File;
use std::io::Write;
use std::time::{Duration, Instant};
use uuid::Uuid;

use fluid::Fluid;
use scene::Scene;
use sdl2::{ event::Event, keyboard::Keycode, mouse::MouseButton };

use config::WINDOW_SIZE;
use config::RESOLUTION;
use config::GLOW_MULTIPLIER;
use config::BLUR_MULTIPLIER;
use config::TIME_STEP;
use config::PARTICLE_ITERATIONS;
use config::PRESSURE_ITERATIONS;
use config::FRAME_RATE;
use config::GRAVITY;
use config::VISCOSITY;
use config::DENSITY;
use config::FLIP_RATIO;
use config::OVER_RELAXATION;
use config::TOTAL_PARTICLES;
use config::PARTICLE_RADIUS;
use config::PARTICLE_COLOR;
use config::OBSTACLE_RADIUS;
use config::TOTAL_CELLS;
use config::CELL_SIZE;

mod scene;
mod objects;
mod windsl;
mod fluid;
mod utils;
mod position;
mod color;
mod velocity;
mod cell;
mod particle;
mod config;
mod obstacle;

fn main() {
    let mut scene: Scene = Scene::new();
    
    let mut fps: u64;
    let mut frames: u64 = 0;
    let mut fps_sum: u64 = 0;
    let simulation_start: Instant = Instant::now();
    let mut last_time = Instant::now();
    let frame_duration = Duration::from_secs_f32(1.0 / scene.frame_rate);

    'running: loop {
        let frame_start = Instant::now();

        frames += 1;
        if last_time.elapsed() >= Duration::from_secs(1) {
            fps = frames;
            fps_sum += fps;
            frames = 0;
            last_time = Instant::now();

            println!("FPS: {}", fps);
        }
        
        for event in scene.windsl.event_pump.poll_iter() {
            match event {
                Event::Quit { .. } => break 'running,
                Event::KeyDown { keycode: Some(Keycode::G), repeat: false, .. } => {
                    scene.is_grid_visible = !scene.is_grid_visible;
                },
                Event::KeyDown { keycode: Some(Keycode::P), repeat: false, .. } => {
                    scene.are_particles_visible = !scene.are_particles_visible;
                },
                Event::KeyDown { keycode: Some(Keycode::S), repeat: false, .. } => {
                    scene.are_particles_blurred = !scene.are_particles_blurred;
                },
                Event::KeyDown { keycode: Some(Keycode::R), repeat: false, .. } => {
                    scene.fluid = Fluid::new();
                },
                Event::KeyDown { keycode: Some(Keycode::V), repeat: false, .. } => {
                    scene.fluid.is_velocity_mapped = !scene.fluid.is_velocity_mapped;
                },
                Event::MouseButtonDown { mouse_btn, .. } => {
                    match mouse_btn {
                        MouseButton::Left => {
                            scene.is_mouse_dragging = true;
                            scene.obstacle.add_obstacle();
                        },
                        _ => {}
                    };
                },
                Event::MouseButtonUp { mouse_btn, .. } => {
                    match mouse_btn {
                        MouseButton::Left => {
                            scene.is_mouse_dragging = false;
                            scene.obstacle.remove_obstacle();
                        },
                        _ => {}
                    };
                },
                Event::MouseMotion { x, y, xrel, yrel, .. } if scene.is_mouse_dragging => {
                    scene.obstacle.update_position(x, y, xrel, yrel);
                },
                _ => { }
            };
        }

        scene.update();

        let elapsed = frame_start.elapsed();
        if elapsed < frame_duration {
            std::thread::sleep(frame_duration - elapsed);
        }
    }

    if simulation_start.elapsed().as_secs() > 0 {
        let err: &str = "Error";
        let path: String = format!("./benchmarks/{}.txt", Uuid::new_v4().to_string());
        let mut f = File::create(path).expect(&err);

        writeln!(f, "Seconds elapsed: {}", simulation_start.elapsed().as_secs())
            .expect(&err);
        writeln!(f, "Average FPS: {}", fps_sum / simulation_start.elapsed().as_secs())
            .expect(&err);
        writeln!(f, "Window size: {}", WINDOW_SIZE)
            .expect(&err);
        writeln!(f, "Resolution: {}", RESOLUTION)
            .expect(&err);
        writeln!(f, "Cell size: {}", CELL_SIZE)
            .expect(&err);
        writeln!(f, "Particle radius: {}", PARTICLE_RADIUS)
            .expect(&err);
        writeln!(f, "Total particles: {}", TOTAL_PARTICLES)
            .expect(&err);
        writeln!(f, "Total cells: {}", TOTAL_CELLS)
            .expect(&err);
        writeln!(f, "Grow multiplier: {}", GLOW_MULTIPLIER)
            .expect(&err);
        writeln!(f, "Blur multiplier: {}", BLUR_MULTIPLIER)
            .expect(&err);
        writeln!(f, "Time step: {}", TIME_STEP)
            .expect(&err);
        writeln!(f, "Particle iterations: {}", PARTICLE_ITERATIONS)
            .expect(&err);
        writeln!(f, "Pressure iterations: {}", PRESSURE_ITERATIONS)
            .expect(&err);
        writeln!(f, "Frame rate : {}", FRAME_RATE)
            .expect(&err);
        writeln!(f, "Gravity: {}", GRAVITY)
            .expect(&err);
        writeln!(f, "Viscosity: {}", VISCOSITY)
            .expect(&err);
        writeln!(f, "Density: {}", DENSITY)
            .expect(&err);
        writeln!(f, "FLIP ratio: {}", FLIP_RATIO)
            .expect(&err);
        writeln!(f, "PIC ratio: {}", 1.0 - FLIP_RATIO)
            .expect(&err);
        writeln!(f, "Over-relaxation: {}", OVER_RELAXATION)
            .expect(&err);
        writeln!(f, "Total particles: {}", TOTAL_PARTICLES)
            .expect(&err);
        writeln!(f, "Particle radius: {}", PARTICLE_RADIUS)
            .expect(&err);
        writeln!(f, "Particle color: {:?}", PARTICLE_COLOR)
            .expect(&err);
        writeln!(f, "Obstacle radius: {}", OBSTACLE_RADIUS)
            .expect(&err);
        writeln!(f, "Is grid visible: {}", scene.is_grid_visible)
            .expect(&err);
        writeln!(f, "Are particles visible: {}", scene.are_particles_visible)
            .expect(&err);
        writeln!(f, "Are particles blurred: {}", scene.are_particles_blurred)
            .expect(&err);
    }
}
