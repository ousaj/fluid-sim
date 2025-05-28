use std::time::{Duration, Instant};

use scene::Scene;
use sdl2::event::Event;

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
    let mut fps: u32;
    let mut frames: u32 = 0;
    let mut last_time = Instant::now();
    let frame_duration = Duration::from_secs_f32(1.0 / scene.frame_rate);

    'running: loop {
        let frame_start = Instant::now();
        
        for event in scene.windsl.event_pump.poll_iter() {
            match event {
                Event::Quit { .. } => break 'running,
                _ => { }
            }
        }

        frames += 1;
        if last_time.elapsed() >= Duration::from_secs(1) {
            fps = frames;
            frames = 0;
            last_time = Instant::now();
            println!("FPS: {}", fps);
        }

        scene.update();

        let elapsed = frame_start.elapsed();
        if elapsed < frame_duration {
            std::thread::sleep(frame_duration - elapsed);
        }
    }
}
