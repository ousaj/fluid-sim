use crate::position::Position;
use crate::velocity::Velocity;
use crate::config::WINDOW_SIZE;

pub struct Obstacle {
    pub position: Position,
    pub velocity: Velocity,
    pub is_enabled: bool,
}

impl Obstacle {
    pub fn new(position: Position, velocity: Velocity) -> Self {
        Obstacle {
            position,
            velocity,
            is_enabled: false,
        }
    }

    pub fn update_position(&mut self, x: i32, y: i32, xrel: i32, yrel: i32) {
        self.position.x = x as f32;
        self.position.y = WINDOW_SIZE - y as f32;
        self.velocity.u = (xrel * 100) as f32;
        self.velocity.v = (-yrel * 100) as f32;
    }

    pub fn add_obstacle(&mut self) {
        self.is_enabled = true;
    }

    pub fn remove_obstacle(&mut self) {
        self.is_enabled = false;
    }
}
