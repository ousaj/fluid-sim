use crate::position::Position;
use crate::velocity::Velocity;

pub struct Obstacle {
    pub position: Position,
    pub radius: f32,
    pub velocity: Velocity,
}

impl Obstacle {
    pub fn new(position: Position, radius: f32, velocity: Velocity) -> Self {
        Obstacle {
            position,
            radius,
            velocity,
        }
    }
}
