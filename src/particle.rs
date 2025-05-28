use crate::position::Position;
use crate::color::Color;
use crate::velocity::Velocity;

#[derive(Copy, Clone)]
pub struct Particle {
    pub position: Position,
    pub color: Color,
    pub velocity: Velocity,
    pub density: f32,
}

impl Particle {
    pub fn new(position: Position, color: Color) -> Self {
        Particle {
            position: Position::fromPosition(position),
            color: Color::from_color(color),
            velocity: Velocity::new(0.0, 0.0),
            density: 0.0
        }
    }
}