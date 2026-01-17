use crate::position::Position;
use crate::color::Color;
use crate::velocity::Velocity;

#[derive(Copy, Clone)]
pub struct Particle {
    pub position: Position,
    pub color: Color,
    pub velocity: Velocity,
    pub cell_index: usize,
}

impl Particle {
    pub fn new(position: Position) -> Self {
        Particle {
            position: Position::from_position(position),
            color: Color::particle_color(),
            velocity: Velocity::new(0.0, 0.0),
            cell_index: 0,
        }
    }
}