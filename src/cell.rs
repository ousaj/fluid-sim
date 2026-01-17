use crate::position::Position;
use crate::color::Color;
use crate::velocity::Velocity;

#[derive(PartialEq, Copy, Clone)]
pub enum CellType {
    SOLID,
    FLUID,
    AIR
}

#[derive(Copy, Clone)]
pub struct Cell {
    pub position: Position,
    pub velocity: Velocity,
    pub previous_velocity: Velocity,
    pub delta: Velocity,
    pub color: Color,
    pub cell_type: CellType,
    pub solid_coeff: f32,
    pub pressure: f32,
    pub density: f32,
}

impl Cell {
    pub fn new(
        position: Position,
        color: Color,
        cell_type: CellType
    ) -> Self {
        Cell {
            position: Position::from_position(position),
            velocity: Velocity::new(0.0, 0.0),
            previous_velocity: Velocity::new(0.0, 0.0),
            color: Color::from_color(color),
            cell_type: cell_type,
            density: 0.0,
            solid_coeff: 0.0,
            pressure: 0.0,
            delta: Velocity::new(0.0, 0.0),
        }
    }
}