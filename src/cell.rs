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
    pub color: Color,
    pub cell_type: CellType,
    pub density: f32,
    pub fluid_coefficient: f32,
    pub pressure: f32,
    pub delta: Velocity,
}

impl Cell {
    pub fn new(position: Position, color: Color, cell_type: CellType) -> Self {
        Cell { 
            position: Position::fromPosition(position), 
            velocity: Velocity::new(0.0, 0.0),
            previous_velocity: Velocity::new(0.0, 0.0),
            color: Color::from_color(color),
            cell_type: cell_type,
            density: 0.0,
            fluid_coefficient: 0.0,
            pressure: 0.0,
            delta: Velocity::new(0.0, 0.0),
        }
    }
}