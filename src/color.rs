use crate::config::PARTICLE_COLOR;

#[derive(Copy, Clone)]
pub struct Color {
    pub r: f32,
    pub g: f32,
    pub b: f32,
}

impl Color {
    pub fn new(r: f32, g: f32, b: f32) -> Self {
        Color {
            r,
            g,
            b
        }
    }
    
    pub fn from_color(color: Color) -> Self {
        Color {
            r: color.r,
            g: color.g,
            b: color.b
        }
    }

    pub fn border_color() -> Self {
        Color {
            r: 0.3,
            g: 0.3,
            b: 0.3,
        }
    }

    pub fn background_color() -> Self {
        Color {
            r: 0.1,
            g: 0.1,
            b: 0.1,
        }
    }

    pub fn particle_color() -> Self {
        Color {
            r: PARTICLE_COLOR[0],
            g: PARTICLE_COLOR[1],
            b: PARTICLE_COLOR[2],
        }
    }
}