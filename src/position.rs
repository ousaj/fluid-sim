#[derive(Copy, Clone)]
pub struct Position {
    pub x: f32,
    pub y: f32,
}

impl Position {
    pub fn new(x: f32, y: f32) -> Self {
        Position {
            x,
            y
        }
    }
    
    pub fn fromPosition(position: Position) -> Self {
        Position {
            x: position.x,
            y: position.y
        }
    }
}