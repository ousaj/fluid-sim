#[derive(Copy, Clone)]
pub struct Velocity {
    pub u: f32,
    pub v: f32,
}

impl Velocity {
    pub fn new(u: f32, v: f32) -> Self {
        Velocity {
            u,
            v            
        }
    }

    pub fn fromVelocity(velocity: Velocity) -> Self {
        Velocity {
            u: velocity.u,
            v: velocity.v
        } 
    }
}