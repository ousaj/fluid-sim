pub struct Utils {}

impl Utils {
    pub fn get_pixel(index: usize, width: f32) -> f32 {
        (index as f32) * width
    }

    pub fn is_border(row: usize, column: usize, resolution: usize) -> bool {
        let is_left: bool = column == 0;
        let is_right: bool = column == resolution - 1;
        let is_bottom: bool = row == 0;
        let is_top: bool = row == resolution - 1;
        
        is_left || is_right || is_bottom || is_top
    }

    pub fn clamp(x: f32, min: f32, max: f32) -> f32 {
        if x < min {
            min
        } else if x > max {
            max 
        } else {
            x 
        }
    }

    pub fn clamp_usize(value: usize, min: usize, max: usize) -> usize {
        if value < min {
            min
        } else if value > max {
            max
        } else {
            value
        }
    }
}
