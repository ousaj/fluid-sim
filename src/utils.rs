pub struct Utils {}

impl Utils {
    pub fn normalize_x(x: f32, width: f32) -> f32 {
        2.0 * (x / width) - 1.0
    }

    pub fn normalize_y(y: f32, height: f32) -> f32 {
        1.0 - (y / height) * 2.0
    }

    pub fn get_pixel(index: f32, width: f32) -> f32 {
        index * width
    }

    pub fn is_border(row: u32, column: u32, total_rows: u32, total_columns: u32) -> bool {
        let is_left: bool = column == 0;
        let is_right: bool = column == total_columns - 1;
        let is_bottom: bool = row == 0;
        let is_top: bool = row == total_rows - 1;
        
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

    pub fn get_row_index(x: f32, normalized_cell_size: f32, rows: u32) -> u32 {
        Utils::clamp(((x + 1.0) / normalized_cell_size).floor().abs(), 0.0, rows as f32 - 1.0) as u32
    }

    pub fn get_column_index(y: f32, normalized_cell_size: f32, columns: u32) -> u32 {
        Utils::clamp(((1.0 - y) / normalized_cell_size).floor().abs(), 0.0, columns as f32 - 1.0) as u32
    }

    pub fn get_cell_index(row: u32, column: u32, rows: u32) -> u32 {
        row * rows + column
    }
}
