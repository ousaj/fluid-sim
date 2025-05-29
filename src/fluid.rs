use crate::utils::Utils;
use crate::position::Position;
use crate::color::Color;
use crate::cell::{Cell, CellType};
use crate::particle::Particle;
use crate::obstacle::Obstacle;

use gl::COLOR;
use rand::Rng;

use crate::config::*;

use rayon::prelude::*;

const FLUID_CELL: u32 = 0;
const AIR_CELL: u32 = 1;
const SOLID_CELL: u32 = 2;

// const U_FIELD: u32 = 0;
// const V_FIELD: u32 = 1;

// const MAX_PARTICLES: usize = CONFIG.particle.total as usize;
// const PARTICLE_RADIUS: f32 = CONFIG.particle.radius;
// const DENSITY: f32 = CONFIG.environment.density;
// const ROWS: f32 = CONFIG.scene.resolution;
// const COLUMNS: f32 = CONFIG.scene.resolution;
// const WIDTH: f32 = CONFIG.scene.width;
// const HEIGHT: f32 = CONFIG.scene.width;
// const SPACING: f32 = HEIGHT / CONFIG.scene.resolution; // Grid cell size.

// const F_NUM_X: f32 = (WIDTH / SPACING) + 1.0;
// const F_NUM_Y: f32 = (HEIGHT / SPACING) + 1.0;
// const F_NUM_CELLS: usize = (F_NUM_X * F_NUM_Y) as usize;

// const P_INV_SPACING: f32 = 1.0 / (2.2 * PARTICLE_RADIUS);

// const P_NUM_X: f32 = (WIDTH * P_INV_SPACING) + 1.0;
// const P_NUM_Y: f32 = (HEIGHT * P_INV_SPACING) + 1.0;
// const P_NUM_CELLS: usize = (P_NUM_X * P_NUM_Y) as usize;

// const GRID_LINES: usize = (CONFIG.scene.resolution * CONFIG.scene.resolution) as usize;

const RESOLUTION: usize = CONFIG.scene.resolution as usize;
const TOTAL_CELLS: usize = RESOLUTION * RESOLUTION as usize;
const TOTAL_PARTICLES: usize = CONFIG.particle.total as usize;
const GRAVITY: f32 = CONFIG.environment.gravity;
const TIME_STEP: f32 = CONFIG.scene.time_step;
const PARTICLE_RADIUS: f32 = CONFIG.particle.radius;
const WINDOW_SIZE: f32 = CONFIG.scene.window_size;
const CELL_SIZE: f32 = WINDOW_SIZE / CONFIG.scene.resolution;
const INVERSE_CELL_SIZE: f32 = 1.0 / CELL_SIZE;
const HALF_CELL_SIZE: f32 = CELL_SIZE * 0.5;
const OVER_RELAXATION: f32 = CONFIG.environment.over_relaxation;
const DENSITY: f32 = CONFIG.environment.density;
const PRESSURE_ITERATIONS: usize = CONFIG.scene.pressure_iterations as usize;
const PARTICLE_ITERATIONS: usize = CONFIG.scene.particle_iterations as usize;
const PRESSURE_SCALING_FACTOR: f32 = DENSITY * CELL_SIZE / TIME_STEP;
const COLOR_DIFFUSION_COEFF: f32 = 0.001;
const FLIP_RATIO: f32 = CONFIG.environment.flip_ratio;

const PARTICLE_DIAMETER: f32 = PARTICLE_RADIUS * 2.0;
const MIN_X: f32 = CELL_SIZE + PARTICLE_RADIUS;
const MIN_Y: f32 = CELL_SIZE + PARTICLE_RADIUS;
const MAX_X: f32 = WINDOW_SIZE - (CELL_SIZE + PARTICLE_RADIUS);
const MAX_Y: f32 = WINDOW_SIZE - (CELL_SIZE + PARTICLE_RADIUS);
const DRIFT_COMPENSATION_FACTOR: f32 = 1.0;

const MIN_DIST: f32 = 2.0 * PARTICLE_RADIUS;
const MIN_DIST_2: f32 = MIN_DIST * MIN_DIST;

fn clamp(x: f32, min: f32, max: f32) -> f32 {
    if x < min {
        min
    } else if x > max {
        max 
    } else {
        x 
    }
}

fn clamp_usize(value: usize, min: usize, max: usize) -> usize {
    if value < min {
        min
    } else if value > max {
        max
    } else {
        value
    }
}

pub struct Fluid {
    pub cell_vertices: [f32; 5 * TOTAL_CELLS],
    pub cell_indices: [u32; TOTAL_CELLS],
    pub particle_vertices: [f32; 5 * TOTAL_PARTICLES],
    pub particle_indices: [u32; TOTAL_PARTICLES],
    pub cell_size: f32,
    pub particle_diameter: f32,
    pub resolution: usize,
    pub particle_velocity: [f32; 2 * TOTAL_PARTICLES],
    pub cell_density: [f32; TOTAL_CELLS],
    pub rest_density: f32,
    pub cell_type: [u32; TOTAL_CELLS],
    pub num_cell_particles: [usize; TOTAL_CELLS],
    pub first_cell_particle: [usize; TOTAL_CELLS + 1],
    pub cell_particle_ids: [usize; TOTAL_PARTICLES],
    pub cell_pressure: [f32; TOTAL_CELLS],
    pub cell_u: [f32; TOTAL_CELLS],
    pub cell_v: [f32; TOTAL_CELLS],
    pub prev_cell_u: [f32; TOTAL_CELLS],
    pub prev_cell_v: [f32; TOTAL_CELLS],
    pub cell_solid_coeff: [f32; TOTAL_CELLS],
    pub du: [f32; TOTAL_CELLS],
    pub dv: [f32; TOTAL_CELLS],
    // pub cell_vertices:    pub min_x: f32,
    // pub particle_indices: [u32; MAX_PARTICLES],
    // pub f_num_y: f32,
    // pub h: f32,
    // pub f_inv_spacing: f32,
    // pub f_num_cells: f32,
    // pub u: [f32; F_NUM_CELLS],
    // pub v: [f32; F_NUM_CELLS],
    // pub du: [f32; F_NUM_CELLS],
    // pub dv: [f32; F_NUM_CELLS],
    // pub prev_u: [f32; F_NUM_CELLS],
    // pub prev_v: [f32; F_NUM_CELLS],
    // pub s: [f32; F_NUM_CELLS],
    // pub cell_type: [u32; F_NUM_CELLS],
    // pub cell_color: [f32; 3 * F_NUM_CELLS],
    // // Particles.

    // pub max_particles: f32,
    // pub particle_pos: [f32; 2 * MAX_PARTICLES],
    // pub particle_color: [f32; 3 * MAX_PARTICLES],
    // pub particle_vel: [f32; 2 * MAX_PARTICLES],
    // pub particle_density: [f32; F_NUM_CELLS],
    // pub particle_rest_density: f32,
    // pub particle_radius: f32,
    // pub p_inv_spacing: f32,
    // pub p_num_x: f32,
    // pub p_num_y: f32,
    // pub p_num_cells: f32,
    // pub num_cell_particles: [usize; P_NUM_CELLS],
    // pub first_cell_particle: [usize; 1 + P_NUM_CELLS],
    // pub cell_particle_ids: [usize; MAX_PARTICLES],
    // pub num_particles: u32,
    // pub particle_diameter: f32,
}


impl Fluid {
    pub fn new() -> Self {
        let resolution: usize = CONFIG.scene.resolution as usize;
        let window_size: f32 = CONFIG.scene.window_size;
        
        let cell_size: f32 = CONFIG.scene.window_size / CONFIG.scene.resolution;
        let mut cell_vertices: [f32; 5 * TOTAL_CELLS] = [0.0; 5 * TOTAL_CELLS];
        let mut cell_indices: [u32; TOTAL_CELLS] = [0; TOTAL_CELLS];
        let mut cell_type: [u32; TOTAL_CELLS] = [0; TOTAL_CELLS];
        let mut cell_solid_coeff: [f32; TOTAL_CELLS] = [0.0; TOTAL_CELLS];
        for i in 0..TOTAL_CELLS {
            let row: usize = i / resolution;
            let column: usize = i % resolution;

            let x: f32 = Utils::get_pixel(row, cell_size);
            let y: f32 = Utils::get_pixel(column,cell_size) + cell_size;

            cell_indices[i] = i as u32;

            cell_vertices[i * 5] = x;
            cell_vertices[i * 5 + 1] = y;

            if Utils::is_border(row, column, resolution) {
                cell_vertices[i * 5 + 2] = 0.3;
                cell_vertices[i * 5 + 3] = 0.3;
                cell_vertices[i * 5 + 4] = 0.3;

                cell_solid_coeff[i] = 0.0;
                cell_type[i] = SOLID_CELL;
            } else {
                cell_vertices[i * 5 + 2] = 0.8;
                cell_vertices[i * 5 + 3] = 0.8;
                cell_vertices[i * 5 + 4] = 0.8;

                cell_solid_coeff[i] = 1.0;
                cell_type[i] = FLUID_CELL; // AIR_CELL
            }
        }

        let particle_radius: f32 = CONFIG.particle.radius;
        let particle_diameter: f32 = CONFIG.particle.radius * 2.0;
        let mut particle_row: usize = 0;
        let mut particle_column: usize = 0;
        let particles_per_row: usize = (((CONFIG.scene.window_size - (cell_size * 2.0)) / particle_diameter).floor() - 1.0) as usize;
        let mut particle_vertices: [f32; 5 * TOTAL_PARTICLES] = [0.0; 5 * TOTAL_PARTICLES];
        let mut particle_indices: [u32; TOTAL_PARTICLES] = [0; TOTAL_PARTICLES];
        for i in 0..TOTAL_PARTICLES {
            if particle_column >= particles_per_row {
                particle_column = 0;
                particle_row += 1;
            }

            let stagger_x: f32 = if particle_row % 2 == 0 { 0.0 } else { particle_diameter / 2.0 };

            let x: f32 = cell_size + particle_radius + (particle_column as f32 * particle_diameter) + stagger_x;
            let y: f32 = CONFIG.scene.window_size - (cell_size + particle_radius + (particle_row as f32 * particle_diameter));

            // Position.
            particle_vertices[i * 5] = x;
            particle_vertices[i * 5 + 1] = y;

            // Velocity.
            particle_vertices[i * 5 + 2] = 0.0;
            particle_vertices[i * 5 + 3] = 0.0;
            particle_vertices[i * 5 + 4] = 3.0;

            particle_indices[i] = i as u32;

            particle_column += 1;
        }

        // let row: usize = 23;
        // let column: usize = 23;
        // let cell_number: usize = row * resolution + column;
        // cell_vertices[cell_number * 5 + 2] = 0.0;
        // cell_vertices[cell_number * 5 + 3] = 0.0;
        // cell_vertices[cell_number * 5 + 4] = 0.0;

        let particle_velocity: [f32; 2 * TOTAL_PARTICLES] = [0.0; 2 * TOTAL_PARTICLES];
        let cell_density: [f32; TOTAL_CELLS] = [0.0; TOTAL_CELLS];
        let rest_density: f32 = 0.0;

        let num_cell_particles: [usize; TOTAL_CELLS] = [0; TOTAL_CELLS];
        let first_cell_particle: [usize; TOTAL_CELLS + 1] = [0; TOTAL_CELLS + 1];
        let cell_particle_ids: [usize; TOTAL_PARTICLES] = [0; TOTAL_PARTICLES];
        let cell_pressure: [f32; TOTAL_CELLS] = [0.0; TOTAL_CELLS];
        let cell_u: [f32; TOTAL_CELLS] = [0.0; TOTAL_CELLS];
        let cell_v: [f32; TOTAL_CELLS] = [0.0; TOTAL_CELLS];
        let prev_cell_u: [f32; TOTAL_CELLS] = [0.0; TOTAL_CELLS];
        let prev_cell_v: [f32; TOTAL_CELLS] = [0.0; TOTAL_CELLS];
        let du: [f32; TOTAL_CELLS] = [0.0; TOTAL_CELLS];
        let dv: [f32; TOTAL_CELLS] = [0.0; TOTAL_CELLS];

        Fluid {
            cell_vertices,
            cell_indices,
            particle_vertices,
            particle_indices,
            cell_size,
            particle_diameter,
            resolution,
            particle_velocity,
            cell_density,
            rest_density,
            cell_type,
            num_cell_particles,
            first_cell_particle,
            cell_particle_ids,
            cell_pressure,
            cell_u,
            cell_v,
            prev_cell_u,
            prev_cell_v,
            cell_solid_coeff,
            du,
            dv,
        }
        // let mut cell_vertices: [f32; 5 * F_NUM_CELLS] = [0.0; 5 * F_NUM_CELLS];
        // let mut cell_indices: [u32; F_NUM_CELLS] = [0; F_NUM_CELLS];

        // let mut grid_vertices: [f32; 5 * GRID_LINES] = [0.0; 5 * GRID_LINES];
        // let mut grid_indices: [u32; 2 * GRID_LINES] = [0; 2 * GRID_LINES];

        // let mut particle_vertices: [f32; 5 * MAX_PARTICLES] = [0.0; 5 * MAX_PARTICLES];
        // let mut particle_indices: [u32; MAX_PARTICLES] = [0; MAX_PARTICLES];

        // let density: f32 = DENSITY;
        // let f_num_x: f32 = (WIDTH / SPACING).floor() + 1.0;
        // let f_num_y: f32 = (HEIGHT / SPACING).floor() + 1.0;
        // let h: f32 = (WIDTH / f_nu                // cell_vertices[i * 5 + 2] = 0.3;
                // cell_vertices[i * 5 + 3] = 0.3;
                // cell_vertices[i * 5 + 4] = 0.3;
        // let u: [f32; F_NUM_CELLS] = [0.0; F_NUM_CELLS];
        // let v: [f32; F_NUM_CELLS] = [0.0; F_NUM_CELLS];
        // let du: [f32; F_NUM_CELLS] = [0.0; F_NUM_CELLS];
        // let dv: [f32; F_NUM_CELLS] = [0.0; F_NUM_CELLS];
        // let prev_u: [f32; F_NUM_CELLS] = [0.0; F_NUM_CELLS];
        // let prev_v: [f32; F_NUM_CELLS] = [0.0; F_NUM_CELLS];
        // let p: [f32; F_NUM_CELLS] = [0.0; F_NUM_CELLS];
        // let mut s: [f32; F_NUM_CELLS] = [0.0; F_NUM_CELLS];
        // let mut cell_type: [u32; F_NUM_CELLS] = [0; F_NUM_CELLS];
        // let cell_color: [f32; 3 * F_NUM_CELLS] = [0.0; 3 * F_NUM_CELLS];

        // let max_particles: f32 = MAX_PARTICLES as f32;

        // let mut particle_pos: [f32; 2 * MAX_PARTICLES] = [0.0; 2 * MAX_PARTICLES];
        // let mut particle_color: [f32; 3 * MAX_PARTICLES] = [0.0; 3 * MAX_PARTICLES];
        // for i in 0..MAX_PARTICLES {
        //     particle_color[3 * i + 2]  = 1.0;total_rowsS];
        // let particle_density: [f32; F_NUM_CELLS] = [0.0; F_NUM_CELLS];
        // let particle_rest_density: f32 =  0.0;
        // let particle_radius: f32 = PARTICLE_RADIUS;
        // let particle_diameter: f32 = particle_radius * 2.0;
        // let p_inv_spacing: f32 = 1.0 / (2.2 * PARTICLE_RADIUS);
        // let p_num_x: f32 = (WIDTH * p_inv_spacing).floor() + 1.0;
        // let p_num_y: f32 Self_ids: [usize; MAX_PARTICLES] = [0; MAX_PARTICLES];

        // let num_particles: u32 = 0;

        // // for i in 0..f_num_x as usize {
        // //     for j in 0..f_num_y as usize {
        // //         let s: f32 = 1.0;
        // //         if (i == 0 || i == f_num_x - 1.0 || j == 0 || )
        // //     }
        // // }

        // for i in 0..F_NUM_CELLS {
        //     let row: u32 = i as u32 / CONFIG.scene.resolution as u32;
        //     let column: u32 = i as u32 % CONFIG.scene.resolution as u32;

        //     let x: f32 = Utils::get_pixel(row as f32, SPACING);
        //     let y: f32 = Utils::get_pixel(column as f32, SPACING) + h;

        //     cell_indices[i] = i as u32;

        //     cell_vertices[i * 5] = x;
        //     cell_vertices[i * 5 + 1] = y;

        //     if Utils::is_border(row, column, ROWS as u32, COLUMNS as u32) {
        //         cell_vertices[i * 5 + 2] = 0.3;
        //         cell_vertices[i * 5 + 3] = 0.3;
        //         cell_vertices[i * 5 + 4] = 0.3;

        //         s[i] = 0.0;
        //         cell_type[i] = SOLID_CELL;
        //     } else {
        //         cell_vertices[i * 5 + 2] = 0.8;
        //         cell_vertices[i * 5 + 3] = 0.8;
        //         cell_vertices[i * 5 + 4] = 0.8;

        //         s[i] = 1.0;
        //         cell_type[i] = FLUID_CELL; // AIR_CELL
        //     }
        // }

        // let cols = ((WIDTH - (h * 2.0)) / particle_diameter).floor() as usize;

        // let mut row: usize = 0;
        // let mut col: usize = 0;
        // for i in 0..MAX_PARTICLES {
        //     if col >= cols {
        //         col = 0;
        //         row += 1;
        //     }

        //     let stagger_x: f32 = if row % 2 == 0 { 0.0 } else { particle_diameter / 2.0 };

        //     let x: f32 = h + CONFIG.particle.radius + (col as f32 * particle_diameter) + stagger_x;
        //     let y: f32 = CONFIG.scene.width - (h + CONFIG.particle.radius + (row as f32 * particle_diameter));

        //     particle_pos[i * 2] = x;
        //     particle_pos[i * 2 + 1] = y;
        //     particle_color[i * 3] = 0.0;
        //     particle_color[i * 3 + 1] = 0.0;
        //     particle_color[i * 3 + 2] = 3.0;

        //     particle_indices[i] = i as u32;

        //     col += 1;
        // }

        // Fluid {
        //     cell_vertices,
        //     cell_indices,
        //     grid_vertices,
        //     grid_indices,
        //     particle_vertices,
        //     particle_indices,
        //     density,
        //     f_num_x,
        //     f_num_y,
        //     h,
        //     f_inv_spacing,
        //     f_num_cells,
        //     u,
        //     v,
        //     du,
        //     dv,
        //     prev_u,
        //     prev_v,
        //     p,
        //     s,
        //     cell_type,
        //     cell_color,
        //     max_particles,
        //     particle_pos,
        //     particle_color,
        //     particle_vel,
        //     particle_density,
        //     particle_rest_density,
        //     particle_radius,
        //     p_inv_spacing,
        //     p_num_x,
        //     p_num_y,
        //     p_num_cells,
        //     num_cell_particles,
        //     first_cell_particle,
        //     cell_particle_ids,
        //     num_particles,
        //     particle_diameter,
        // }
    }

    pub fn integrate_particles(&mut self) {
        for i in 0..TOTAL_PARTICLES {
            self.particle_velocity[2 * i + 1] += TIME_STEP * GRAVITY;

            self.particle_vertices[5 * i] += self.particle_velocity[2 * i] * TIME_STEP;
            self.particle_vertices[5 * i + 1] += self.particle_velocity[2 * i + 1] * TIME_STEP;
        }
    }

    pub fn push_particles_apart(&mut self) {
        self.num_cell_particles.fill(0);

        for i in 0..TOTAL_PARTICLES {
            let x: f32 = self.particle_vertices[5 * i];
            let y: f32 = self.particle_vertices[5 * i + 1];

            let xi: usize = clamp_usize((x * INVERSE_CELL_SIZE).floor() as usize, 1, RESOLUTION - 2);
            let yi: usize = clamp_usize((y * INVERSE_CELL_SIZE).floor() as usize, 1, RESOLUTION - 2);

            let cell_nr: usize = xi * RESOLUTION + yi;
            self.num_cell_particles[cell_nr] += 1;
        }

        let mut first: usize = 0;
        for i in 0..TOTAL_CELLS { // TODO: Ajustar para que no entre a las celdas sólidas porque no hace falta.
            first += self.num_cell_particles[i];
            self.first_cell_particle[i] = first;
        }

        self.first_cell_particle[TOTAL_CELLS] = first;

        for i in 0..TOTAL_PARTICLES {
            let x: f32 = self.particle_vertices[5 * i];
            let y: f32 = self.particle_vertices[5 * i + 1];

            let xi: usize = clamp_usize((x * INVERSE_CELL_SIZE).floor() as usize, 1, RESOLUTION - 2);
            let yi: usize = clamp_usize((y * INVERSE_CELL_SIZE).floor() as usize, 1, RESOLUTION - 2);

            let cell_nr: usize = xi * RESOLUTION + yi;
            self.first_cell_particle[cell_nr] -= 1;
            self.cell_particle_ids[self.first_cell_particle[cell_nr]] = i;
        }

        for _iter in 0..PARTICLE_ITERATIONS {
            for i in 0..TOTAL_PARTICLES {
                let px: f32 = self.particle_vertices[5 * i];
                let py: f32 = self.particle_vertices[5 * i + 1];

                let pxi: usize = (px * INVERSE_CELL_SIZE).floor() as usize;
                let pyi: usize = (py * INVERSE_CELL_SIZE).floor() as usize;

                // TODO: Puede que aquí puedas aplicar el límite 1,1 23,23.
                let x0: usize = pxi.saturating_sub(1).max(1);
                let y0: usize = pyi.saturating_sub(1).max(1);
                let x1: usize = (pxi + 1).min(RESOLUTION - 1);
                let y1: usize = (pyi + 1).min(RESOLUTION - 1);

                for xi in x0..=x1 {
                    for yi in y0..=y1 {
                        let cell_nr: usize = xi * RESOLUTION + yi;

                        let first: usize = self.first_cell_particle[cell_nr];
                        let last: usize = self.first_cell_particle[cell_nr + 1];

                        for j in first..last {
                            let id: usize = self.cell_particle_ids[j];

                            if id == i {
                                continue;
                            }

                            let qx: f32 = self.particle_vertices[5 * id];
                            let qy: f32 = self.particle_vertices[5 * id + 1];

                            let mut dx: f32 = qx - px;
                            let mut dy: f32 = qy - py;

                            let d2: f32 = dx * dx + dy * dy;

                            if d2 > MIN_DIST_2 || d2 == 0.0 {
                                continue;
                            }

                            let d: f32 = d2.sqrt();

                            let s: f32 = 0.5 * (MIN_DIST - d) / d;

                            dx *= s;
                            dy *= s;

                            self.particle_vertices[5 * i] -= dx;
                            self.particle_vertices[5 * i + 1] -= dy;
                            self.particle_vertices[5 * id] += dx;
                            self.particle_vertices[5 * id + 1] += dy;

                            for k in 0..3 {
                                let color0: f32 = self.particle_vertices[5 * i + k];
                                let color1: f32 = self.particle_vertices[5 * id + k];
                                let color: f32 = (color0 + color1) * 0.5;

                                self.particle_vertices[5 * i + k] = color0 + (color - color0) * COLOR_DIFFUSION_COEFF;
                                self.particle_vertices[5 * id + k] = color1 + (color - color1) * COLOR_DIFFUSION_COEFF;
                            }
                        }
                    }
                }
            }
        }
    }

    pub fn handle_particle_collisions(&mut self, obstacle: &Obstacle) {
        for i in 0..TOTAL_PARTICLES {
            let mut x: f32 = self.particle_vertices[i * 5];
            let mut y: f32 = self.particle_vertices[i * 5 + 1];

            if x < MIN_X {
                x = MIN_X;
                self.particle_velocity[2 * i] = 0.0;
            }

            if x > MAX_X {
                x = MAX_X;
                self.particle_velocity[2 * i] = 0.0;
            }


            if y < MIN_Y {
                y = MIN_Y;
                self.particle_velocity[2 * i + 1] = 0.0;
            }

            if y > MAX_Y {
                y = MAX_Y;
                self.particle_velocity[2 * i + 1] = 0.0;
            }

            self.particle_vertices[i * 5] = x;
            self.particle_vertices[i * 5 + 1] = y;
        }
    }

    pub fn transfer_velocities(&mut self, to_grid: bool) {
        // let n: usize = (self.f_num_y - 1.0) as usize;
        // let h: f32 = self.h;
        // let h1: f32 = self.f_inv_spacing;
        // let h2: f32 = 0.5 * h;

        if to_grid {
            self.prev_cell_u.copy_from_slice(&self.cell_u);
            self.prev_cell_v.copy_from_slice(&self.cell_v);

            self.du.fill(0.0);
            self.dv.fill(0.0);

            self.cell_u.fill(0.0);
            self.cell_v.fill(0.0);

            // Reseteamos los tipos de celdas fluidas.
            for i in 0..TOTAL_CELLS {
                self.cell_type[i] = if self.cell_solid_coeff[i] == 0.0 { SOLID_CELL } else { AIR_CELL };
            }

            for i in 0..TOTAL_PARTICLES {
                let x: f32 = self.particle_vertices[5 * i];
                let y: f32 = self.particle_vertices[5 * i + 1];

                let xi: usize = clamp_usize((x * INVERSE_CELL_SIZE).floor() as usize, 1 ,RESOLUTION - 2);
                let yi: usize = clamp_usize((y * INVERSE_CELL_SIZE).floor() as usize, 1, RESOLUTION - 2);
                let cell_nr: usize = xi * RESOLUTION + yi;

                // Pasamos a fluidas las celdas con partículas.
                if self.cell_type[cell_nr] == AIR_CELL {
                    self.cell_type[cell_nr] = FLUID_CELL;
                }
            }
        }

        for component in 0..2 {
            let dx: f32 = if component == 0 { 0.0 } else { HALF_CELL_SIZE };
            let dy: f32 = if component == 0 { HALF_CELL_SIZE } else { 0.0 };

            let (f, prev_f, d) = match component {
                0 => (&mut self.cell_u, &mut self.prev_cell_u, &mut self.du),
                1 => (&mut self.cell_v, &mut self.prev_cell_v, &mut self.dv),
                _ => unreachable!(),
            };

            for i in 0..TOTAL_PARTICLES {
                let mut x: f32 = self.particle_vertices[5 * i];
                let mut y: f32 = self.particle_vertices[5 * i + 1];

                // Mantenemos la partícula dentro del dominio.
                x = clamp(x, MIN_X, MAX_X);
                y = clamp(y, MIN_Y, MAX_Y);

                let x0: usize = (((x - dx) * INVERSE_CELL_SIZE).floor() as usize).min(RESOLUTION - 2);
                let y0: usize = (((y - dy) * INVERSE_CELL_SIZE).floor() as usize).min(RESOLUTION - 2);
                let x1: usize = (x0 + 1).min(RESOLUTION - 2);
                let y1: usize = (y0 + 1).min(RESOLUTION - 2);
                
                let tx: f32 = ((x - dx) - x0 as f32 * CELL_SIZE) * INVERSE_CELL_SIZE;
                let ty: f32 = ((y - dy) - y0 as f32 * CELL_SIZE) * INVERSE_CELL_SIZE;
                let sx: f32 = 1.0 - tx;
                let sy: f32 = 1.0 - ty;

                let d0: f32 = sx * sy; // Peso para la celda inferior izquierda.
                let d1: f32 = tx * sy; // Peso para la celda inferior derecha.
                let d2: f32 = tx * ty; // Peso para la celda superior derecha.
                let d3: f32 = sx * ty; // Peso para la celda superior izquierda.

                let nr0: usize = x0 * RESOLUTION + y0; // Celda inferior izquierda.
                let nr1: usize = x1 * RESOLUTION + y0; // Celda inferior derecha.
                let nr2: usize = x1 * RESOLUTION + y1; // Celda superior derecha.
                let nr3: usize = x0 * RESOLUTION + y1; // Celda superior izquierda.

                if to_grid {
                    let particle_velocity: f32 = self.particle_velocity[2 * i + component];

                    f[nr0] += particle_velocity * d0;
                    f[nr1] += particle_velocity * d1;
                    f[nr2] += particle_velocity * d2;
                    f[nr3] += particle_velocity * d3;

                    d[nr0] += d0;
                    d[nr1] += d1;
                    d[nr2] += d2;
                    d[nr3] += d3;
                } else { // To particle.
                    let offset: usize = if component == 0 { RESOLUTION } else { 1 };
                    let valid0: f32 = if self.cell_type[nr0] != AIR_CELL || self.cell_type[nr0 - offset] != AIR_CELL { 1.0 } else { 0.0 };
                    let valid1: f32 = if self.cell_type[nr1] != AIR_CELL || self.cell_type[nr1 - offset] != AIR_CELL { 1.0 } else { 0.0 };
                    let valid2: f32 = if self.cell_type[nr2] != AIR_CELL || self.cell_type[nr2 - offset] != AIR_CELL { 1.0 } else { 0.0 };
                    let valid3: f32 = if self.cell_type[nr3] != AIR_CELL || self.cell_type[nr3 - offset] != AIR_CELL { 1.0 } else { 0.0 };

                    let v: f32 = self.particle_velocity[2 * i + component];
                    let d: f32 = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

                    if d > 0.0 {
                        let pic_v = (valid0 * d0 * f[nr0] + valid1 * d1 * f[nr1] + valid2 * d2 * f[nr2] + valid3 * d3 * f[nr3]) / d;
                        let corr = (valid0 * d0 * (f[nr0] - prev_f[nr0]) + valid1 * d1 * (f[nr1] - prev_f[nr1])
                            + valid2 * d2 * (f[nr2] - prev_f[nr2]) + valid3 * d3 * (f[nr3] - prev_f[nr3])) / d;
                        let flip_v = v + corr;

                        self.particle_velocity[2 * i + component] = (1.0 - FLIP_RATIO) * pic_v + FLIP_RATIO * flip_v;
                    }
                }
            }

            if to_grid {
                for i in 0..TOTAL_CELLS {
                // for i in 0..f.len() {
                    if d[i] > 0.0 {
                        f[i] /= d[i];
                    }
                }

                // Restaurar velocidades en celdas sólidas.
                // TODO: Por qué hace esto dos veces?
                for i in 0..RESOLUTION {
                    for j in 0..RESOLUTION {
                        let solid: bool = self.cell_type[i * RESOLUTION + j] == SOLID_CELL;

                        // TODO: Siempre va a ser mayor a cero.
                        if solid || (i > 0 && self.cell_type[(i - 1) * RESOLUTION + j] == SOLID_CELL) {
                            self.cell_u[i * RESOLUTION + j] = self.prev_cell_u[i * RESOLUTION + j];
                        }

                        if solid || (j > 0 && self.cell_type[i * RESOLUTION + (j - 1)] == SOLID_CELL) {
                            self.cell_v[i * RESOLUTION + j] = self.prev_cell_v[i * RESOLUTION + j];
                        }
                    }
                }
            }
        }

        // for component in 0..2 {
        //     let dx: f32 = if component == 0 { 0.0 } else { h2 };
        //     let dy: f32 = if component == 0 { h2 } else { 0.0 };

        //     let (f, prev_f, d) = match component {
        //         0 => (&mut self.u, &self.prev_u, &mut self.du),
        //         1 => (&mut self.v, &self.prev_v, &mut self.dv),
        //         _ => unreachable!(),
        //     };

        //     for i in 0..MAX_PARTICLES {
        //         let mut x: f32 = self.particle_pos[2 * i];
        //         let mut y: f32 = self.particle_pos[2 * i + 1];

        //         x = clamp(x, h, (self.f_num_x - 1.0) * h);
        //         y = clamp(y, h, (self.f_num_y - 1.0) * h);

        //         let x0: usize = ((x - dx) * h1).floor().min(self.f_num_x) as usize;
        //         let tx: f32 = ((x - dx) - x0 as f32 * h) * h1;
        //         let x1: usize = (x0 + 1).min((self.f_num_x - 2.0) as usize);

        //         let y0: usize = ((y - dy) * h1).floor().min(self.f_num_y) as usize;
        //         let ty: f32 = ((y - dy) - y0 as f32 * h) * h1;
        //         let y1: usize = (y0 + 1).min((self.f_num_y) as usize);

        //         let sx: f32 = 1.0 - tx;
        //         let sy: f32 = 1.0 - ty;

        //         let d0: f32 = sx * sy;
        //         let d1: f32 = tx * sy;
        //         let d2: f32 = tx * ty;
        //         let d3: f32 = sx * ty;

        //         let nr0: usize = x0 * n + y0;
        //         let nr1: usize = x1 * n + y0;
        //         let nr2: usize = x1 * n + y1;
        //         let nr3: usize = x0 * n + y1;

        //         if to_grid {
        //             let pv = self.particle_vel[2 * i + component];

        //             f[nr0] += pv * d0; d[nr0] += d0;
        //             f[nr1] += pv * d1; d[nr1] += d1;
        //             f[nr2] += pv * d2; d[nr2] += d2;
        //             f[nr3] += pv * d3; d[nr3] += d3;
        //         } else {
        //             let offset: usize = if component == 0 { n as usize } else { 1 };

        //             let valid0: f32 = if self.cell_type[nr0] != AIR_CELL || self.cell_type[nr0 - offset] != AIR_CELL { 1.0 } else { 0.0 };
        //             let valid1: f32 = if self.cell_type[nr1] != AIR_CELL || self.cell_type[nr1 - offset] != AIR_CELL { 1.0 } else { 0.0 };
        //             let valid2: f32 = if self.cell_type[nr2] != AIR_CELL || self.cell_type[nr2 - offset] != AIR_CELL { 1.0 } else { 0.0 };
        //             let valid3: f32 = if self.cell_type[nr3] != AIR_CELL || self.cell_type[nr3 - offset] != AIR_CELL { 1.0 } else { 0.0 };

        //             let v: f32 = self.particle_vel[2 * i + component];
        //             let d: f32 = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

        //             if d > 0.0 {
        //                 let pic_v = (valid0 * d0 *f[nr0] + valid1 * d1 * f[nr1] + valid2 * d2 * f[nr2] + valid3 * d3 * f[nr3]) / d;
        //                 let corr = (valid0 * d0 * (f[nr0] - prev_f[nr0]) + valid1 * d1 * (f[nr1] - prev_f[nr1]) + valid2 * d2 * (f[nr2] - prev_f[nr2]) + valid3 * d3 * (f[nr3] - prev_f[nr3])) / d;
        //                 let flip_v = v + corr;

        //                 self.particle_vel[2 * i + component] = (1.0 - flip_ratio) * pic_v + flip_ratio * flip_v;
        //             }
        //         }
        //     }

        //     if to_grid {
        //         for i in 0..f.len() {
        //             if d[i] > 0.0 {
        //                 f[i] /= d[i];
        //             }
        //         }

        //         let f_num_x: usize = (self.f_num_x) as usize;
        //         let f_num_y: usize = (self.f_num_y) as usize;
        //         for i in 0..f_num_x {
        //             for j in 0..f_num_y {
        //                 let solid = self.cell_type[i * n + j] == SOLID_CELL;

        //                 if solid || (i > 0 && self.cell_type[(i - 1) * n + j] == SOLID_CELL) {
        //                     self.u[i * n + j] = self.prev_u[i * n + j];
        //                 }

        //                 if solid || (j > 0 && self.cell_type[i * n + j - 1] == SOLID_CELL) {
        //                     self.v[i * n + j] = self.prev_v[i * n + j];
        //                 }
        //             }
        //         }
        //     }
        // }
    }

    pub fn update_density(&mut self) {
        self.cell_density.fill(0.0);

        for i in 0..TOTAL_PARTICLES {
            let mut x: f32 = self.particle_vertices[5 * i];
            let mut y: f32 = self.particle_vertices[5 * i + 1];

            x = clamp(x, MIN_X, MAX_X);
            y = clamp(y, MIN_Y, MAX_Y);

            // Bottom row of 2x2 grid.
            let x0: usize = ((x - HALF_CELL_SIZE) * INVERSE_CELL_SIZE).floor() as usize;
            // Top row of 2x2 grid.
            let x1: usize = (x0 + 1).min(RESOLUTION - 1); // TODO: Puede que sea 2.
            // Left column of 2x2 grid.
            let y0: usize = ((y - HALF_CELL_SIZE) * INVERSE_CELL_SIZE).floor() as usize;
            // Right column of 2x2 grid.
            let y1: usize = (y0 + 1).min(RESOLUTION - 1);
            // Percentage of the way the shifted particle is from the left cell to the right cell.
            let tx: f32 = ((x - HALF_CELL_SIZE) - x0 as f32 * CELL_SIZE) * INVERSE_CELL_SIZE;
            // Percentage of the way the shifted particle is from the bottom cell to the top cell.
            let ty: f32 = ((y - HALF_CELL_SIZE) - y0 as f32 * CELL_SIZE) * INVERSE_CELL_SIZE;

            let sx: f32 = 1.0 - tx;
            let sy: f32 = 1.0 - ty;

            // TODO: Creo que estos ifs no hacen falta.

            // Bottom-left.
            if x0 < RESOLUTION && y0 < RESOLUTION {
                self.cell_density[x0 * RESOLUTION + y0] += sx * sy;
            }

            // Bottom-right.
            if x1 < RESOLUTION && y0 < RESOLUTION {
                self.cell_density[x1 * RESOLUTION + y0] += tx * sy;
            }

            // Top-right.
            if x1 < RESOLUTION && y1 < RESOLUTION {
                self.cell_density[x1 * RESOLUTION + y1] += tx * ty;
            }

            // Top-left.
            if x0 < RESOLUTION && y1 < RESOLUTION {
                self.cell_density[x0 * RESOLUTION + y1] += sx * ty;
            }
        }

        // Compute average density (rest density) if not defined.
        // This becomes the reference density to compare against in pressure correction,
        // helping maintain the fluid's natural spacing.
        if self.rest_density == 0.0 {
            let mut density_sum: f32 = 0.0;
            let mut fluid_cells: f32 = 0.0;
            for i in 0..TOTAL_CELLS {
                if self.cell_type[i] == FLUID_CELL {
                    density_sum += self.cell_density[i];
                    fluid_cells += 1.0;
                }
            }

            if fluid_cells > 0.0 {
                self.rest_density = density_sum / fluid_cells;
            }
        }
    }

    // Corrección de presión.
    pub fn solve_incompressibility(&mut self) {
        self.cell_pressure.fill(0.0);
        self.prev_cell_u.copy_from_slice(&self.cell_u);
        self.prev_cell_v.copy_from_slice(&self.cell_v);

        for _iter in 0..PRESSURE_ITERATIONS {
            // Esto lo hacemos para evitar computar las celdas de borde.
            for i in 1..(RESOLUTION-1) {
                for j in 1..(RESOLUTION-1) {
                    if self.cell_type[i * RESOLUTION + j] != FLUID_CELL {
                        continue;
                    }

                    let cell: usize = i * RESOLUTION + j;
                    let left_cell: usize = (i - 1) * RESOLUTION + j;
                    let right_cell: usize = (i + 1) * RESOLUTION + j;
                    let bottom_cell: usize = i * RESOLUTION + (j - 1);
                    let top_cell: usize = i * RESOLUTION + j + 1;

                    let left_cell_solid_coeff: f32 = self.cell_solid_coeff[left_cell];
                    let right_cell_solid_coeff: f32 = self.cell_solid_coeff[right_cell];
                    let bottom_cell_solid_coeff: f32 = self.cell_solid_coeff[bottom_cell];
                    let top_cell_solid_coeff: f32 = self.cell_solid_coeff[top_cell];

                    let total_solid_coeff: f32 = left_cell_solid_coeff + right_cell_solid_coeff + bottom_cell_solid_coeff + top_cell_solid_coeff;

                    if total_solid_coeff == 0.0 {
                        continue;
                    }

                    // Cálculo de divergencia.
                    // Cuánto fluido entra o sale de una celda.
                    // Se mide con las diferencias de velocidad entre caras.
                    let mut divergence: f32 = 
                          self.cell_u[right_cell]
                        - self.cell_u[cell]
                        + self.cell_v[top_cell]
                        - self.cell_v[cell];
                    
                    // Corrección de divergencia.
                    // Esto ayuda a que las partículas no se amontonen.
                    // Drift control.
                    if self.rest_density > 0.0 {
                        let compression = self.cell_density[i * RESOLUTION + j] - self.rest_density;
                        if compression > 0.0 {
                            divergence = divergence - DRIFT_COMPENSATION_FACTOR * compression;
                        }
                    }

                    // Cálculo de la corrección de presión.
                    // Se calcula cuánto corregir la presión.
                    // OVER_RELAXATION acelera la convergencia.
                    let mut pressure: f32 = -divergence / total_solid_coeff;
                    pressure *= OVER_RELAXATION;
                    self.cell_pressure[cell] += pressure * PRESSURE_SCALING_FACTOR;

                    self.cell_u[cell] -= left_cell_solid_coeff * pressure;
                    self.cell_u[right_cell] += right_cell_solid_coeff * pressure;
                    self.cell_v[cell] -= bottom_cell_solid_coeff * pressure;
                    self.cell_v[top_cell] += top_cell_solid_coeff * pressure;
                }
            }    
        }
    }

    pub fn update_particle_colors(&mut self) {
        for i in 0..TOTAL_PARTICLES {
            self.particle_vertices[5 * i + 2] = clamp(self.particle_vertices[5 * i + 2] - 0.01, 0.0, 1.0);
            self.particle_vertices[5 * i + 3] = clamp(self.particle_vertices[5 * i + 3] - 0.01, 0.0, 1.0);
            self.particle_vertices[5 * i + 4] = clamp(self.particle_vertices[5 * i + 4] - 0.01, 0.0, 1.0);

            let x: f32 = self.particle_vertices[5 * i];
            let y: f32 = self.particle_vertices[5 * i + 1];

            // RESOLUTION - 1 te da el índice de celda máximo del dominio (borde). RESOLUTION - 2 te da el índice de celda máximo válido.
            // Con esto sacas la celda en la que se encuentra una partícula.
            let particle_row: usize = clamp_usize((x * INVERSE_CELL_SIZE).floor() as usize, 1, RESOLUTION - 2); 
            let particle_column: usize = clamp_usize((y * INVERSE_CELL_SIZE).floor() as usize, 1, RESOLUTION - 2);
            let cell_number: usize = particle_row * RESOLUTION + particle_column;

            if self.rest_density > 0.0 {
                let relative_density: f32 = self.cell_density[cell_number] / self.rest_density;
                if relative_density > 0.7 {
                    // self.particle_vertices[5 * i + 2] = 0.8;
                    // self.particle_vertices[5 * i + 3] = 0.8;
                    // self.particle_vertices[5 * i + 4] = 1.0;
                }
            }
        }
    }

    pub fn update_cell_colors(&mut self) {
        for i in 0..TOTAL_CELLS {
            let cell_type: u32 = self.cell_type[i];
            match cell_type {
                SOLID_CELL => {
                    self.cell_vertices[5 * i + 2] = 0.5;
                    self.cell_vertices[5 * i + 3] = 0.5;
                    self.cell_vertices[5 * i + 4] = 0.5;
                },
                AIR_CELL => {
                    self.cell_vertices[5 * i + 2] = 0.0;
                    self.cell_vertices[5 * i + 3] = 0.0;
                    self.cell_vertices[5 * i + 4] = 0.0;
                },
                FLUID_CELL => {
                    self.cell_vertices[5 * i + 2] = 0.8;
                    self.cell_vertices[5 * i + 3] = 0.8;
                    self.cell_vertices[5 * i + 4] = 1.0;
                },
                _ => unreachable!()
            }
        }

        // for i in 0..F_NUM_CELLS {
        //     if self.cell_type[i] == SOLID_CELL {
        //         self.cell_color[3 * i] = 0.5;
        //         self.cell_color[3 * i + 1] = 0.5;
        //         self.cell_color[3 * i + 2] = 0.5;
        //     } else if self.cell_type[i] == FLUID_CELL {
        //         if self.particle_rest_density > 0.0 {
        //             self.particle_density[i] /= self.particle_rest_density;
        //         }

        //         let min_val: f32 = 0.0;
        //         let max_val: f32 = 2.0;

        //         self.particle_density[i] = (self.particle_density[i].max(min_val)).min(max_val - 0.0001);

        //         let d = max_val - min_val;
        //         self.particle_density[i] = if d == 0.0 { 0.5 } else { (self.particle_density[i] - min_val) / d };

        //         let m = 0.25;
        //         let num = (self.particle_density[i] / m).floor() as u32;
        //         let s = (self.particle_density[i] - (num as f32) * m) / m;

        //         let (r, g, b) = match num {
        //             0 => (0.0, s, 1.0),
        //             1 => (0.0, 1.0, 1.0 - s),
        //             2 => (s, 1.0, 0.0),
        //             3 => (1.0, 1.0 - s, 0.0),
        //             _ => (1.0, 0.0, 0.0), // Fallback for unexpected input
        //         };

        //         self.cell_color[3 * i] = r;
        //         self.cell_color[3 * i + 1] = g;
        //         self.cell_color[3 * i + 2] = b;
        //     }
        // }
    }
}
