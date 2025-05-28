use crate::utils::Utils;
use crate::position::Position;
use crate::color::Color;
use crate::cell::{Cell, CellType};
use crate::particle::Particle;
use crate::obstacle::Obstacle;

use rand::Rng;

use crate::config::*;

use rayon::prelude::*;

const FLUID_CELL: u32 = 0;
const AIR_CELL: u32 = 1;
const SOLID_CELL: u32 = 2;

const U_FIELD: u32 = 0;
const V_FIELD: u32 = 1;

const MAX_PARTICLES: usize = CONFIG.particle.total as usize;
const PARTICLE_RADIUS: f32 = CONFIG.particle.radius;
const DENSITY: f32 = CONFIG.environment.density;
const ROWS: f32 = CONFIG.scene.resolution;
const COLUMNS: f32 = CONFIG.scene.resolution;
const WIDTH: f32 = CONFIG.scene.width;
const HEIGHT: f32 = CONFIG.scene.width;
const SPACING: f32 = HEIGHT / CONFIG.scene.resolution; // Grid cell size.

const F_NUM_X: f32 = (WIDTH / SPACING) + 1.0;
const F_NUM_Y: f32 = (HEIGHT / SPACING) + 1.0;
const F_NUM_CELLS: usize = (F_NUM_X * F_NUM_Y) as usize;

const P_INV_SPACING: f32 = 1.0 / (2.2 * PARTICLE_RADIUS);

const P_NUM_X: f32 = (WIDTH * P_INV_SPACING) + 1.0;
const P_NUM_Y: f32 = (HEIGHT * P_INV_SPACING) + 1.0;
const P_NUM_CELLS: usize = (P_NUM_X * P_NUM_Y) as usize;

const GRID_LINES: usize = (CONFIG.scene.resolution * CONFIG.scene.resolution) as usize;

fn clamp(x: f32, min: f32, max: f32) -> f32 {
    if x < min {
        min
    } else if x > max {
        max 
    } else {
        x 
    }
}

pub struct Fluid {
    pub cell_vertices: [f32; 5 * F_NUM_CELLS],
    pub cell_indices: [u32; F_NUM_CELLS],
    pub grid_vertices: [f32; 5 * GRID_LINES],
    pub grid_indices: [u32; 2 * GRID_LINES],
    pub particle_vertices: [f32; 5 * MAX_PARTICLES],
    pub particle_indices: [u32; MAX_PARTICLES],

    pub density: f32,
    pub f_num_x: f32,
    pub f_num_y: f32,
    pub h: f32,
    pub f_inv_spacing: f32,
    pub f_num_cells: f32,
    pub u: [f32; F_NUM_CELLS],
    pub v: [f32; F_NUM_CELLS],
    pub du: [f32; F_NUM_CELLS],
    pub dv: [f32; F_NUM_CELLS],
    pub prev_u: [f32; F_NUM_CELLS],
    pub prev_v: [f32; F_NUM_CELLS],
    pub p: [f32; F_NUM_CELLS],
    pub s: [f32; F_NUM_CELLS],
    pub cell_type: [u32; F_NUM_CELLS],
    pub cell_color: [f32; 3 * F_NUM_CELLS],
    // Particles.

    pub max_particles: f32,
    pub particle_pos: [f32; 2 * MAX_PARTICLES],
    pub particle_color: [f32; 3 * MAX_PARTICLES],
    pub particle_vel: [f32; 2 * MAX_PARTICLES],
    pub particle_density: [f32; F_NUM_CELLS],
    pub particle_rest_density: f32,
    pub particle_radius: f32,
    pub p_inv_spacing: f32,
    pub p_num_x: f32,
    pub p_num_y: f32,
    pub p_num_cells: f32,
    pub num_cell_particles: [usize; P_NUM_CELLS],
    pub first_cell_particle: [usize; 1 + P_NUM_CELLS],
    pub cell_particle_ids: [usize; MAX_PARTICLES],
    pub num_particles: u32,
    pub particle_diameter: f32,
}


impl Fluid {
    pub fn new() -> Self {
        let mut cell_vertices: [f32; 5 * F_NUM_CELLS] = [0.0; 5 * F_NUM_CELLS];
        let mut cell_indices: [u32; F_NUM_CELLS] = [0; F_NUM_CELLS];

        let mut grid_vertices: [f32; 5 * GRID_LINES] = [0.0; 5 * GRID_LINES];
        let mut grid_indices: [u32; 2 * GRID_LINES] = [0; 2 * GRID_LINES];

        let mut particle_vertices: [f32; 5 * MAX_PARTICLES] = [0.0; 5 * MAX_PARTICLES];
        let mut particle_indices: [u32; MAX_PARTICLES] = [0; MAX_PARTICLES];

        let density: f32 = DENSITY;
        let f_num_x: f32 = (WIDTH / SPACING).floor() + 1.0;
        let f_num_y: f32 = (HEIGHT / SPACING).floor() + 1.0;
        let h: f32 = (WIDTH / f_num_x).max(HEIGHT / f_num_y);
        let f_inv_spacing: f32 = 1.0 / h;
        let f_num_cells: f32 = f_num_x * f_num_y;

        let u: [f32; F_NUM_CELLS] = [0.0; F_NUM_CELLS];
        let v: [f32; F_NUM_CELLS] = [0.0; F_NUM_CELLS];
        let du: [f32; F_NUM_CELLS] = [0.0; F_NUM_CELLS];
        let dv: [f32; F_NUM_CELLS] = [0.0; F_NUM_CELLS];
        let prev_u: [f32; F_NUM_CELLS] = [0.0; F_NUM_CELLS];
        let prev_v: [f32; F_NUM_CELLS] = [0.0; F_NUM_CELLS];
        let p: [f32; F_NUM_CELLS] = [0.0; F_NUM_CELLS];
        let mut s: [f32; F_NUM_CELLS] = [0.0; F_NUM_CELLS];
        let mut cell_type: [u32; F_NUM_CELLS] = [0; F_NUM_CELLS];
        let cell_color: [f32; 3 * F_NUM_CELLS] = [0.0; 3 * F_NUM_CELLS];

        let max_particles: f32 = MAX_PARTICLES as f32;

        let mut particle_pos: [f32; 2 * MAX_PARTICLES] = [0.0; 2 * MAX_PARTICLES];
        let mut particle_color: [f32; 3 * MAX_PARTICLES] = [0.0; 3 * MAX_PARTICLES];
        for i in 0..MAX_PARTICLES {
            particle_color[3 * i + 2]  = 1.0;
        }

        let particle_vel: [f32; 2 * MAX_PARTICLES] = [0.0; 2 * MAX_PARTICLES];
        let particle_density: [f32; F_NUM_CELLS] = [0.0; F_NUM_CELLS];
        let particle_rest_density: f32 =  0.0;
        let particle_radius: f32 = PARTICLE_RADIUS;
        let particle_diameter: f32 = particle_radius * 2.0;
        let p_inv_spacing: f32 = 1.0 / (2.2 * PARTICLE_RADIUS);
        let p_num_x: f32 = (WIDTH * p_inv_spacing).floor() + 1.0;
        let p_num_y: f32 = (HEIGHT * p_inv_spacing).floor() + 1.0;
        let p_num_cells: f32 = p_num_x * p_num_y;

        let num_cell_particles: [usize; P_NUM_CELLS] = [0; P_NUM_CELLS];
        let first_cell_particle: [usize; 1 + P_NUM_CELLS] = [0; 1 + P_NUM_CELLS];
        let cell_particle_ids: [usize; MAX_PARTICLES] = [0; MAX_PARTICLES];

        let num_particles: u32 = 0;

        // for i in 0..f_num_x as usize {
        //     for j in 0..f_num_y as usize {
        //         let s: f32 = 1.0;
        //         if (i == 0 || i == f_num_x - 1.0 || j == 0 || )
        //     }
        // }

        for i in 0..F_NUM_CELLS {
            let row: u32 = i as u32 / CONFIG.scene.resolution as u32;
            let column: u32 = i as u32 % CONFIG.scene.resolution as u32;

            let x: f32 = Utils::get_pixel(row as f32, SPACING);
            let y: f32 = Utils::get_pixel(column as f32, SPACING) + h;

            cell_indices[i] = i as u32;

            cell_vertices[i * 5] = x;
            cell_vertices[i * 5 + 1] = y;

            if Utils::is_border(row, column, ROWS as u32, COLUMNS as u32) {
                cell_vertices[i * 5 + 2] = 0.3;
                cell_vertices[i * 5 + 3] = 0.3;
                cell_vertices[i * 5 + 4] = 0.3;

                s[i] = 0.0;
                cell_type[i] = SOLID_CELL;
            } else {
                cell_vertices[i * 5 + 2] = 0.8;
                cell_vertices[i * 5 + 3] = 0.8;
                cell_vertices[i * 5 + 4] = 0.8;

                s[i] = 1.0;
                cell_type[i] = FLUID_CELL; // AIR_CELL
            }
        }

        let cols = ((WIDTH - (h * 2.0)) / particle_diameter).floor() as usize;

        let mut row: usize = 0;
        let mut col: usize = 0;
        for i in 0..MAX_PARTICLES {
            if col >= cols {
                col = 0;
                row += 1;
            }

            let stagger_x: f32 = if row % 2 == 0 { 0.0 } else { particle_diameter / 2.0 };

            let x: f32 = h + CONFIG.particle.radius + (col as f32 * particle_diameter) + stagger_x;
            let y: f32 = CONFIG.scene.width - (h + CONFIG.particle.radius + (row as f32 * particle_diameter));

            particle_pos[i * 2] = x;
            particle_pos[i * 2 + 1] = y;
            particle_color[i * 3] = 0.0;
            particle_color[i * 3 + 1] = 0.0;
            particle_color[i * 3 + 2] = 3.0;

            particle_indices[i] = i as u32;

            col += 1;
        }

        Fluid {
            cell_vertices,
            cell_indices,
            grid_vertices,
            grid_indices,
            particle_vertices,
            particle_indices,
            density,
            f_num_x,
            f_num_y,
            h,
            f_inv_spacing,
            f_num_cells,
            u,
            v,
            du,
            dv,
            prev_u,
            prev_v,
            p,
            s,
            cell_type,
            cell_color,
            max_particles,
            particle_pos,
            particle_color,
            particle_vel,
            particle_density,
            particle_rest_density,
            particle_radius,
            p_inv_spacing,
            p_num_x,
            p_num_y,
            p_num_cells,
            num_cell_particles,
            first_cell_particle,
            cell_particle_ids,
            num_particles,
            particle_diameter,
        }
    }

    pub fn integrate_particles(&mut self, sdt: f32, gravity: f32) {
        for i in 0..MAX_PARTICLES {
            self.particle_vel[2 * i + 1] += sdt * gravity;
            self.particle_pos[2 * i] += self.particle_vel[2 * i] * sdt;
            self.particle_pos[2 * i + 1] += self.particle_vel[2 * i + 1] * sdt;
        }
    }

    pub fn push_particles_apart(&mut self, num_iters: u32) {
        let color_diffusion_coeff: f32 = 0.001;

        self.num_cell_particles.fill(0);

        for i in 0..MAX_PARTICLES {
            let x: f32 = self.particle_pos[2 * i];
            let y: f32 = self.particle_pos[2 * i + 1];

            let xi: f32 = clamp((x * self.p_inv_spacing).floor(), 0.0, self.p_num_x - 1.0);
            let yi: f32 = clamp((y * self.p_inv_spacing).floor(), 0.0, self.p_num_y - 1.0);

            let cell_nr: usize = (xi * (self.p_num_y - 1.0) + yi) as usize;
            self.num_cell_particles[cell_nr] += 1;
        }

        let mut first: usize = 0;
        for i in 0..P_NUM_CELLS {
            first += self.num_cell_particles[i];
            self.first_cell_particle[i] = first;
        }

        self.first_cell_particle[P_NUM_CELLS] = first;

        for i in 0..MAX_PARTICLES {
            let x: f32 = self.particle_pos[2 * i];
            let y: f32 = self.particle_pos[2 * i + 1];

            let xi: f32 = clamp((x * self.p_inv_spacing).floor(), 0.0, self.p_num_x - 1.0);
            let yi: f32 = clamp((y * self.p_inv_spacing).floor(), 0.0, self.p_num_y - 1.0);

            let cell_nr: usize = (xi * (self.p_num_y - 1.0) + yi) as usize;
            self.first_cell_particle[cell_nr] -= 1;
            self.cell_particle_ids[self.first_cell_particle[cell_nr]] = i;
        }

        let min_dist: f32 = 2.0 * self.particle_radius;
        let min_dist_2: f32 = min_dist * min_dist;

        for iter in 0..num_iters {
            for i in 0..MAX_PARTICLES {
                let px: f32 = self.particle_pos[2 * i];
                let py: f32 = self.particle_pos[2 * i + 1];

                let pxi: f32 = (px * self.p_inv_spacing).floor();
                let pyi: f32=  (py * self.p_inv_spacing).floor();
                let x0: usize = (pxi - 1.0).max(0.0) as usize;
                let y0: usize = (pyi - 1.0).max(0.0) as usize;
                let x1: usize = (pxi + 1.0).min(self.p_num_x - 1.0) as usize;
                let y1: usize = (pyi + 1.0).min(self.p_num_y - 1.0) as usize;

                for xi in x0..=x1 {
                    for yi in y0..=y1 {
                        let cell_nr: usize = (xi * (self.p_num_y - 1.0) as usize + yi) as usize;

                        let first: usize = self.first_cell_particle[cell_nr];
                        let last: usize = self.first_cell_particle[cell_nr + 1];

                        for j in first..last {
                            let id: usize = self.cell_particle_ids[j];

                            if id == i {
                                continue;
                            }

                            let qx: f32 = self.particle_pos[2 * id];
                            let qy: f32 = self.particle_pos[2 * id + 1];

                            let mut dx: f32 = qx - px;
                            let mut dy: f32 = qy - py;

                            let d2: f32 = dx * dx + dy * dy;

                            if d2 > min_dist_2 || d2 == 0.0 {
                                continue;
                            }

                            let d: f32 = d2.sqrt();

                            let s: f32 = 0.5 * (min_dist - d) / d;

                            dx *= s;
                            dy *= s;

                            self.particle_pos[2 * i] -= dx;
                            self.particle_pos[2 * i + 1] -= dy;
                            self.particle_pos[2 * id] += dx;
                            self.particle_pos[2 * id + 1] += dy;

                            for k in 0..3 {
                                let color0: f32 = self.particle_color[3 * i + k];
                                let color1: f32 = self.particle_color[3 * id + k];
                                let color: f32 = (color0 + color1) * 0.5;

                                self.particle_color[3 * i + k] = color0 + (color - color0) * color_diffusion_coeff;
                                self.particle_color[3 * id + k] = color1 + (color - color1) * color_diffusion_coeff;
                            }
                        }
                    }
                }
            }
        }
    }

    pub fn handle_particle_collisions(&mut self, obstacle: &Obstacle) {
        let h: f32 = 1.0 / self.f_inv_spacing;
        let r: f32 = self.particle_radius;
        let or: f32 = obstacle.radius;
        let or2: f32 = or * or;

        let min_dist: f32 = or + r;
        let min_dist2: f32 = min_dist * min_dist;

        let min_x: f32 = h + r;
        let max_x: f32 = (self.f_num_x - 1.0) * h - r;
        let min_y: f32 = h + r;
        let max_y: f32 = (self.f_num_y - 1.0) * h - r;

        for i in 0..MAX_PARTICLES {
            let mut x: f32 = self.particle_pos[2 * i];
            let mut y: f32 = self.particle_pos[2 * i + 1];

            let dx: f32 = x - obstacle.position.x;
            let dy: f32 = y - obstacle.position.y;
            let d2: f32 = dx * dx + dy * dy;

            if d2 < min_dist2 {
                self.particle_vel[2 * i] = obstacle.velocity.u;
                self.particle_vel[2 * i + 1] = obstacle.velocity.v;
            }

            if x < min_x {
                x = min_x;
                self.particle_vel[2 * i] = 0.0;
            }

            if x > max_x {
                x = max_x;
                self.particle_vel[2 * i] = 0.0;
            }

            if y < min_y {
                y = min_y;
                self.particle_vel[2 * i + 1] = 0.0;
            }

            if y > max_y {
                y = max_y;
                self.particle_vel[2 * i + 1] = 0.0;
            }

            self.particle_pos[2 * i] = x;
            self.particle_pos[2 * i + 1] = y;
        }
    }

    pub fn transfer_velocities(&mut self, to_grid: bool, flip_ratio: f32) {
        let n: usize = (self.f_num_y - 1.0) as usize;
        let h: f32 = self.h;
        let h1: f32 = self.f_inv_spacing;
        let h2: f32 = 0.5 * h;

        if to_grid {
            self.prev_u.copy_from_slice(&self.u);
            self.prev_v.copy_from_slice(&self.v);

            self.du.fill(0.0);
            self.dv.fill(0.0);

            self.u.fill(0.0);
            self.v.fill(0.0);

            for i in 0..F_NUM_CELLS {
                self.cell_type[i] = if self.s[i] == 0.0 { SOLID_CELL } else { AIR_CELL };
            }

            for i in 0..MAX_PARTICLES {
                let x: f32 = self.particle_pos[2 * i];
                let y: f32 = self.particle_pos[2 * i + 1];

                let xi: usize = clamp((x * h1).floor(), 0.0, self.f_num_x) as usize;
                let yi: usize = clamp((y * h1).floor(), 0.0, self.f_num_y) as usize;

                let cell_nr: usize = xi * n + yi;

                if self.cell_type[cell_nr] == AIR_CELL {
                    self.cell_type[cell_nr] = FLUID_CELL;
                }
            }
        }

        for component in 0..2 {
            let dx: f32 = if component == 0 { 0.0 } else { h2 };
            let dy: f32 = if component == 0 { h2 } else { 0.0 };

            let (f, prev_f, d) = match component {
                0 => (&mut self.u, &self.prev_u, &mut self.du),
                1 => (&mut self.v, &self.prev_v, &mut self.dv),
                _ => unreachable!(),
            };

            for i in 0..MAX_PARTICLES {
                let mut x: f32 = self.particle_pos[2 * i];
                let mut y: f32 = self.particle_pos[2 * i + 1];

                x = clamp(x, h, (self.f_num_x - 1.0) * h);
                y = clamp(y, h, (self.f_num_y - 1.0) * h);

                let x0: usize = ((x - dx) * h1).floor().min(self.f_num_x) as usize;
                let tx: f32 = ((x - dx) - x0 as f32 * h) * h1;
                let x1: usize = (x0 + 1).min((self.f_num_x - 2.0) as usize);

                let y0: usize = ((y - dy) * h1).floor().min(self.f_num_y) as usize;
                let ty: f32 = ((y - dy) - y0 as f32 * h) * h1;
                let y1: usize = (y0 + 1).min((self.f_num_y) as usize);

                let sx: f32 = 1.0 - tx;
                let sy: f32 = 1.0 - ty;

                let d0: f32 = sx * sy;
                let d1: f32 = tx * sy;
                let d2: f32 = tx * ty;
                let d3: f32 = sx * ty;

                let nr0: usize = x0 * n + y0;
                let nr1: usize = x1 * n + y0;
                let nr2: usize = x1 * n + y1;
                let nr3: usize = x0 * n + y1;

                if to_grid {
                    let pv = self.particle_vel[2 * i + component];

                    f[nr0] += pv * d0; d[nr0] += d0;
                    f[nr1] += pv * d1; d[nr1] += d1;
                    f[nr2] += pv * d2; d[nr2] += d2;
                    f[nr3] += pv * d3; d[nr3] += d3;
                } else {
                    let offset: usize = if component == 0 { n as usize } else { 1 };

                    let valid0: f32 = if self.cell_type[nr0] != AIR_CELL || self.cell_type[nr0 - offset] != AIR_CELL { 1.0 } else { 0.0 };
                    let valid1: f32 = if self.cell_type[nr1] != AIR_CELL || self.cell_type[nr1 - offset] != AIR_CELL { 1.0 } else { 0.0 };
                    let valid2: f32 = if self.cell_type[nr2] != AIR_CELL || self.cell_type[nr2 - offset] != AIR_CELL { 1.0 } else { 0.0 };
                    let valid3: f32 = if self.cell_type[nr3] != AIR_CELL || self.cell_type[nr3 - offset] != AIR_CELL { 1.0 } else { 0.0 };

                    let v: f32 = self.particle_vel[2 * i + component];
                    let d: f32 = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

                    if d > 0.0 {
                        let pic_v = (valid0 * d0 *f[nr0] + valid1 * d1 * f[nr1] + valid2 * d2 * f[nr2] + valid3 * d3 * f[nr3]) / d;
                        let corr = (valid0 * d0 * (f[nr0] - prev_f[nr0]) + valid1 * d1 * (f[nr1] - prev_f[nr1]) + valid2 * d2 * (f[nr2] - prev_f[nr2]) + valid3 * d3 * (f[nr3] - prev_f[nr3])) / d;
                        let flip_v = v + corr;

                        self.particle_vel[2 * i + component] = (1.0 - flip_ratio) * pic_v + flip_ratio * flip_v;
                    }
                }
            }

            if to_grid {
                for i in 0..f.len() {
                    if d[i] > 0.0 {
                        f[i] /= d[i];
                    }
                }

                let f_num_x: usize = (self.f_num_x) as usize;
                let f_num_y: usize = (self.f_num_y) as usize;
                for i in 0..f_num_x {
                    for j in 0..f_num_y {
                        let solid = self.cell_type[i * n + j] == SOLID_CELL;

                        if solid || (i > 0 && self.cell_type[(i - 1) * n + j] == SOLID_CELL) {
                            self.u[i * n + j] = self.prev_u[i * n + j];
                        }

                        if solid || (j > 0 && self.cell_type[i * n + j - 1] == SOLID_CELL) {
                            self.v[i * n + j] = self.prev_v[i * n + j];
                        }
                    }
                }
            }
        }
    }

    pub fn update_particle_density(&mut self) {
        let n: usize = (self.f_num_y - 1.0) as usize;
        let h: f32 = self.h;
        let h1: f32 = self.f_inv_spacing;
        let h2: f32 = 0.5 * h;

        self.particle_density.fill(0.0);

        for i in 0..MAX_PARTICLES {
            let mut x: f32 = self.particle_pos[2 * i];
            let mut y: f32 = self.particle_pos[2 * i + 1];

            x = clamp(x, h, (self.f_num_x - 1.0) * h);
            y = clamp(y, h, (self.f_num_y - 1.0) * h);

            let x0: usize = ((x - h2) * h1).floor() as usize;
            let tx: f32 = ((x - h2) - x0 as f32 * h) * h1;
            let x1: usize = (x0 + 1).min(n);

            let y0: usize = ((y - h2) * h1).floor() as usize;
            let ty: f32 = ((y - h2) - y0 as f32 * h) * h1;
            let y1: usize = (y0 + 1).min(n);

            let sx: f32 = 1.0 - tx;
            let sy: f32 = 1.0 - ty;

            if x0 < n && y0 < n {
                self.particle_density[x0 * n + y0] += sx * sy;
            }
            if x1 < n && y0 < n {
                self.particle_density[x1 * n + y0] += tx * sy;
            }
            if x1 < n && y1 < n {
                self.particle_density[x1 * n + y1] += tx * ty;
            }
            if x0 < n && y1 < n {
                self.particle_density[x0 * n + y1] += sx * ty;
            }
        }

        if self.particle_rest_density == 0.0 {
            let mut sum: f32 = 0.0;
            let mut num_fluid_cells: f32 = 0.0;

            for i in 0..F_NUM_CELLS {
                if self.cell_type[i] == FLUID_CELL { // FLUID_CELL
                    sum += self.particle_density[i];
                    num_fluid_cells += 1.0;
                }
            }

            if num_fluid_cells > 0.0 {
                self.particle_rest_density = sum / num_fluid_cells;
            }
        }
    }

    pub fn solve_incompressibility(&mut self, num_iters: u32, dt: f32, over_relaxation: f32, compensate_drift: bool) {
        self.p.fill(0.0);

        self.prev_u.copy_from_slice(&self.u);
        self.prev_v.copy_from_slice(&self.v);

        let n: usize = (self.f_num_y - 1.0) as usize;
        let cp = self.density * self.h / dt;

        let f_num_x: usize = (self.f_num_x) as usize;
        let f_num_y: usize = (self.f_num_y) as usize;

        for iter in 0..num_iters {
            for i in 1..f_num_x {
                for j in 1..f_num_y {
                    if self.cell_type[i * n + j] != FLUID_CELL {
                        continue;
                    }

                    let center: usize = i * n + j;    
                    let left: usize = (i - 1) * n + j;
                    let right: usize = (i + 1) * n + j;
                    let bottom: usize = i * n + j - 1;
                    let top: usize = i * n + j + 1;

                    let sx0: f32 = self.s[left];
                    let sx1: f32 = self.s[right];
                    let sy0: f32 = self.s[bottom];
                    let sy1: f32 = self.s[top];

                    let s: f32 = sx0 + sx1 + sy0 + sy1;

                    if s == 0.0 {
                        continue;
                    }

                    let mut div: f32 = self.u[right] - self.u[center] + self.v[top] - self.v[center];

                    if self.particle_rest_density > 0.0 && compensate_drift {
                        let k: f32 = 1.0;
                        let compression = self.particle_density[i * n + j] - self.particle_rest_density;
                        if compression > 0.0 {
                            div = div - k * compression;
                        }
                    }

                    let mut p: f32 = -div / s;
                    p *= over_relaxation;
                    self.p[center] += cp * p;

                    self.u[center] -= sx0 * p;
                    self.u[right] += sx1 * p;
                    self.v[center] -= sy0 * p;
                    self.v[top] += sy1 * p;
                }
            }
        }
    }

    pub fn update_particle_colors(&mut self) {
        let h1: f32 = self.f_inv_spacing;

        for i in 0..MAX_PARTICLES {
            let s: f32 = 0.01;

            self.particle_color[3 * i] = clamp(self.particle_color[3 * i] - s, 0.0, 1.0);
            self.particle_color[3 * i + 1] = clamp(self.particle_color[3 * i + 1] - s, 0.0, 1.0);
            self.particle_color[3 * i + 2] = clamp(self.particle_color[3 * i + 2] - s, 0.0, 1.0);

            let x: f32 = self.particle_pos[2 * i];
            let y: f32 = self.particle_pos[2 * i + 1];

            let xi: f32 = clamp((x * h1).floor(), 1.0, self.f_num_x - 1.0);
            let yi: f32 = clamp((y * h1).floor(), 1.0, self.f_num_y - 1.0);
            let cell_nr: usize = (xi * (self.f_num_y - 1.0) + yi) as usize;

            let d0: f32 = self.particle_rest_density;

            if d0 > 0.0 {
                let rel_density: f32 = self.particle_density[cell_nr] / d0;
                if rel_density > 0.7 {
                    let s: f32 = 0.8;
                    self.particle_color[3 * i] = s;
                    self.particle_color[3 * i + 1] = s;
                    self.particle_color[3 * i + 2] = 1.0;
                }
            }

        }
    }

    pub fn update_cell_colors(&mut self) {
        self.cell_color.fill(0.0);

        for i in 0..F_NUM_CELLS {
            if self.cell_type[i] == SOLID_CELL {
                self.cell_color[3 * i] = 0.5;
                self.cell_color[3 * i + 1] = 0.5;
                self.cell_color[3 * i + 2] = 0.5;
            } else if self.cell_type[i] == FLUID_CELL {
                if self.particle_rest_density > 0.0 {
                    self.particle_density[i] /= self.particle_rest_density;
                }

                let min_val: f32 = 0.0;
                let max_val: f32 = 2.0;

                self.particle_density[i] = (self.particle_density[i].max(min_val)).min(max_val - 0.0001);

                let d = max_val - min_val;
                self.particle_density[i] = if d == 0.0 { 0.5 } else { (self.particle_density[i] - min_val) / d };

                let m = 0.25;
                let num = (self.particle_density[i] / m).floor() as u32;
                let s = (self.particle_density[i] - (num as f32) * m) / m;

                let (r, g, b) = match num {
                    0 => (0.0, s, 1.0),
                    1 => (0.0, 1.0, 1.0 - s),
                    2 => (s, 1.0, 0.0),
                    3 => (1.0, 1.0 - s, 0.0),
                    _ => (1.0, 0.0, 0.0), // Fallback for unexpected input
                };

                self.cell_color[3 * i] = r;
                self.cell_color[3 * i + 1] = g;
                self.cell_color[3 * i + 2] = b;
            }
        }
    }
}
