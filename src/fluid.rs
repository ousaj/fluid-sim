use crate::utils::Utils;
use crate::obstacle::Obstacle;

use crate::particle::Particle;
use crate::position::Position;
use crate::cell::Cell;
use crate::cell::CellType;
use crate::color::Color;

use crate::config::FLIP_RATIO;
use crate::config::GRAVITY;
use crate::config::IS_VELOCITY_MAPPED;
use crate::config::OVER_RELAXATION;
use crate::config::PARTICLE_ITERATIONS;
use crate::config::PRESSURE_ITERATIONS;
use crate::config::TOTAL_PARTICLES;
use crate::config::WINDOW_SIZE;
use crate::config::VISCOSITY;
use crate::config::PARTICLE_RADIUS;
use crate::config::DENSITY;
use crate::config::TIME_STEP;
use crate::config::OBSTACLE_RADIUS;
use crate::config::RESOLUTION;
use crate::config::CELL_SIZE;
use crate::config::TOTAL_CELLS;
use crate::config::PARTICLE_VELOCITY_RANGE;
use crate::config::CELL_DENSITY_RANGE;

use rayon::prelude::*;

const HALF_CELL_SIZE: f32 = CELL_SIZE * 0.5;
const INVERSE_CELL_SIZE: f32 = 1.0 / CELL_SIZE;
const MIN_X: f32 = CELL_SIZE + PARTICLE_RADIUS;
const MIN_Y: f32 = CELL_SIZE + PARTICLE_RADIUS;
const MAX_X: f32 = WINDOW_SIZE - (CELL_SIZE + PARTICLE_RADIUS);
const MAX_Y: f32 = WINDOW_SIZE - (CELL_SIZE + PARTICLE_RADIUS);
const DRIFT_COMPENSATION_FACTOR: f32 = 2.2;
const PRESSURE_SCALING_FACTOR: f32 = DENSITY * CELL_SIZE / TIME_STEP;

const MIN_DIST: f32 = 2.0 * PARTICLE_RADIUS;
const MIN_DIST_2: f32 = MIN_DIST * MIN_DIST;

const OBSTACLE_MIN_DIST: f32 = OBSTACLE_RADIUS + PARTICLE_RADIUS;
const OBSTACLE_MIN_DIST_2: f32 = OBSTACLE_MIN_DIST * OBSTACLE_MIN_DIST;
// const COLOR_DIFFUSION_COEFF: f32 = 0.001;

const PARTICLE_DIAMETER: f32 = PARTICLE_RADIUS * 2.0;

pub struct Fluid {
    pub num_cell_particles: [usize; TOTAL_CELLS],
    pub first_cell_particle: [usize; TOTAL_CELLS + 1],
    pub cell_particle_ids: [usize; TOTAL_PARTICLES],
    pub is_velocity_mapped: bool,
    pub particles: [Particle; TOTAL_PARTICLES],
    pub cells: [Cell; TOTAL_CELLS],
    pub rest_density: f32,
}

impl Fluid {
    pub fn new() -> Self {
        let resolution: usize = RESOLUTION as usize;
        let cell_size: f32 = WINDOW_SIZE / RESOLUTION as f32;
        let num_cell_particles: [usize; TOTAL_CELLS] = [0; TOTAL_CELLS];

        let cell: Cell = Cell::new(Position::new(0.0, 0.0), Color::new(0.0, 0.0, 0.0), CellType::SOLID);
        let mut cells: [Cell; TOTAL_CELLS] = [cell; TOTAL_CELLS];
        cells
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, cell)| {
                let row: usize = i / resolution;
                let column: usize = i % resolution;

                let x: f32 = Utils::get_pixel(row, cell_size);
                let y: f32 = Utils::get_pixel(column,cell_size) + cell_size;

                cell.position.x = x;
                cell.position.y = y;

                if Utils::is_border(row, column, resolution) {
                    cell.color = Color::border_color();
                    cell.solid_coeff = 0.0;
                    cell.cell_type = CellType::SOLID;
                } else {
                    cell.color = Color::background_color();
                    cell.solid_coeff = 1.0;
                    cell.cell_type = CellType::AIR;
                }
                
            });

        let particle_diameter: f32 = PARTICLE_RADIUS * 2.0;
        let particles_per_row: usize = (((WINDOW_SIZE - (cell_size * 2.0)) / particle_diameter).floor() - 1.0) as usize;

        let first_cell_particle: [usize; TOTAL_CELLS + 1] = [0; TOTAL_CELLS + 1];
        let cell_particle_ids: [usize; TOTAL_PARTICLES] = [0; TOTAL_PARTICLES];
        let is_velocity_mapped: bool = IS_VELOCITY_MAPPED;
        let rest_density: f32 = 0.0;

        let mut particles: [Particle; TOTAL_PARTICLES] = [Particle::new(Position::new(0.0, 0.0)); TOTAL_PARTICLES];
        particles
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, particle)| {
                let row: usize = i / particles_per_row;
                let column: usize = i % particles_per_row;

                let stagger_x: f32 = if row % 2 == 0 { 0.0 } else { PARTICLE_RADIUS };
                let x: f32 = CELL_SIZE + PARTICLE_RADIUS + (column as f32 * PARTICLE_DIAMETER) + stagger_x;
                let y: f32 = WINDOW_SIZE - (CELL_SIZE + PARTICLE_RADIUS + (row as f32 * PARTICLE_DIAMETER));

                particle.position.x = x;
                particle.position.y = y;
                particle.color.r = 0.0;
                particle.color.g = 0.0;
                particle.color.b = 1.0;
            });

        Fluid {
            num_cell_particles,
            first_cell_particle,
            cell_particle_ids,
            is_velocity_mapped,
            particles,
            cells,
            rest_density,
        }
    }

    pub fn integrate_particles(&mut self) {
        self.particles
            .par_iter_mut()
            .for_each(|particle| {
                particle.velocity.v += TIME_STEP * GRAVITY * (1.0 + 40.0 * VISCOSITY);
                particle.position.x += particle.velocity.u * TIME_STEP;
                particle.position.y += particle.velocity.v * TIME_STEP;
            });
    }

    pub fn push_particles_apart(&mut self) {
        self.num_cell_particles.fill(0);

        self.particles
            .iter()
            .for_each(|particle| {
                let x: f32 = particle.position.x;
                let y: f32 = particle.position.y;

                let xi: usize = Utils::clamp_usize((x * INVERSE_CELL_SIZE).floor() as usize, 1, RESOLUTION - 2);
                let yi: usize = Utils::clamp_usize((y * INVERSE_CELL_SIZE).floor() as usize, 1, RESOLUTION - 2);

                let cell_nr: usize = xi * RESOLUTION + yi;
                self.num_cell_particles[cell_nr] += 1;
            });

        let mut first: usize = 0;
        self.cells
            .iter()
            .enumerate()
            .for_each(|(i, cell)| {
                first += self.num_cell_particles[i];
                self.first_cell_particle[i] = first;
            });

        self.first_cell_particle[TOTAL_CELLS] = first;

        self.particles
            .iter()
            .enumerate()
            .for_each(|(i, particle)| {
                let x: f32 = particle.position.x;
                let y: f32 = particle.position.y;

                let xi: usize = Utils::clamp_usize((x * INVERSE_CELL_SIZE).floor() as usize, 1, RESOLUTION - 2);
                let yi: usize = Utils::clamp_usize((y * INVERSE_CELL_SIZE).floor() as usize, 1, RESOLUTION - 2);

                let cell_nr: usize = xi * RESOLUTION + yi;
                self.first_cell_particle[cell_nr] -= 1;
                self.cell_particle_ids[self.first_cell_particle[cell_nr]] = i;
            });

        // Smoothing iterations.
        for _iter in 0..PARTICLE_ITERATIONS {
            // Para cada particula.
            for i in 0..TOTAL_PARTICLES {
                let px: f32 = self.particles[i].position.x;
                let py: f32 = self.particles[i].position.y;

                let pxi: usize = (px * INVERSE_CELL_SIZE).floor() as usize;
                let pyi: usize = (py * INVERSE_CELL_SIZE).floor() as usize;

                // TODO: Puede que aquí puedas aplicar el límite 1,1 23,23.
                let x0: usize = pxi.saturating_sub(1).max(1);
                let y0: usize = pyi.saturating_sub(1).max(1);
                let x1: usize = (pxi + 1).min(RESOLUTION - 1);
                let y1: usize = (pyi + 1).min(RESOLUTION - 1);

                // Las 9 celdas vecinas.
                for xi in x0..=x1 {
                    for yi in y0..=y1 {
                        let cell_nr: usize = xi * RESOLUTION + yi;

                        let first: usize = self.first_cell_particle[cell_nr];
                        let last: usize = self.first_cell_particle[cell_nr + 1];

                        if first == last {
                            continue;
                        }

                        // Cada particula en esa celda.
                        for j in first..last {
                            let id: usize = self.cell_particle_ids[j];

                            if id == i {
                                continue;
                            }

                            let qx: f32 = self.particles[id].position.x;
                            let qy: f32 = self.particles[id].position.y;


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

                            self.particles[i].position.x -= dx;
                            self.particles[i].position.y -= dy;

                            self.particles[id].position.x += dx;
                            self.particles[id].position.y += dy;
                        }
                    }
                }
            }
        }
    }

    pub fn handle_particle_collisions(&mut self, obstacle: &Obstacle) {
        self.particles
            .par_iter_mut()
            .for_each(|particle| {
                let mut x: f32 = particle.position.x;
                let mut y: f32 = particle.position.y;

                let dx: f32 = x - obstacle.position.x;
                let dy: f32 = y - obstacle.position.y;
                let d2: f32 = dx * dx + dy * dy;

                if obstacle.is_enabled && d2 < OBSTACLE_MIN_DIST_2 {
                    particle.velocity.u = obstacle.velocity.u;
                    particle.velocity.v = obstacle.velocity.v;
                }

                if x < MIN_X {
                    x = MIN_X;
                    particle.velocity.u = 0.0;
                }

                if x > MAX_X {
                    x = MAX_X;
                    particle.velocity.u = 0.0;
                }

                if y < MIN_Y {
                    y = MIN_Y;
                    particle.velocity.v = 0.0;
                }

                if y > MAX_Y {
                    y = MAX_Y;
                    particle.velocity.v = 0.0;
                }

                particle.position.x = x;
                particle.position.y = y;
            });
    }

    pub fn transfer_velocities(&mut self, to_grid: bool) {
        if to_grid {
            self.cells
                .par_iter_mut()
                .for_each(|cell| {
                    cell.previous_velocity = cell.velocity;
                    cell.delta.u = 0.0;
                    cell.delta.v = 0.0;
                    cell.velocity.u = 0.0;
                    cell.velocity.v = 0.0;

                    cell.cell_type = if cell.solid_coeff == 0.0 { CellType::SOLID } else { CellType::AIR };
                });

            self.particles
                .iter()
                .for_each(|particle| {
                    let x: f32 = particle.position.x;
                    let y: f32 = particle.position.y;

                    let xi: usize = Utils::clamp_usize((x * INVERSE_CELL_SIZE).floor() as usize, 1 ,RESOLUTION - 2);
                    let yi: usize = Utils::clamp_usize((y * INVERSE_CELL_SIZE).floor() as usize, 1, RESOLUTION - 2);
                    let cell_nr: usize = xi * RESOLUTION + yi;

                    // Pasamos a fluidas las celdas con partículas.
                    if self.cells[cell_nr].cell_type == CellType::AIR {
                        self.cells[cell_nr].cell_type = CellType::FLUID;
                    }
                });
        }

        for component in 0..2 {
            let dx: f32 = if component == 0 { 0.0 } else { HALF_CELL_SIZE };
            let dy: f32 = if component == 0 { HALF_CELL_SIZE } else { 0.0 };

            for i in 0..TOTAL_PARTICLES {
                let mut x: f32 = self.particles[i].position.x;
                let mut y: f32 = self.particles[i].position.y;

                // Mantenemos la partícula dentro del dominio.
                x = Utils::clamp(x, MIN_X, MAX_X);
                y = Utils::clamp(y, MIN_Y, MAX_Y);

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
                    if component == 0 {
                        self.cells[nr0].velocity.u += self.particles[i].velocity.u * d0;
                        self.cells[nr1].velocity.u += self.particles[i].velocity.u * d1;
                        self.cells[nr2].velocity.u += self.particles[i].velocity.u * d2;
                        self.cells[nr3].velocity.u += self.particles[i].velocity.u * d3;

                        self.cells[nr0].delta.u += d0;
                        self.cells[nr1].delta.u += d1;
                        self.cells[nr2].delta.u += d2;
                        self.cells[nr3].delta.u += d3;
                    } else {
                        self.cells[nr0].velocity.v += self.particles[i].velocity.v * d0;
                        self.cells[nr1].velocity.v += self.particles[i].velocity.v * d1;
                        self.cells[nr2].velocity.v += self.particles[i].velocity.v * d2;
                        self.cells[nr3].velocity.v += self.particles[i].velocity.v * d3;

                        self.cells[nr0].delta.v += d0;
                        self.cells[nr1].delta.v += d1;
                        self.cells[nr2].delta.v += d2;
                        self.cells[nr3].delta.v += d3;
                    }
                } else { // To particle.
                    let offset: usize = if component == 0 { RESOLUTION } else { 1 };
                    let valid0: f32 = if self.cells[nr0].cell_type != CellType::AIR || self.cells[nr0 - offset].cell_type != CellType::AIR { 1.0 } else { 0.0 };
                    let valid1: f32 = if self.cells[nr1].cell_type != CellType::AIR || self.cells[nr1 - offset].cell_type != CellType::AIR { 1.0 } else { 0.0 };
                    let valid2: f32 = if self.cells[nr2].cell_type != CellType::AIR || self.cells[nr2 - offset].cell_type != CellType::AIR { 1.0 } else { 0.0 };
                    let valid3: f32 = if self.cells[nr3].cell_type != CellType::AIR || self.cells[nr3 - offset].cell_type != CellType::AIR { 1.0 } else { 0.0 };

                    let v: f32 = if component == 0 { self.particles[i].velocity.u } else { self.particles[i].velocity.v };
                    let d: f32 = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

                    if d > 0.0 {
                        let pic_v: f32;
                        let corr: f32;
                        let flip_v: f32;

                        if component == 0 {
                            pic_v = (
                                valid0 * d0 * self.cells[nr0].velocity.u 
                              + valid1 * d1 * self.cells[nr1].velocity.u 
                              + valid2 * d2 * self.cells[nr2].velocity.u 
                              + valid3 * d3 * self.cells[nr3].velocity.u
                            ) / d;

                            corr = (
                                valid0 * d0 * (self.cells[nr0].velocity.u - self.cells[nr0].previous_velocity.u) 
                              + valid1 * d1 * (self.cells[nr1].velocity.u - self.cells[nr1].previous_velocity.u) 
                              + valid2 * d2 * (self.cells[nr2].velocity.u - self.cells[nr2].previous_velocity.u) 
                              + valid3 * d3 * (self.cells[nr3].velocity.u - self.cells[nr3].previous_velocity.u)
                            ) / d;

                            flip_v = v + corr;
                            self.particles[i].velocity.u = ((1.0 - FLIP_RATIO) * pic_v + FLIP_RATIO * flip_v) * (1.0 - VISCOSITY); 
                        } else {
                            pic_v = (
                                valid0 * d0 * self.cells[nr0].velocity.v 
                              + valid1 * d1 * self.cells[nr1].velocity.v 
                              + valid2 * d2 * self.cells[nr2].velocity.v 
                              + valid3 * d3 * self.cells[nr3].velocity.v
                            ) / d;

                            corr = (
                                valid0 * d0 * (self.cells[nr0].velocity.v - self.cells[nr0].previous_velocity.v) 
                              + valid1 * d1 * (self.cells[nr1].velocity.v - self.cells[nr1].previous_velocity.v) 
                              + valid2 * d2 * (self.cells[nr2].velocity.v - self.cells[nr2].previous_velocity.v) 
                              + valid3 * d3 * (self.cells[nr3].velocity.v - self.cells[nr3].previous_velocity.v)
                            ) / d;

                            flip_v = v + corr;
                            self.particles[i].velocity.v = ((1.0 - FLIP_RATIO) * pic_v + FLIP_RATIO * flip_v) * (1.0 - VISCOSITY); 
                        }
                    }
                }
            }

        }

        if to_grid {
            self.cells
                .par_iter_mut()
                .for_each(|cell| {
                    if cell.delta.u > 0.0 {
                        cell.velocity.u /= cell.delta.u;
                    }

                    if cell.delta.v > 0.0 {
                        cell.velocity.v /= cell.delta.v;
                    }
                });
        }
    }

    pub fn update_density(&mut self) {
        self.cells
            .par_iter_mut()
            .for_each(|cell| {
                cell.density = 0.0;
            });

        // TODO: Paralelizar.
        self.particles
            .iter()
            .for_each(|particle| {
                let mut x: f32 = particle.position.x;
                let mut y: f32 = particle.position.y;

                x = Utils::clamp(x, MIN_X, MAX_X);
                y = Utils::clamp(y, MIN_Y, MAX_Y);

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

                // Bottom-left.
                if x0 < RESOLUTION && y0 < RESOLUTION {
                    self.cells[x0 * RESOLUTION + y0].density += (sx * sy);
                }

                // Bottom-right.
                if x1 < RESOLUTION && y0 < RESOLUTION {
                    self.cells[x1 * RESOLUTION + y0].density += tx * sy;
                }

                // Top-right.
                if x1 < RESOLUTION && y1 < RESOLUTION {
                    self.cells[x1 * RESOLUTION + y1].density += tx * ty;
                }

                // Top-left.
                if x0 < RESOLUTION && y1 < RESOLUTION {
                    self.cells[x0 * RESOLUTION + y1].density += sx * ty;
                }
            });

        if self.rest_density == 0.0 {
            let (density_sum, fluid_cells) = self.cells
                .par_iter()
                .filter_map(|cell| {
                    if cell.cell_type == CellType::FLUID {
                        Some((cell.density, 1.0))
                    } else {
                        None
                    }
                })
                .reduce(
                    || (0.0, 0.0),
                    |(d1, c1), (d2, c2)| (d1 + d2, c1 + c2),
                );

            if fluid_cells > 0.0 {
                self.rest_density = density_sum / fluid_cells;
            }
        }
    }

    // // Corrección de presión.
    pub fn solve_incompressibility(&mut self) {
        self.cells
            .par_iter_mut()
            .for_each(|cell| {
                cell.pressure = 0.0;
                cell.previous_velocity = cell.velocity;
            });

        for _iter in 0..PRESSURE_ITERATIONS {
            for cell in 0..TOTAL_CELLS {
                let i: usize = cell / RESOLUTION;
                let j: usize = cell % RESOLUTION;
                if self.cells[i * RESOLUTION + j].cell_type != CellType::FLUID {
                    continue;
                }

                let cell: usize = i * RESOLUTION + j;
                let left_cell: usize = (i - 1) * RESOLUTION + j;
                let right_cell: usize = (i + 1) * RESOLUTION + j;
                let bottom_cell: usize = i * RESOLUTION + (j - 1);
                let top_cell: usize = i * RESOLUTION + j + 1;

                let left_cell_solid_coeff: f32 = self.cells[left_cell].solid_coeff;
                let right_cell_solid_coeff: f32 = self.cells[right_cell].solid_coeff;
                let bottom_cell_solid_coeff: f32 = self.cells[bottom_cell].solid_coeff;
                let top_cell_solid_coeff: f32 = self.cells[top_cell].solid_coeff;

                let total_solid_coeff: f32 = left_cell_solid_coeff + right_cell_solid_coeff + bottom_cell_solid_coeff + top_cell_solid_coeff;

                if total_solid_coeff == 0.0 {
                    continue;
                }

                // Cálculo de divergencia.
                // Cuánto fluido entra o sale de una celda.
                // Se mide con las diferencias de velocidad entre caras.
                let mut divergence: f32 = 
                      self.cells[right_cell].velocity.u
                    - self.cells[cell].velocity.u
                    + self.cells[top_cell].velocity.v
                    - self.cells[cell].velocity.v;
                
                // Corrección de divergencia.
                // Esto ayuda a que las partículas no se amontonen.
                // Drift control.
                if self.rest_density > 0.0 {
                    let compression = self.cells[i * RESOLUTION + j].density - self.rest_density;
                    if compression > 0.0 {
                        divergence = divergence - DRIFT_COMPENSATION_FACTOR * compression;
                    }
                }

                // Cálculo de la corrección de presión.
                // Se calcula cuánto corregir la presión.
                // OVER_RELAXATION acelera la convergencia.
                let mut pressure: f32 = -divergence / total_solid_coeff;
                pressure *= OVER_RELAXATION;
                self.cells[cell].pressure += pressure * PRESSURE_SCALING_FACTOR;

                self.cells[cell].velocity.u -= left_cell_solid_coeff * pressure;
                self.cells[right_cell].velocity.u += right_cell_solid_coeff * pressure;
                self.cells[cell].velocity.v -= bottom_cell_solid_coeff * pressure;
                self.cells[top_cell].velocity.v += top_cell_solid_coeff * pressure;
                // }
            }
        }
    }

    pub fn update_particle_colors(&mut self) {
        self.particles
            .par_iter_mut()
            .for_each(|particle| {
                if !self.is_velocity_mapped {
                    particle.color = Color::particle_color();
                } else {
                    let v: f32 = particle.velocity.u.abs() + particle.velocity.v.abs();
                    let norm: f32 = (v - PARTICLE_VELOCITY_RANGE[0]) / (PARTICLE_VELOCITY_RANGE[1] - PARTICLE_VELOCITY_RANGE[0]);

                    particle.color.r = norm;
                    particle.color.g = 0.8;
                    particle.color.b = 1.0 - norm;
                }
            });
    }

    pub fn update_cell_colors(&mut self) {
        self.cells
            .par_iter_mut()
            .for_each(|cell| {
                match cell.cell_type {
                    CellType::SOLID => cell.color = Color::border_color(),
                    CellType::AIR => cell.color = Color::background_color(),
                    CellType::FLUID => {
                        let d: f32 = cell.density.min(CELL_DENSITY_RANGE[1]);
                        let norm: f32 = (d - CELL_DENSITY_RANGE[0]) / (CELL_DENSITY_RANGE[1] - CELL_DENSITY_RANGE[0]);

                        cell.color.r = norm;
                        cell.color.g = 0.0;
                        cell.color.b = 1.0 - norm;
                    },
                };
            });
    }
}
