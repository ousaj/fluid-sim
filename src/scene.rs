use crate::windsl::Windsl;
use crate::objects::*;
use crate::fluid::Fluid;
use crate::config::*;
use crate::obstacle::Obstacle;
use crate::position::Position;
use crate::velocity::Velocity;

use rayon::prelude::*;
const WIDTH: f32 = CONFIG.scene.width;
const HEIGHT: f32 = CONFIG.scene.width;
const SPACING: f32 = HEIGHT / CONFIG.scene.resolution; // Grid cell size.

const F_NUM_X: f32 = (WIDTH / SPACING) + 1.0;
const F_NUM_Y: f32 = (HEIGHT / SPACING) + 1.0;
const F_NUM_CELLS: usize = (F_NUM_X * F_NUM_Y) as usize;

pub struct Scene {
    pub windsl: Windsl,
    pub frame_rate: f32,
    obstacle: Obstacle,
    fluid: Fluid,
    vbo: Vbo,
    vao: Vao,
    ibo: Ibo,
    program: Program,
}

impl Scene {
    pub fn new() -> Self {
        let frame_rate: f32 = CONFIG.scene.frame_rate;
        let obstacle: Obstacle = Obstacle::new(Position::new(0.0, 0.0), CONFIG.environment.obstacle_radius, Velocity::new(0.0, 0.0));
        let fluid: Fluid = Fluid::new();
        let windsl = Windsl::new(CONFIG.scene.width as usize, CONFIG.scene.width as usize).unwrap();
        let program = create_program().unwrap();

        program.set();
        program.set_resolution(CONFIG.scene.width, CONFIG.scene.width);

        unsafe { 
            gl::Viewport(0, 0, CONFIG.scene.width as i32, CONFIG.scene.width as i32);
            gl::Enable(gl::PROGRAM_POINT_SIZE);
        }


        let vbo = Vbo::gen();
        vbo.set(&Vec::new());

        let vao = Vao::gen();
        vao.set();

        let ibo = Ibo::gen();
        ibo.set(&Vec::new());

        Scene {
            obstacle,
            fluid,
            windsl,
            vbo,
            vao,
            ibo,
            program,
            frame_rate
        }
    }

    fn clear_window(&mut self) {
        unsafe {
            gl::ClearColor(0.0, 0.0, 0.0, 1.0);
            gl::Clear(gl::COLOR_BUFFER_BIT);
        }
    }

    fn draw_grid(&mut self) {
        // self.fluid.cell_vertices
        //     .par_chunks_mut(5)
        //     .zip(self.fluid.cells.par_iter())
        //     .for_each(|(chunk, cell)| {
        //         chunk[0] = cell.position.x;
        //         chunk[1] = cell.position.y;
        //         chunk[2] = cell.color.r;
        //         chunk[3] = cell.color.g;
        //         chunk[4] = cell.color.b;
        //     });

        for i in 0..F_NUM_CELLS {
            self.fluid.cell_vertices[5 * i + 2] = self.fluid.cell_color[3 * i];
            self.fluid.cell_vertices[5 * i + 3] = self.fluid.cell_color[3 * i + 1];
            self.fluid.cell_vertices[5 * i + 4] = self.fluid.cell_color[3 * i + 2];
        }

        self.program.set_draw_disk(false);
        self.program.set_point_size(self.fluid.h);

        self.vbo.set(&self.fluid.cell_vertices);
        self.ibo.set(&self.fluid.cell_indices);

        unsafe {
            gl::DrawElements(
                gl::POINTS,
                self.fluid.cell_indices.len() as i32,
                gl::UNSIGNED_INT,
                std::ptr::null(),
            );
        }
    }

    fn draw_particles(&mut self) {
        for i in 0..CONFIG.particle.total {
            self.fluid.particle_vertices[5 * i]     = self.fluid.particle_pos[2 * i];
            self.fluid.particle_vertices[5 * i + 1] = self.fluid.particle_pos[2 * i + 1];
            self.fluid.particle_vertices[5 * i + 2] = self.fluid.particle_color[3 * i];
            self.fluid.particle_vertices[5 * i + 3] = self.fluid.particle_color[3 * i + 1];
            self.fluid.particle_vertices[5 * i + 4] = self.fluid.particle_color[3 * i + 2];
        }

        self.program.set_draw_disk(true);
        self.program.set_point_size(self.fluid.particle_diameter);

        self.vbo.set(&self.fluid.particle_vertices);
        self.ibo.set(&self.fluid.particle_indices);

        unsafe {
            gl::DrawElements(
                gl::POINTS,
                CONFIG.particle.total as i32,
                gl::UNSIGNED_INT,
                std::ptr::null(),
            );
        }
    }

    // Renders the output.
    pub fn draw(&mut self) {
        self.clear_window();
        self.draw_grid();
        self.draw_particles();
    }

    // Recalculates the simulation.
    fn simulate(&mut self) {
        // Particle simulation.
        self.fluid.integrate_particles(CONFIG.scene.dt, CONFIG.environment.gravity);
        self.fluid.push_particles_apart(CONFIG.scene.particle_iterations);
        self.fluid.handle_particle_collisions(&self.obstacle);

        // El problema lo tengo aquí.
        self.fluid.transfer_velocities(true, CONFIG.environment.flip_ratio); // From particle, to grid.
        self.fluid.update_particle_density();
        self.fluid.solve_incompressibility(CONFIG.scene.pressure_iterations, CONFIG.scene.dt, CONFIG.environment.over_relaxation, CONFIG.environment.compensate_drift);
        self.fluid.transfer_velocities(false, CONFIG.environment.flip_ratio); // From grid, to particle.

        // Esto es más cosmético que otra cosa.
        self.fluid.update_particle_colors();
        self.fluid.update_cell_colors();
    }
    
    // It executes each frame.
    pub fn update(&mut self) {
        self.simulate();
        self.draw();

        self.windsl.window.gl_swap_window();
    }
}
