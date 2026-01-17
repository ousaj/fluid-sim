use crate::windsl::Windsl;
use crate::objects::*;
use crate::fluid::Fluid;
use crate::obstacle::Obstacle;
use crate::position::Position;
use crate::velocity::Velocity;

use crate::config::FRAME_RATE;
use crate::config::WINDOW_SIZE;
use crate::config::IS_GRID_VISIBLE;
use crate::config::ARE_PARTICLES_BLURRED;
use crate::config::ARE_PARTICLES_VISIBLE;
use crate::config::GLOW_MULTIPLIER;
use crate::config::BLUR_MULTIPLIER;
use crate::config::TOTAL_PARTICLES;
use crate::config::TOTAL_CELLS;
use crate::config::PARTICLE_DIAMETER;
use crate::config::CELL_SIZE;

extern crate sdl2;

use rayon::prelude::*;

pub struct Scene {
    pub windsl: Windsl,
    pub frame_rate: f32,
    pub obstacle: Obstacle,
    pub fluid: Fluid,
    vbo: Vbo,
    _vao: Vao,
    ibo: Ibo,
    program: Program,
    pub is_grid_visible: bool,
    pub are_particles_visible: bool,
    pub are_particles_blurred: bool,
    pub is_mouse_dragging: bool,
    pub particle_vertices: [f32; 5 * TOTAL_PARTICLES],
    pub particle_indices: [u32; TOTAL_PARTICLES],
    pub cell_vertices: [f32; 5 * TOTAL_CELLS],
    pub cell_indices: [u32; TOTAL_CELLS],
}

impl Scene {
    pub fn new() -> Self {
        let frame_rate: f32 = FRAME_RATE;
        let is_mouse_dragging = false;
        let obstacle: Obstacle = Obstacle::new(Position::new(0.0, 0.0), Velocity::new(0.0, 0.0));
        let windsl = Windsl::new(WINDOW_SIZE as usize, WINDOW_SIZE as usize).unwrap();
        let program = create_program().unwrap();
        let fluid: Fluid = Fluid::new();

        program.set();
        program.set_resolution(WINDOW_SIZE, WINDOW_SIZE);

        unsafe { 
            gl::Viewport(0, 0, WINDOW_SIZE as i32, WINDOW_SIZE as i32);
            gl::Enable(gl::PROGRAM_POINT_SIZE);

            // Experimental
            gl::Enable(gl::BLEND);
            gl::BlendFunc(gl::ONE, gl::ONE);
        }

        let vbo = Vbo::gen();
        vbo.set(&Vec::new());

        let _vao = Vao::gen();
        _vao.set();

        let ibo = Ibo::gen();
        ibo.set(&Vec::new());

        let is_grid_visible: bool = IS_GRID_VISIBLE;
        let are_particles_visible: bool = ARE_PARTICLES_VISIBLE;
        let are_particles_blurred: bool = ARE_PARTICLES_BLURRED;
        program.set_glow_multiplier(GLOW_MULTIPLIER);
        program.set_blur_multiplier(BLUR_MULTIPLIER);

        let particle_vertices: [f32; 5 * TOTAL_PARTICLES] = [0.0; 5 * TOTAL_PARTICLES];
        let mut particle_indices: [u32; TOTAL_PARTICLES] = [0; TOTAL_PARTICLES];
        particle_indices
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, particle_index)| {
                *particle_index = i as u32;
            });

        let cell_vertices: [f32; 5 * TOTAL_CELLS] = [0.0; 5 * TOTAL_CELLS];
        let mut cell_indices: [u32; TOTAL_CELLS] = [0; TOTAL_CELLS];
        cell_indices
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, cell_index)| {
                *cell_index = i as u32;
            });

        Scene {
            obstacle,
            fluid,
            windsl,
            vbo,
            _vao,
            ibo,
            program,
            frame_rate,
            is_grid_visible,
            are_particles_visible,
            are_particles_blurred,
            is_mouse_dragging,
            particle_vertices,
            particle_indices,
            cell_vertices,
            cell_indices,
        }
    }

    fn clear_window(&mut self) {
        unsafe {
            gl::ClearColor(0.0, 0.0, 0.0, 1.0);
            gl::Clear(gl::COLOR_BUFFER_BIT);
        }
    }

    fn draw_grid(&mut self) {
        self.program.set_draw_disk(false);
        self.program.set_point_size(CELL_SIZE);

        self.cell_vertices
            .par_chunks_mut(5)
            .zip(self.fluid.cells.par_iter())
            .for_each(|(cell_vertex, cell)| {
                cell_vertex[0] = cell.position.x;
                cell_vertex[1] = cell.position.y; 
                cell_vertex[2] = cell.color.r; 
                cell_vertex[3] = cell.color.g; 
                cell_vertex[4] = cell.color.b; 
            });

        self.vbo.set(&self.cell_vertices);
        self.ibo.set(&self.cell_indices);

        unsafe {
            gl::DrawElements(
                gl::POINTS,
                TOTAL_CELLS as i32,
                gl::UNSIGNED_INT,
                std::ptr::null(),
            );
        }
    }

    fn draw_particles(&mut self) {
        self.program.set_draw_disk(true);
        self.program.set_point_size(PARTICLE_DIAMETER);
        self.program.set_blur_particles(self.are_particles_blurred);

        self.particle_vertices
            .par_chunks_mut(5)
            .zip(self.fluid.particles.par_iter())
            .for_each(|(particle_vertex, particle)| {
                particle_vertex[0] = particle.position.x;
                particle_vertex[1] = particle.position.y;
                particle_vertex[2] = particle.color.r;
                particle_vertex[3] = particle.color.g;
                particle_vertex[4] = particle.color.b;
            });

        self.vbo.set(&self.particle_vertices);
        self.ibo.set(&self.particle_indices);

        unsafe {
            gl::DrawElements(
                gl::POINTS,
                TOTAL_PARTICLES as i32,
                gl::UNSIGNED_INT,
                std::ptr::null(),
            );
        }
    }

    // Renders the output.
    pub fn draw(&mut self) {
        self.clear_window();

        if self.is_grid_visible {
            self.draw_grid();
        }

        if self.are_particles_visible {
            self.draw_particles();
        }
    }

    // Recalculates the simulation.
    fn simulate(&mut self) {
        // // Particle simulation.
        self.fluid.integrate_particles();
        self.fluid.push_particles_apart();
        self.fluid.handle_particle_collisions(&self.obstacle);

        self.fluid.transfer_velocities(true); // From particle, to grid.
        self.fluid.update_density();
        self.fluid.solve_incompressibility();
        self.fluid.transfer_velocities(false); // From grid, to particle.

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
