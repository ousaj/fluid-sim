use sdl2::{ video::{ GLContext, SwapInterval, Window }, EventPump, Sdl };

pub struct Windsl {
    pub _sdl: Sdl,
    pub window: Window,
    pub _gl_context: GLContext,
    pub _gl: (),
    pub event_pump: EventPump,
}

impl Windsl {
    pub fn new(width: usize, height: usize) -> Result<Self, &'static str> {
        let _sdl = sdl2::init().unwrap();
        let video_subsystem = _sdl.video().unwrap();

        let gl_attr = video_subsystem.gl_attr();
        gl_attr.set_context_profile(sdl2::video::GLProfile::Core);
        gl_attr.set_context_version(3, 3);

        let window: Window = video_subsystem
            .window("Simulación de fluidos", width as u32, height as u32)
            .opengl()
            .build()
            .unwrap();

        let _gl_context = window.gl_create_context().unwrap();
        let _gl = gl::load_with(|s| {
            video_subsystem.gl_get_proc_address(s) as *const std::os::raw::c_void
        });

        window
            .subsystem()
            .gl_set_swap_interval(SwapInterval::Immediate)
            .unwrap();

        let event_pump: sdl2::EventPump = _sdl.event_pump().unwrap();

        Ok(Windsl {
            _sdl,
            window,
            _gl_context,
            _gl,
            event_pump
        })
    }
}
