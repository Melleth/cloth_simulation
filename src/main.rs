extern crate kiss3d;
extern crate nalgebra as na;

use kiss3d::camera::{ArcBall, FirstPerson};
use kiss3d::event::{Action, Key, WindowEvent};
use kiss3d::light::Light;
use kiss3d::window::Window;
use kiss3d::resource::{Mesh, MeshManager};
use na::{Point3, UnitQuaternion, Vector3};
use std::cell::RefCell;
use std::rc::Rc;


mod cloth {
    use kiss3d::camera::{ArcBall, FirstPerson};
    use kiss3d::event::{Action, Key, WindowEvent};
    use kiss3d::light::Light;
    use kiss3d::window::Window;
    use kiss3d::resource::{Mesh, MeshManager};
    use na::{Point3, UnitQuaternion, Vector3};
    use std::cell::RefCell;
    use std::rc::Rc;
    
    pub struct Cloth {
        vertices: Vec<Point3<f32>>,
        indices: Vec<Point3<u16>>,
        widht: usize,
        height: usize,
        mesh: Option<Mesh>,
    }

    impl Cloth {
        pub fn new(width: usize, height: usize) -> Cloth {
            let mut vertices: Vec<Point3<f32>> =  vec![Point3::origin(); width*height];
            let mut indices: Vec<Point3<u16>> = vec![Point3::new(0u16, 0, 0); width * height * 2];
            for i in 0..height {
                for j in 0..width {
                    // Construct in z plane for now :)
                    vertices[j + i * width] = Point3::new(j as f32, i as f32, 0.0);
                }
            }

            // Helper var to avoid casting multiple times.
            let w = width as u16;

            // For each quad
            for i in 0..(width * height) {
                let y = (i / width) as u16;
                if y == height as u16 - 1 { break; }
                let x = (i as u16) % (w - 1);
                // Define 2 triangles.
                indices[i * 2 + 0] = Point3::new(x + y * w, x + 1 + y * w, x + (y+1) * w);
                indices[i * 2 + 1] = Point3::new(x + (y + 1) * w, x + 1 + (y + 1) * w, x + 1 + y * w);
                
                // Uncomment these for debugging purposes.
                //let indice = indices[i * 2 + 0];
                //println!("i: {}, y: {}, indice: {}, vertices: {}, {}, {}", i, y, indice, vertices[indice.coords.x as usize], vertices[indice.coords.y as usize], vertices[indice.coords.z as usize]);
            }

            // Return the cloth struct.
            Cloth { widht: width, 
                    height: height,
                    vertices: vertices,
                    indices: indices, 
                    mesh: None }
        }

        // Either returns the mesh or instantiates it and returns it.
        pub fn get_mesh(mut self) -> Mesh {
            if let Some(m) = self.mesh {
                m
            } else {
                self.mesh = Some(Mesh::new(self.vertices, self.indices, None, None, false));
                self.mesh.unwrap()
            }
        }
    }
}

use cloth::Cloth;

fn main() {
    // Camera initialization stuff.
    let eye = Point3::new(10.0f32, 10.0, 10.0);
    let at = Point3::origin();
    let mut first_person = FirstPerson::new(eye, at);
    let mut arc_ball = ArcBall::new(eye, at);
    let mut use_arc_ball = true;


    // Window initialization stuff.
    let mut window = Window::new("Cloth Simulation");
    window.set_light(Light::StickToCamera);
    window.set_background_color(0.0, 0.21, 0.53);

    // Define the cloth, get the mesh from it and add it to the scene.
    let mut cloth = cloth::Cloth::new(10,10);
    let mesh = Rc::new(RefCell::new(cloth.get_mesh()));
    MeshManager::get_global_manager(|mm| mm.add(mesh.clone(), "cloth_mesh"));
    let mut cloth_object = window
        .add_geom_with_name("cloth_mesh", Vector3::new(1.0, 1.0, 1.0))
        .unwrap();
    cloth_object.enable_backface_culling(false);

    // Set wireframe rendering.
    cloth_object.set_color(96.4, 0.0, 0.61);
    cloth_object.set_points_size(10.0);
    cloth_object.set_lines_width(2.0);
    cloth_object.set_surface_rendering_activation(false);


    // Update loop
    while !window.should_close() {
        // rotate the arc-ball camera.
        let curr_yaw = arc_ball.yaw();
        arc_ball.set_yaw(curr_yaw + 0.001);
        window.set_light(Light::StickToCamera);

        // update the current camera.
        for event in window.events().iter() {
            match event.value {
                WindowEvent::Key(key, Action::Release, _) => {
                    if key == Key::Key1 {
                        use_arc_ball = true
                    } else if key == Key::Key2 {
                        use_arc_ball = false
                    }
                }
                _ => {}
            }
        }

        if use_arc_ball {
            window.render_with_camera(&mut arc_ball);
        } else {
            window.render_with_camera(&mut first_person);
        }
    }
}