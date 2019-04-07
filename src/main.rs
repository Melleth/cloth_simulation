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
        pub vertices: Vec<Point3<f32>>,
        pub indices: Vec<Point3<u16>>,
        width: usize,
        height: usize,
        mesh: Option<Mesh>,
        springs: Vec<Spring>,

        velocities: Vec<Vector3<f32>>,
        masses: Vec<f32>,
        is_fixed: Vec<bool>,
    }

    pub struct Spring {
        rest_length: f32,
        k: f32,
        k_damping: f32,
        distance_scale: f32,
        vert1: usize,
        vert2: usize,
    }

    impl Spring {
        pub fn new(rest_length: f32, force_constant: f32, k_damping: f32, distance_scale: f32, vert1: usize, vert2: usize) -> Spring {

            Spring {
                rest_length: rest_length,
                k: force_constant,
                k_damping: k_damping,
                distance_scale: distance_scale,
                vert1: vert1,
                vert2: vert2, 
            }
        }
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

            //Setup mass points
            let mass : f32 = 1.0;
            let mut masses: Vec<f32> = vec![mass; width * height];
            let mut velocities: Vec<Vector3<f32>> = vec![Vector3::new(0.0, 0.0, 0.0); width * height];
            let mut is_fixed: Vec<bool> = vec![false; width * height];

            //Fix top corners
            is_fixed[width * height - 1] = true;
            is_fixed[width * height - width] = true;

            let mut springs: Vec<Spring> = Vec::with_capacity(4 * (width * height) - 3 * (width - 1) - 3 * (height - 1) - 4);

            let rest_length : f32 = 1.0;
            let diagonal_rest_length : f32 =  (2.0 * (rest_length * rest_length)).sqrt();
            let force_constant = 1.0;
            let k_damping = 0.1;
            let distance_scale = 10.0;

            //Connect top row
            for x in 0..(width-1)  {
                springs.push(Spring::new(rest_length, force_constant, k_damping, distance_scale, x, x + 1)); //Right

                springs.push(Spring::new(rest_length, force_constant, k_damping, distance_scale, x, x + width)); //Bottom

                springs.push(Spring::new(diagonal_rest_length, force_constant, k_damping, distance_scale, x, (x + 1) + width)); //Bottom right
            }

            //Connect right column
            for y in 0..(height - 1) {
                springs.push(Spring::new(rest_length, force_constant, k_damping, distance_scale, (width-1) + (y * width), (width-1) + ((y+1) * width))); //Bottom
            }

            //Connect middle
            for y in 1..(height - 1) {
                for x in 0..(width - 1){
                    springs.push(Spring::new(diagonal_rest_length, force_constant, k_damping, distance_scale, x + (y * width), (x+1) + ((y-1) * width))); //Top right

                    springs.push(Spring::new(rest_length, force_constant, k_damping, distance_scale, x + (y * width), (x+1) + (y * width))); //Right

                    springs.push(Spring::new(diagonal_rest_length, force_constant, k_damping, distance_scale, x + (y * width), (x+1) + ((y+1) * width))); //Bottom Right

                    springs.push(Spring::new(rest_length, force_constant, k_damping, distance_scale, x + (y * width), x + ((y+1) * width))); //Bottom
                }
            }

            //Connect bottom row
            for x in 0..(width - 1) {
                springs.push(Spring::new(diagonal_rest_length, force_constant, k_damping, distance_scale, x + ((height-1) * width), (x+1) + ((height-2) * width))); //Top right
                
                springs.push(Spring::new(rest_length, force_constant, k_damping, distance_scale, x + ((height-1) * width), (x+1) + ((height-1) * width))); //Right
            }
            
            //Bottom right corner is connected through the other loops

            //TODO: Add bending springs? (skips 1 over)

            Cloth { width: width, 
                    height: height,
                    vertices: vertices,
                    indices: indices,
                    springs: springs,
                    velocities: velocities,
                    masses: masses,
                    is_fixed: is_fixed,
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

        pub fn update_velocity(&mut self, dt : f32) {
            let gravity : Vector3<f32> = Vector3::new(0.0, -9.81, 0.0);

            //Reset velocity
            for vel in self.velocities.iter_mut()  {
                *vel = gravity * dt;
            }

            //Add up spring forces
            for spring in self.springs.iter()  {
                let vert1_index : usize = spring.vert1;
                let vert2_index : usize = spring.vert2;

                //Calculate stretch term
                let distance : f32 = na::distance(&self.vertices[spring.vert1], &self.vertices[spring.vert2]);
                let stretch : f32 = spring.k * (distance - spring.rest_length) * spring.distance_scale;
                
                //Calculate damping term
                let direction_n = (self.vertices[vert2_index] - self.vertices[vert1_index]).normalize();
                let velocity_diff = self.velocities[vert2_index] - self.velocities[vert1_index];
                let damping = spring.k_damping * velocity_diff.dot(&direction_n);

                //Calculate and apply force to the connected points
                let force = (stretch + damping) * direction_n;

                self.velocities[vert1_index] += force / self.masses[vert1_index];
                self.velocities[vert2_index] -= force / self.masses[vert2_index];
            }
        }
        
        pub fn update_position(&mut self, dt : f32) {
            for i in 0..self.vertices.len() {
                //Fixed point doesnt move
                if self.is_fixed[i] {
                    continue;
                }
                self.vertices[i] += self.velocities[i] * dt;
            }
        }
        
        pub fn update(&mut self, dt : f32) {
            self.update_velocity(dt);
            self.update_position(dt);
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
    
    window.set_framerate_limit(None);
    //window.set_framerate_limit(Some(60));

    // Define the cloth, get the mesh from it and add it to the scene.
    let mut cloth = cloth::Cloth::new(30,30);
    let mut old_group = window.add_group();

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
        old_group.unlink();
        
        let mut mesh = Rc::new(RefCell::new(Mesh::new(cloth.vertices.clone(), cloth.indices.clone(), None, None, true)));

        /*MeshManager::get_global_manager(|mm| mm.add(mesh.clone(), "cloth_mesh"));
        let mut cloth_object = window
            .add_geom_with_name("cloth_mesh", Vector3::new(1.0, 1.0, 1.0))
            .unwrap();*/
        let mut new_group = window.add_group();
        let mut cloth_object = new_group.add_mesh(mesh, Vector3::new(1.0, 1.0, 1.0));
        old_group = new_group;
        cloth_object.enable_backface_culling(false);

        // Set wireframe rendering.
        cloth_object.set_color(96.4, 0.0, 0.61);
        cloth_object.set_points_size(5.0);
        cloth_object.set_lines_width(2.0);
        cloth_object.set_surface_rendering_activation(false);


        cloth.update(0.02);
    }
}