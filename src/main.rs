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
        width: usize,
        height: usize,
        mesh: Option<Mesh>,
        springs: Vec<Spring>,
        cloth_points: Vec<ClothPoint>,
    }

    #[derive(Clone)]
    pub struct ClothPoint {
        mass: f32,
        velocity: Vector3<f32>,
        spring_indices: Vec<usize>,
        is_fixed: bool,
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

    impl ClothPoint {
        pub fn new(mass: f32, velocity: Vector3<f32>, is_fixed: bool) -> ClothPoint {
            
            ClothPoint {
                mass: mass,
                velocity: velocity,
                spring_indices: Vec::new(),
                is_fixed: is_fixed,
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

            let mass : f32 = 1.0;
            
            //TODO: Maybe just tear this apart and have vectors for each value, might be faster for resseting the velocity every t?
            let mut cloth_points: Vec<ClothPoint> = vec![ClothPoint::new(mass, Vector3::new(0.0, 0.0, 0.0), false); width * height];

            //Fix top corners
            cloth_points[0].is_fixed = true;
            cloth_points[width-1].is_fixed = true;

            let mut springs: Vec<Spring> = Vec::with_capacity(4 * (width * height) - 3 * (width - 1) - 3 * (height - 1) - 4);

            let rest_length : f32 = 1.0;
            let diagonal_rest_length : f32 =  (2.0 * (rest_length * rest_length)).sqrt();
            let force_constant = 1.0;
            let k_damping = 0.1;
            let distance_scale = 10.0;

            //Connect top row
            for x in 0..width  {
                cloth_points[x].spring_indices.push(springs.len());
                springs.push(Spring::new(rest_length, force_constant, k_damping, distance_scale, x, x + 1)); //Right

                cloth_points[x].spring_indices.push(springs.len());
                springs.push(Spring::new(rest_length, force_constant, k_damping, distance_scale, x, x + width)); //Bottom

                cloth_points[x].spring_indices.push(springs.len());
                springs.push(Spring::new(diagonal_rest_length, force_constant, k_damping, distance_scale, x, (x + 1) + width)); //Bottom right
            }

            //Connect right column
            for y in 0..(height - 1) {
                cloth_points[(width-1) + (y * width)].spring_indices.push(springs.len());
                springs.push(Spring::new(rest_length, force_constant, k_damping, distance_scale, (width-1) + (y * width), (width-1) + ((y+1) * width))); //Bottom
            }

            //Connect middle
            for y in 1..(height - 1) {
                for x in 0..(width - 1){
                    cloth_points[x + (y * width)].spring_indices.push(springs.len());
                    springs.push(Spring::new(diagonal_rest_length, force_constant, k_damping, distance_scale, x + (y * width), (x+1) + ((y-1) * width))); //Top right

                    cloth_points[x + (y * width)].spring_indices.push(springs.len());
                    springs.push(Spring::new(rest_length, force_constant, k_damping, distance_scale, x + (y * width), (x+1) + (y * width))); //Right

                    cloth_points[x + (y * width)].spring_indices.push(springs.len());
                    springs.push(Spring::new(diagonal_rest_length, force_constant, k_damping, distance_scale, x + (y * width), (x+1) + ((y+1) * width))); //Bottom Right

                    cloth_points[x + (y * width)].spring_indices.push(springs.len());
                    springs.push(Spring::new(rest_length, force_constant, k_damping, distance_scale, x + (y * width), x + ((y+1) * width))); //Bottom
                }
            }

            //Connect bottom row
            for x in 0..(width - 1) {
                cloth_points[x + ((height-1) * width)].spring_indices.push(springs.len());
                springs.push(Spring::new(diagonal_rest_length, force_constant, k_damping, distance_scale, x + ((height-1) * width), (x+1) + ((height-2) * width))); //Top right
                
                cloth_points[x + ((height-1) * width)].spring_indices.push(springs.len());
                springs.push(Spring::new(rest_length, force_constant, k_damping, distance_scale, x + ((height-1) * width), (x+1) + ((height-1) * width))); //Right
            }
            
            //Bottom right corner is connected through the other loops

            //TODO: Add bending springs? (skips 1 over)

            Cloth { width: width, 
                    height: height,
                    vertices: vertices,
                    indices: indices,
                    springs: springs,
                    cloth_points: cloth_points,
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
            for cloth_point in self.cloth_points.iter_mut()  {
                cloth_point.velocity = gravity * dt;
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
                let velocity_diff = self.cloth_points[vert2_index].velocity - self.cloth_points[vert1_index].velocity;
                let damping = spring.k_damping * velocity_diff.dot(&direction_n);

                //Calculate and apply force to the connected points
                let force = (stretch + damping) * direction_n;

                self.cloth_points[vert1_index].velocity += force / self.cloth_points[vert1_index].mass;
                self.cloth_points[vert2_index].velocity -= force / self.cloth_points[vert2_index].mass;
            }
        }
        
        pub fn update_position(&mut self, dt : f32) {
            for i in 0..self.vertices.len() {
                //Fixed point doesnt move
                if self.cloth_points[i].is_fixed {
                    continue;
                }
                self.vertices[i] += self.cloth_points[i].velocity * dt;
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