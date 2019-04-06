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
        widht: usize,
        height: usize,
        springs: Vec<Spring>,
        clothPoints: Vec<ClothPoint>,
    }

    #[derive(Clone)]
    pub struct ClothPoint {
        mass: f32,
        velocity: Vector3<f32>,
        springIndices: Vec<usize>,
        isFixed: bool,
    }

    pub struct Spring {
        
        restLength: f32,
        k: f32,
        kDamping: f32,
        distanceScale: f32,

        vert1: usize,
        vert2: usize,
    }

    impl Spring {
        pub fn new(restLength: f32, forceConstant: f32, kDamping: f32, distanceScale: f32, vert1: usize, vert2: usize) -> Spring {

            Spring {
                restLength: restLength,
                k: forceConstant,
                kDamping: kDamping,
                distanceScale: distanceScale,
                vert1: vert1,
                vert2: vert2, 
            }
        }
    }

    impl ClothPoint {
        pub fn new( mass: f32, velocity: Vector3<f32>, isFixed: bool) -> ClothPoint {
            
            ClothPoint {
                mass: mass,
                velocity: velocity,
                springIndices: Vec::new(),
                isFixed: isFixed,
            }
        }
    }

    impl Cloth {
        pub fn new(width: usize, height: usize) -> Cloth {
            let mut vertices: Vec<Point3<f32>> =  vec![Point3::origin(); width*height];
            for i in 0..height {
                for j in 0..width {
                    // construct in z plane
                    vertices[j + i * width] = Point3::new(j as f32, i as f32, 0.0);
                }
            }

            let mass : f32 = 1.0;
            
            //TODO: Maybe just tear this apart and have vectors for each value, might be faster for resseting the velocity every t?
            let mut clothPoints: Vec<ClothPoint> = vec![ClothPoint::new(mass, Vector3::new(0.0, 0.0, 0.0), false); width * height];

            //Fix top corners
            clothPoints[0].isFixed = true;
            clothPoints[width-1].isFixed = true;

            let mut springs: Vec<Spring> = Vec::with_capacity(4 * (width * height) - 3 * (width - 1) - 3 * (height - 1) - 4);

            let restLength : f32 = 1.0;
            let diagonalRestLength : f32 =  (2.0 * (restLength * restLength)).sqrt();
            let forceConstant = 1.0;
            let kDamping = 0.1;
            let distanceScale = 10.0;

            //Connect top row
            for x in 0..width  {
                clothPoints[x].springIndices.push(springs.len());
                springs.push(Spring::new(restLength, forceConstant, kDamping, distanceScale, x, x + 1)); //Right

                clothPoints[x].springIndices.push(springs.len());
                springs.push(Spring::new(restLength, forceConstant, kDamping, distanceScale, x, x + width)); //Bottom

                clothPoints[x].springIndices.push(springs.len());
                springs.push(Spring::new(diagonalRestLength, forceConstant, kDamping, distanceScale, x, (x + 1) + width)); //Bottom right
            }

            //Connect right column
            for y in 0..(height - 1) {
                clothPoints[(width-1) + (y * width)].springIndices.push(springs.len());
                springs.push(Spring::new(restLength, forceConstant, kDamping, distanceScale, (width-1) + (y * width), (width-1) + ((y+1) * width))); //Bottom
            }

            //Connect middle
            for y in 1..(height - 1) {
                for x in 0..(width - 1){
                    clothPoints[x + (y * width)].springIndices.push(springs.len());
                    springs.push(Spring::new(diagonalRestLength, forceConstant, kDamping, distanceScale, x + (y * width), (x+1) + ((y-1) * width))); //Top right

                    clothPoints[x + (y * width)].springIndices.push(springs.len());
                    springs.push(Spring::new(restLength, forceConstant, kDamping, distanceScale, x + (y * width), (x+1) + (y * width))); //Right

                    clothPoints[x + (y * width)].springIndices.push(springs.len());
                    springs.push(Spring::new(diagonalRestLength, forceConstant, kDamping, distanceScale, x + (y * width), (x+1) + ((y+1) * width))); //Bottom Right

                    clothPoints[x + (y * width)].springIndices.push(springs.len());
                    springs.push(Spring::new(restLength, forceConstant, kDamping, distanceScale, x + (y * width), x + ((y+1) * width))); //Bottom
                }
            }

            //Connect bottom row
            for x in 0..(width - 1) {
                clothPoints[x + ((height-1) * width)].springIndices.push(springs.len());
                springs.push(Spring::new(diagonalRestLength, forceConstant, kDamping, distanceScale, x + ((height-1) * width), (x+1) + ((height-2) * width))); //Top right
                
                clothPoints[x + ((height-1) * width)].springIndices.push(springs.len());
                springs.push(Spring::new(restLength, forceConstant, kDamping, distanceScale, x + ((height-1) * width), (x+1) + ((height-1) * width))); //Right
            }
            
            //Bottom right corner is connected through the other loops

            //TODO: Add bending springs? (skips 1 over)

            Cloth { widht: width, 
                    height: height,
                    vertices: vertices,
                    springs: springs,
                    clothPoints: clothPoints }
        }

        pub fn UpdateVelocity(&mut self, dt : f32) {
            
            let gravity : Vector3<f32> = Vector3::new(0.0, -9.81, 0.0);

            //Reset velocity
            for clothPoint in self.clothPoints.iter_mut()  {
                clothPoint.velocity = gravity * dt;
            }

            //Add up spring forces
            for spring in self.springs.iter()  {
                let vert1Index : usize = spring.vert1;
                let vert2Index : usize = spring.vert2;

                //Calculate stretch term
                let distance : f32 = na::distance(&self.vertices[spring.vert1], &self.vertices[spring.vert2]);
                let stretch : f32 = spring.k * (distance - spring.restLength) * spring.distanceScale;
                
                //Calculate damping term
                let directionN = (self.vertices[vert2Index] - self.vertices[vert1Index]).normalize();
                let velocityDiff = self.clothPoints[vert2Index].velocity - self.clothPoints[vert1Index].velocity;
                let damping = spring.kDamping * velocityDiff.dot(&directionN);

                //Calculate and apply force to the connected points
                let force = (stretch + damping) * directionN;

                self.clothPoints[vert1Index].velocity += force / self.clothPoints[vert1Index].mass;
                self.clothPoints[vert2Index].velocity -= force / self.clothPoints[vert2Index].mass;
            }
        }
        
        pub fn UpdatePosition(&mut self, dt : f32) {
            
            for i in 0..self.vertices.len() {
                
                //Fixed point doesnt move
                if self.clothPoints[i].isFixed {
                    continue;
                }

                self.vertices[i] += self.clothPoints[i].velocity * dt;
            }
        }
        
        pub fn Update(&mut self, dt : f32) {

            self.UpdateVelocity(dt);
            self.UpdatePosition(dt);
        }
    }
}

use cloth::Cloth;




fn main() {
    let eye = Point3::new(10.0f32, 10.0, 10.0);
    let at = Point3::origin();
    let mut first_person = FirstPerson::new(eye, at);
    let mut arc_ball = ArcBall::new(eye, at);
    let mut use_arc_ball = true;

    let mut cloth = cloth::Cloth::new(10,10);

    let mut window = Window::new("Cloth Simulation");
    window.set_light(Light::StickToCamera);

    while !window.should_close() {
        // rotate the arc-ball camera.
        let curr_yaw = arc_ball.yaw();
        arc_ball.set_yaw(curr_yaw + 0.05);

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

        // Draw origin axes
        window.draw_line(
            &Point3::origin(),
            &Point3::new(1.0, 0.0, 0.0),
            &Point3::new(1.0, 0.0, 0.0),
        );
        window.draw_line(
            &Point3::origin(),
            &Point3::new(0.0, 1.0, 0.0),
            &Point3::new(0.0, 1.0, 0.0),
        );
        window.draw_line(
            &Point3::origin(),
            &Point3::new(0.0, 0.0, 1.0),
            &Point3::new(0.0, 0.0, 1.0),
        );

        // Draw shaded mesh.
        let a = Point3::new(-1.0, -1.0, 0.0);
        let b = Point3::new(1.0, -1.0, 0.0);
        let c = Point3::new(0.0, 1.0, 0.0);

        let vertices = vec![a, b, c];
        let indices = vec![Point3::new(0u16, 1, 2)];

        let mesh = Rc::new(RefCell::new(Mesh::new(
            vertices, indices, None, None, false,
        )));

        // XXX: it would be better to do: MeshManager::add(Rc....) directly.
        MeshManager::get_global_manager(|mm| mm.add(mesh.clone(), "custom_mesh"));

        let mut c1 = window
            .add_geom_with_name("custom_mesh", Vector3::new(1.0, 1.0, 1.0))
            .unwrap();
        let mut c2 = window
            .add_geom_with_name("custom_mesh", Vector3::new(1.0, 1.0, 1.0))
            .unwrap();
        if use_arc_ball {
            window.render_with_camera(&mut arc_ball);
        } else {
            window.render_with_camera(&mut first_person);
        }
    }
}