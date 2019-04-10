extern crate kiss3d;
extern crate nalgebra as na;
extern crate rayon;

mod cloth;

use kiss3d::camera::{ArcBall, FirstPerson};
use kiss3d::event::{Action, Key, WindowEvent};
use kiss3d::light::Light;
use kiss3d::window::Window;
use kiss3d::resource::{Mesh, MeshManager};
use na::{Translation3, Point3, UnitQuaternion, Vector3};
use std::cell::RefCell;
use std::rc::Rc;
use std::time::{Duration, Instant};

use cloth::Cloth;

fn main() {

    let width = 33;
    let height = 33;

    // Camera initialization stuff.
    let eye = Point3::new(-12.0f32, 20.0, -12.0);
    //let at = Point3::origin();
    let at = Point3::new(width as f32,-20.0,height as f32);
    let mut first_person = FirstPerson::new(eye, at);
    let mut arc_ball = ArcBall::new(eye, at);
    let mut use_arc_ball = false;


    // Window initialization stuff.
    let mut window = Window::new("Cloth Simulation");
    window.set_light(Light::StickToCamera);
    window.set_background_color(0.0, 0.21, 0.53);

    // Define the cloth, get the mesh from it and add it to the scene.
    let mut cloth = Cloth::new(width,height);
    //cloth.fix_edges();
    let mut old_group = window.add_group();
    let mut simulate = false;

    //TODO: Add sphere pos/size settings nicely 
    //Draw a nice colored sphere
    let sphere_size : f32 = 5.0;
    let mut s = window.add_sphere(sphere_size);
    s.set_color(0.9506297, 0.9983103, 0.95816237);
    s.append_translation(&Translation3::new((cloth.width as f32)/2.0, -(sphere_size * 1.5), (cloth.width as f32)/3.0));
    s.set_texture_from_file(&std::path::Path::new("./ball.png"), "ball");
    
    
    // Update loop
    let mut num_iter = 1;
    let mut show_wireframe = false;
    while !window.should_close() {
        // rotate the arc-ball camera.
        let curr_yaw = arc_ball.yaw();
        arc_ball.set_yaw(curr_yaw + 0.001);
        window.set_light(Light::StickToCamera);

        // update the current camera.
        for event in window.events().iter() {
            match event.value {
                WindowEvent::Key(key, Action::Release, _) => {
                    match key {
                        Key::Key1 => use_arc_ball = true,
                        Key::Key2 => use_arc_ball = false,
                        Key::Space => simulate = !simulate,
                        Key::W => show_wireframe = !show_wireframe,
                        _ => ()
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
        
        let mut mesh = Rc::new(RefCell::new(Mesh::new(cloth.vertices.clone(), cloth.indices.clone(), None, Some(cloth.uv.clone()), true)));

        let mut new_group = window.add_group();
        let mut cloth_object = new_group.add_mesh(mesh, Vector3::new(1.0, 1.0, 1.0));
        old_group = new_group;
        cloth_object.enable_backface_culling(false);
        //cloth_object.set_color(0.964, 0.0, 0.61);
        cloth_object.set_texture_from_file(&std::path::Path::new("./sPikachuHD.png"), "memes");

        
        if show_wireframe {
            cloth_object.set_points_size(10.0);
            cloth_object.set_lines_width(2.0);
            cloth_object.set_surface_rendering_activation(false);
        }

        let time = Instant::now();
        if simulate {
            for i in 0..num_iter {
                cloth.update(0.02);
            }

            let duration = time.elapsed();
            if duration < Duration::from_millis(16) {
                num_iter += 1;
            } else {
                num_iter -= 1;
            }
            println!("Did {} iterations at {:?}", num_iter, duration);
        }
    }
}