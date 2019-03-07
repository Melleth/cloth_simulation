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

fn main() {
    let eye = Point3::new(10.0f32, 10.0, 10.0);
    let at = Point3::origin();
    let mut first_person = FirstPerson::new(eye, at);
    let mut arc_ball = ArcBall::new(eye, at);
    let mut use_arc_ball = true;

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