use na::{Point2, Point3, Vector3};
use rayon::prelude::*;

mod grid;
use grid::Grid;

const USE_MULTITHREADING: bool = true;


pub struct Cloth {
    pub vertices: Vec<Point3<f32>>,
    pub indices: Vec<Point3<u16>>,
    pub uv: Vec<Point2<f32>>,
    pub width: usize,
    pub height: usize,
    springs: Vec<Spring>,
    velocities: Vec<Vector3<f32>>,
    masses: Vec<f32>,
    is_fixed: Vec<bool>,
    grid : Grid,
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
        let mut UV: Vec<Point2<f32>> = vec![Point2::origin(); width*height];
        let mut indices: Vec<Point3<u16>> = vec![Point3::new(0u16, 0, 0); width * height * 2];
        //Grid
        let grid_size = std::cmp::max(width, height) * 5;
        let mut grid = Grid::new(grid_size, 1.0, Vector3::new((width / 2) as f32, (height / 2) as f32, 0.0));

        for i in 0..height {
            for j in 0..width {
                // Construct in z plane for now :)
                vertices[j + i * width] = Point3::new(j as f32, 0.0, i as f32);
                UV[j+i*width] = Point2::new(1.0 - (j as f32) / (width as f32), 1.0 - (i as f32) / (height as f32));
                //Add verts to grid
                grid.add_vert(j + i * width, vertices[j + i * width]);
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
        let mass : f32 = 0.5;
        let masses: Vec<f32> = vec![mass; width * height];
        let velocities: Vec<Vector3<f32>> = vec![Vector3::new(0.0, 0.0, 0.0); width * height];
        let mut is_fixed: Vec<bool> = vec![false; width * height];

        //Fix top corners
        is_fixed[width * height - 1] = true;
        is_fixed[width * height - width] = true;

        let mut springs: Vec<Spring> = Vec::with_capacity(4 * (width * height) - 3 * (width - 1) - 3 * (height - 1) - 4);

        let rest_length : f32 = 1.0;
        let diagonal_rest_length : f32 =  (2.0 * (rest_length * rest_length)).sqrt();
        let force_constant = 10.0; //k
        let k_damping = 0.1;
        let distance_scale = 10.0;

        //Connect top row
        for x in 0..(width - 1)  {
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

        Cloth { width: width, 
                height: height,
                vertices: vertices,
                uv: UV,
                indices: indices,
                springs: springs,
                velocities: velocities,
                masses: masses,
                is_fixed: is_fixed,
                grid: grid }
    }

    //Calculates the normal for a triangle defined by the given points
    fn triangle_normal(&self, p1 : Point3<f32>, p2 :  Point3<f32>, p3 :  Point3<f32>) -> Vector3<f32>
    {
        let v1 = p2-p1;
        let v2 = p3-p1;

        return v1.cross(&v2);
    }

    //Get the force the wind enacts upon a triangles surface
    fn get_wind_force(&self, p1 : Point3<f32>, p2 :  Point3<f32>, p3 :  Point3<f32>, wind_vector : Vector3<f32>) -> Vector3<f32>
    {
        let normal = self.triangle_normal(p1,p2,p3);
        let normal_n = normal.normalize();
        let wind_force = normal * (normal_n.dot(&wind_vector));
        
        return wind_force;
    }
    //Adds up the gravity and spring forces, and sets the resulting velocities.
    pub fn update_velocity(&mut self, dt : f32) {
        let gravity : Vector3<f32> = Vector3::new(0.0, -9.81, 0.0);

        //Reset velocity
        for v in 0..self.velocities.len()  {
            self.velocities[v] = gravity * self.masses[v] * dt;
        }

        if USE_MULTITHREADING {
            //Add up spring forces
            let mut parallel_results: Vec<(usize, usize, Vector3<f32>)>;
            parallel_results = (0..self.springs.len()).into_par_iter().map(|i| {
                let spring = &self.springs[i];
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
                (vert1_index, vert2_index, force * dt / self.masses[vert1_index])
            }).collect();
            
            for (v1, v2, result) in parallel_results {
                self.velocities[v1] += result;
                self.velocities[v2] -= result;
            }
        } else {
            for spring in &self.springs {
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
                self.velocities[vert1_index] += force * dt / self.masses[vert1_index];
                self.velocities[vert2_index] -= force * dt / self.masses[vert2_index];
            }
        }

        let wind_vector : Vector3<f32> = Vector3::new(0.5, 0.0, 0.5);

        //Add up wind forces
        for x in 0..(self.width - 1)  {
            for y in 0..(self.height - 1){
                
                //Get wind forces
                let wind_force1 = self.get_wind_force(self.vertices[x + y * self.width], self.vertices[x + 1 + y * self.width], self.vertices[x + (y + 1) * self.width], wind_vector);
                let wind_force2 = self.get_wind_force(self.vertices[x + (y + 1) * self.width], self.vertices[x + 1 + (y + 1) * self.width], self.vertices[x + 1 + y * self.width], wind_vector);

                //Triangle 1
                self.velocities[x + y * self.width] += wind_force1 * dt / self.masses[x + y * self.width];
                self.velocities[x + 1 + y * self.width] += wind_force1 * dt / self.masses[x + 1 + y * self.width];
                self.velocities[x + (y+1) * self.width] += wind_force1 * dt / self.masses[x + (y+1) * self.width];
                //Triangle2
                self.velocities[x + (y + 1) * self.width] += wind_force2 * dt / self.masses[x + (y + 1) * self.width];
                self.velocities[x + 1 + (y + 1) * self.width] += wind_force2 * dt / self.masses[x + 1 + (y + 1) * self.width];
                self.velocities[x + 1 + y * self.width] += wind_force2 * dt / self.masses[x + 1 + y * self.width];
            }
        }
    }
    
    //Updates the positions according to the velocities.
    pub fn update_position(&mut self, dt: f32) {
        for i in 0..self.vertices.len() {
            //Fixed point doesnt move
            if self.is_fixed[i] {
                continue;
            }
            let old_pos = self.vertices[i];
            self.vertices[i] += self.velocities[i] * dt;

            self.grid.move_vert(i, old_pos, self.vertices[i]);
        }
    }
    
    pub fn update(&mut self, dt : f32) {
        //TODO: Add sphere pos/size settings nicely
        self.sphere_collision(Point3::new((self.width as f32)/2.0, -(5.0 * 1.5), (self.width as f32)/3.0), 5.0, 0.2);
        self.self_intersection(1.0);
        self.update_velocity(dt);
        self.update_position(dt);
    }

    //Directional rejection sphere - point collision detection and correction.
    pub fn sphere_collision(&mut self, sphere_position: Point3<f32>, sphere_radius: f32, delta: f32) {
        let sphere_radius_delta = sphere_radius + delta;

        for i in 0..self.vertices.len()  {
            //Calculate distance from sphere
            let distance_vec = self.vertices[i] - sphere_position;
            let distance_norm = distance_vec.norm();

            //If inside sphere, correct
            if distance_norm < sphere_radius_delta {

                let mut pen_normal = distance_vec.normalize();
                let d_dot_n = distance_vec.dot(&pen_normal);

                //Directional rejection
                if d_dot_n < 0.0 {
                    pen_normal = distance_vec - 2.0 * d_dot_n * pen_normal;
                }

                let old_pos = self.vertices[i];
                self.vertices[i] += pen_normal * (sphere_radius_delta - distance_norm);

                self.grid.move_vert(i, old_pos, self.vertices[i]);
            }
        }
    }

    //self intersection resolution by putting a collision sphere around all the verts, don't correct for neighbours.
    pub fn self_intersection(&mut self, sphere_radius : f32) {

        let double_radius = 2.0 * sphere_radius;

        for y in 0..(self.height) {
            for x in 0..(self.width) {
                let vert_index = x + y * self.width;
                        
                        //Retrieve the potential collision candidates from the grid
                        let near_indices = self.grid.get_near_indices(self.vertices[vert_index], double_radius);
                        
                        //Skip over the neighbours of this point
                        let mut neighbours : Vec<usize> = Vec::new();
                        if x != 0                                       {neighbours.push((x-1) + ( y    * self.width));} //left
                        if x != 0 && y != 0                             {neighbours.push((x-1) + ((y-1) * self.width));} //top-left
                        if x != 0 && y != (self.height-1)               {neighbours.push((x-1) + ((y+1) * self.width));} //bottom-left
                        if x != (self.width-1) && y != 0                {neighbours.push((x+1) + ((y-1) * self.width));} //top-right
                        if x != (self.width-1) && y != (self.height-1)  {neighbours.push((x+1) + ((y+1) * self.width));} //bottom-right
                        if x != (self.width-1)                          {neighbours.push((x+1) +  (y    * self.width));} //right
                        if y != 0                                       {neighbours.push( x    + ((y-1) * self.width));} //top
                        if y != (self.height - 1)                       {neighbours.push( x    + ((y+1) * self.width));} //bottom

                        for &collision_index in near_indices.iter() {
                                //skip self or neighbour
                            if !(vert_index == collision_index || neighbours.contains(&collision_index)) {
                                    //Calculate distance between spheres
                                    let distance_vec = self.vertices[collision_index] - self.vertices[vert_index];
                                    let distance_norm = distance_vec.norm();

                                    //If inside sphere, correct
                                    if distance_norm < (double_radius) {

                                        let mut pen_normal = distance_vec / distance_norm;
                                        let d_dot_n = distance_vec.dot(&pen_normal);

                                        //Directional rejection
                                        if d_dot_n < 0.0 {
                                            pen_normal = distance_vec - 2.0 * d_dot_n * pen_normal;
                                        }
                                        let pen_depth = (2.0 * sphere_radius) - distance_norm;

                                        let old_pos = self.vertices[vert_index];
                                        let old_col_pos = self.vertices[collision_index];

                                        self.vertices[vert_index] -= pen_normal * (pen_depth / 2.0);
                                        self.vertices[collision_index] += pen_normal * (pen_depth / 2.0);

                                        //Move in grid
                                        self.grid.move_vert(vert_index, old_pos, self.vertices[vert_index]);
                                        self.grid.move_vert(collision_index, old_col_pos, self.vertices[collision_index]);
                                }
                            }
                        
                    }
                
            }
        }
    }
}
