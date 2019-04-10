use na::{Vector3, Point3};

pub struct Grid {
    cells : Vec<Vec<Vec<Vec<usize>>>>,
    cell_dimensions : f32,
    grid_position : Vector3<f32>, //min corner
    cell_count : usize,
}

impl Grid {
    pub fn new(grid_dimension : usize, cell_dimension : f32, grid_center : Vector3<f32>) -> Grid {
        
        let cell_count = (grid_dimension as f32/cell_dimension + 1.0).ceil()as usize;
        let cells = vec![vec![vec![Vec::new();cell_count];cell_count];cell_count];
        let half_dim = (grid_dimension as f32) / 2.0;
        let half_dim_vec : Vector3<f32> = Vector3::new(half_dim, half_dim, half_dim);

        Grid {
            cells : cells,
            cell_dimensions : cell_dimension,
            grid_position : grid_center - half_dim_vec,
            cell_count : cell_count,
        }
    }

    pub fn add_vert(&mut self, vert_index : usize, vert_pos : Point3<f32>) {
        
        //Move to grid space.
        let grid_vert_pos = vert_pos - self.grid_position;

        //Calculate grid cell and push
        let x = (grid_vert_pos.coords[0] / self.cell_dimensions) as usize;
        let y = (grid_vert_pos.coords[1] / self.cell_dimensions) as usize;
        let z = (grid_vert_pos.coords[2] / self.cell_dimensions) as usize;
        self.cells[x][y][z].push(vert_index);
    }

    pub fn remove_vert(&mut self, vert_index : usize, vert_old_pos : Point3<f32>) {

        //Move origin to center of grid.
        let grid_vert_old_pos = vert_old_pos - self.grid_position;

        //Calculate grid cell and push
        let x_o = (grid_vert_old_pos.coords[0] / self.cell_dimensions) as usize;
        let y_o = (grid_vert_old_pos.coords[1] / self.cell_dimensions) as usize;
        let z_o = (grid_vert_old_pos.coords[2] / self.cell_dimensions) as usize;
        
        //Remove from old
        self.cells[x_o][y_o][z_o].retain(|&e| e != vert_index)
    }

    pub fn move_vert(&mut self, vert_index : usize, vert_old_pos : Point3<f32>, vert_new_pos : Point3<f32>) {
        
        //Move origin to center of grid.
        let grid_vert_old_pos = vert_old_pos - self.grid_position;
        let grid_vert_new_pos = vert_new_pos - self.grid_position;

        //Calculate grid cell and push
        let x_o = (grid_vert_old_pos.coords[0] / self.cell_dimensions) as usize;
        let y_o = (grid_vert_old_pos.coords[1] / self.cell_dimensions) as usize;
        let z_o = (grid_vert_old_pos.coords[2] / self.cell_dimensions) as usize;

        //Calculate grid cell and push
        let x = (grid_vert_new_pos.coords[0] / self.cell_dimensions) as usize;
        let y = (grid_vert_new_pos.coords[1] / self.cell_dimensions) as usize;
        let z = (grid_vert_new_pos.coords[2] / self.cell_dimensions) as usize;
        
        if (x != x_o) || (y != y_o) || (z != z_o) {
            //Remove from old and add to new
            self.cells[x_o][y_o][z_o].retain(|&e| e != vert_index);
            self.cells[x][y][z].push(vert_index);
        }
    }

    pub fn get_near_indices(&self, position : Point3<f32>, radius : f32) -> Vec<usize>
    {
        //Calculate collision box by contructing min and max corners from position in grid space and radius
        let corner_unit : Vector3<f32> = Vector3::new(1.0, 1.0, 1.0);
        let corner_vec = corner_unit * radius; 

        let min = ((position - self.grid_position) - corner_vec) / self.cell_dimensions;
        let max = ((position - self.grid_position) + corner_vec) / self.cell_dimensions;
        
        //Min corner
        let x_min = min.coords[0] as usize;
        let y_min = min.coords[1] as usize;
        let z_min = min.coords[2] as usize;

        //Max corner
        let x_max = max.coords[0].ceil() as usize;
        let y_max = max.coords[1].ceil() as usize;
        let z_max = max.coords[2].ceil() as usize;

        //Gather all the points in the overlapping cells, and return
        let mut near_indices : Vec<usize> = Vec::new();

        for x in x_min..x_max {
            for y in y_min..y_max {
                for z in z_min..z_max {
                    for n in 0..self.cells[x][y][z].len() {
                        near_indices.push(self.cells[x][y][z][n]);
                    }
                }
            }
        }

        return near_indices;
    }
}