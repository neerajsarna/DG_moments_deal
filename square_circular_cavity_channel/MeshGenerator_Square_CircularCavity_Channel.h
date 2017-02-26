	// generate a ring shaped mesh internally
	// template<int dim>
	// void
	// MeshGenerator::Base_MeshGenerator<dim>::mesh_internal_square_circular_cavity_channel(const unsigned int part_x,const unsigned int part_y)
	// {
 //            triangulation.clear();
 //            Point<dim> p1;
 //            Point<dim> p2;
 //            const double inner_radius = 0.25;
 //            double outer_radius;
 //            const double length_in_z = 0.0;
 //            const int repetitions_in_z = 0;
 //            std::vector<unsigned int > repetitions(dim);

 //            p1(0) = constants.xl;
 //            p1(1) = constants.yb;

 //            p2(0) = constants.xr;
 //            p2(1) = constants.yt;

 //            repetitions[0] = part_x;
 //            repetitions[1] = part_y;

 //            // the outer radius will be half of the edge length
 //            outer_radius = (constants.xr - constants.xl)/2;


 //            //The diagonal of the rectangle is the line joining p1 and p2
 //            GridGenerator::hyper_cube_with_cylindrical_hole (triangulation,
 //                                                            inner_radius,
 //                                                             outer_radius,
 //                                                            length_in_z, 
 //                                                            repetitions_in_z);

 //            //set the boundary id
 //            // prescribe the ids for the channel
 //            set_square_circular_cavity_channel_bid();
 //            triangulation.set_manifold(156,boundary);
 //            triangulation.refine_global(2);
	// }

    template<int dim>
    void
    MeshGenerator::Base_MeshGenerator<dim>::mesh_internal_square_circular_cavity_channel(const unsigned int part_x,const unsigned int part_y)
    {
            triangulation.clear();
            gridin.attach_triangulation(triangulation);
            std::ifstream f("../../square_circular_cavity_channel/square_circular_cavity.msh");
            gridin.read_msh(f);
    }


    template<int dim>
    void
    MeshGenerator::Base_MeshGenerator<dim>::set_square_circular_cavity_channel_bid() const
    {
        typename Triangulation<dim>::cell_iterator cell = triangulation.begin(),
                                                    endc = triangulation.end();

        for (; cell != endc ; cell++)
        {
                for (unsigned int face = 0 ; face < GeometryInfo<dim>::faces_per_cell ; face++)
                  if (cell->face(face)->at_boundary())
                  { 
                    bool square_edge = false;
                    double x_cord = cell->face(face)->center()(0);
                    double y_cord = cell->face(face)->center()(1);

                    // Following is the boundary ids:
                    // Left Wall = 0
                    // Bottom Wall = 1
                    // Right Wall = 2
                    // Top Wall = 3
                    // the periodic faces get the id 100 and 101
                    // right edge
                    
                    if (x_cord == constants.xr)
                    {
                      cell->face(face)->set_boundary_id(102);
                      square_edge = true;
                    }

                    // left edge, the inflow boundary
                    if (x_cord == constants.xl)
                    {
                      square_edge = true;  
                      cell->face(face)->set_boundary_id(101);
                    }

                    // this is the bottom wall
                    if (y_cord == constants.yb)
                    {
                        square_edge = true;
                      cell->face(face)->set_boundary_id(0);
                    }

                    // top edge, This is the second wall
                     if (y_cord == constants.yt)
                     {
                         square_edge = true;
                        cell->face(face)->set_boundary_id(1);
                     }

                     // if the square_edge has not been found then the face is at the circular region
                     if(square_edge==false)
                        {
                            cell->face(face)->set_boundary_id(2);    
                            cell->face(face)->set_manifold_id(156);
                        }
                        
                   }
        }

  }

        
