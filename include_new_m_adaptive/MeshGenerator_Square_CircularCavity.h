	// generate a ring shaped mesh internally
	template<int dim>
	void
	MeshGenerator::Base_MeshGenerator<dim>::mesh_internal_square_circular_cavity()
	{
            triangulation.clear();
            const double inner_radius = 0.5;
            double outer_radius = 2.0;
            const double length_in_z = 0.0;
            const int repetitions_in_z = 0;
            std::vector<unsigned int > repetitions(dim);



            repetitions[0] = part_x;
            repetitions[1] = part_y;

            //The diagonal of the rectangle is the line joining p1 and p2
            GridGenerator::hyper_cube_with_cylindrical_hole (triangulation,
                                                            inner_radius,
                                                             outer_radius,
                                                            length_in_z, 
                                                            repetitions_in_z);

            switch(problem_type)
                {
                    case heat_conduction:
                    {
                        set_square_circular_cavity_bid();
                        triangulation.set_manifold(100,boundary);
                        triangulation.refine_global(initial_refinement);
                        break;
                    }
                    case inflow_outflow:
                    {
                        break;
                    }

                    case lid_driven_cavity:
                    {
                        AssertThrow(1 == 0 ,ExcNotImplemented());
                        break;
                    }

                    case periodic:
                    {
                        AssertThrow(1 == 0,ExcNotImplemented());
                        break;
                    }

                    default:
                    {
                        AssertThrow(1 == 0,ExcMessage("Should not have reached here"));
                        break;
                    }
                }

                
            
	}
    
    template<int dim>
    void
    MeshGenerator::Base_MeshGenerator<dim>::set_square_circular_cavity_bid() const
    {
        typename Triangulation<dim>::cell_iterator cell = triangulation.begin(),
                                                    endc = triangulation.end();

        const double left_edge = -0.5;
        const double right_edge = 0.5;

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
                    
                    if (x_cord == right_edge)
                    {
                      cell->face(face)->set_boundary_id(2);
                      square_edge = true;
                    }

                    // left edge
                    if (x_cord == left_edge)
                    {
                      square_edge = true;  
                      cell->face(face)->set_boundary_id(0);
                    }

                    // this is the bottom wall
                    if (y_cord == left_edge)
                    {
                        square_edge = true;
                      cell->face(face)->set_boundary_id(1);
                    }

                    // top edge, This is the second wall
                     if (y_cord == right_edge)
                     {
                         square_edge = true;
                        cell->face(face)->set_boundary_id(3);
                     }

                     // if the square_edge has not been found then the face is at the circular region
                     if(square_edge==false)
                        {
                            cell->face(face)->set_boundary_id(4);    
                            cell->face(face)->set_manifold_id(100);
                        }
                        
                   }
        }

  }

        
