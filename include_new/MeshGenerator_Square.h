	// generate a ring shaped mesh internally
	template<int dim>
	void
	Base_MeshGenerator<dim>::mesh_internal_square(const unsigned int parts_x,const unsigned int parts_y)
	{
            triangulation.clear();
            Point<dim> p1;
            Point<dim> p2;
            std::vector<unsigned int > repetitions(dim);

            p1(0) = constants.xl;
            p1(1) = constants.yb;

            p2(0) = constants.xr;
            p2(1) = constants.yt;

            repetitions[0] = parts_x;
            repetitions[1] = parts_y;

            //The diagonal of the rectangle is the time joining p1 and p2
            GridGenerator::subdivided_hyper_rectangle(triangulation,
                                                      repetitions,
                                                       p1,p2);


            AssertThrow(constants.mesh_type == square_domain,ExcMessage("Incorrect mesh type"));

            switch(constants.problem_type)
            {
                case periodic:
                {
                    set_periodic_bid();
                    break;
                }

                case heat_conduction:
                {
                    set_square_bid();
                    break;
                }

                case lid_driven_cavity:
                {

                    set_square_bid();
                    break;
                }
                    // in this case the boundary id has to be developed by gmsh
                case inflow_outflow:
                {
                    std::cout << "Prescribing the boundary ids through gmsh" << std::endl;
                    break;
                }

                default:
                {
                    AssertThrow(1 == 0, ExcMessage("Should not have reached here"));
                    break;
                }

            }
            
	}
    
    // set boundary ids for all the walls involved in a square
    template<int dim>
    void
    Base_MeshGenerator<dim>::set_square_bid() const
    {
        typename Triangulation<dim>::cell_iterator cell = triangulation.begin(),
                                                    endc = triangulation.end();

        for (; cell != endc ; cell++)
        {
                for (unsigned int face = 0 ; face < GeometryInfo<dim>::faces_per_cell ; face++)
              
                  if (cell->face(face)->at_boundary())
                  { 
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
                      cell->face(face)->set_boundary_id(2);

                    // left edge
                    if (x_cord == constants.xl)
                      cell->face(face)->set_boundary_id(0);

                    // this is the bottom wall
                    if (y_cord == constants.yb)
                      cell->face(face)->set_boundary_id(1);

                    // top edge, This is the second wall
                    if (y_cord == constants.yt)
                      cell->face(face)->set_boundary_id(3);
                   }
        }
  }
