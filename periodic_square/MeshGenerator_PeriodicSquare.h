	// generate a ring shaped mesh internally
	template<int dim>
	void
	Base_MeshGenerator<dim>::mesh_internal_periodic_square()
	{
            Point<dim> p1;
            Point<dim> p2;
            std::vector<unsigned int > repetitions(dim);

            p1(0) = constants.xl;
            p1(1) = constants.yb;

            p2(0) = constants.xr;
            p2(1) = constants.yt;

            repetitions[0] = constants.part_x;
            repetitions[1] = constants.part_y;

            //The diagonal of the rectangle is the time joining p1 and p2
            GridGenerator::subdivided_hyper_rectangle(triangulation,
                                      	              repetitions,
                                        	            p1,
                                            	        p2);

            // set the boundary id
            set_periodic_bid();
            
	}

	// read a ring shaped mesh generated from gmsh
	template<int dim>
	void
	Base_MeshGenerator<dim>::mesh_gmsh_periodic_square()
	{

		// first we read the mesh
		gridin.attach_triangulation(triangulation);
        std::ifstream f(mesh_file_name);
        gridin.read_msh(f);

        // then we se the boudnary id
        set_periodic_bid();
	}

	template<int dim>
	void
	Base_MeshGenerator<dim>::set_periodic_bid() const
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

                    // the periodic faces get the id 100 and 101
                    // right edge
                    if (x_cord == constants.xr)
                      cell->face(face)->set_boundary_id(100);

                    // left edge
                    if (x_cord == constants.xl)
                      cell->face(face)->set_boundary_id(101);

                    // This is the first wall
                    if (y_cord == constants.yb)
                      cell->face(face)->set_boundary_id(0);

                    // top edge, This is the second wall
                    if (y_cord == constants.yt)
                      cell->face(face)->set_boundary_id(1);
                   }
        }
  }


