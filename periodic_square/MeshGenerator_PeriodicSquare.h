
	template<int dim>
	void
	Base_MeshGenerator<dim>::set_periodic_bid() const
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
                    double x_cord = cell->face(face)->center()(0);
                    double y_cord = cell->face(face)->center()(1);

                    // the periodic faces get the id 100 and 101
                    // right edge
                    if (x_cord == right_edge)
                      cell->face(face)->set_boundary_id(100);

                    // left edge
                    if (x_cord == left_edge)
                      cell->face(face)->set_boundary_id(101);

                    // This is the first wall
                    if (y_cord == left_edge)
                      cell->face(face)->set_boundary_id(0);

                    // top edge, This is the second wall
                    if (y_cord == right_edge)
                      cell->face(face)->set_boundary_id(1);
                   }
        }
  }


