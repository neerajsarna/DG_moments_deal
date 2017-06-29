template<int dim>
void
Base_MeshGenerator<dim>::mesh_internal_box()
{
	 const unsigned int repetitions = 10;
	 GridGenerator::subdivided_hyper_cube(triangulation,repetitions);

	 // set the boundary id for the cube
	 set_bid_cube_inflow();
}

template<int dim>
void
Base_MeshGenerator<dim>::set_bid_cube_inflow()
{
	const double half_edge = 1.0;				// edge length of half of the cube
	const double left_plane = 0;		    	//x coordinate, yz plane 
	const double right_plane = 1;			    // x coordinate, yz plane
	const double bottom_plane = 0; 			   // y coordinate, xz plane
	const double top_plane = 1;				   // y coordinate, xz plane
	const double front_plane = 15.0;		  // z coordinate, xy plane
	const double back_plane = 0.0; 			  // z coordinate, xy plane


	typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(),
													  endc = triangulation.end();


	for (; cell != endc ; cell++)
	{
		for (unsigned int face = 0 ; face < GeometryInfo<dim>::faces_per_cell ; face++)
			if (cell->face(face)->at_boundary())
			{ 
				bool square_edge = false;

				// x y and z coordinates of the face center
				double x_cord = cell->face(face)->center()(0);
				double y_cord = cell->face(face)->center()(1);
				double z_cord = cell->face(face)->center()(2);

				if (fabs(x_cord - left_plane) < 1e-16)
					cell->face(face)->set_boundary_id(101);


				if (fabs(y_cord - bottom_plane) < 1e-16)
					cell->face(face)->set_boundary_id(50);
					

				if (fabs(x_cord - right_plane) < 1e-16)
					cell->face(face)->set_boundary_id(102);


				if (fabs(y_cord - top_plane) < 1e-16)
					cell->face(face)->set_boundary_id(50);


				// specular walls on the front and on the back
				if (fabs(z_cord - back_plane) < 1e-16)
					cell->face(face)->set_boundary_id(50);
					

				if (fabs(z_cord - front_plane) < 1e-16)
					cell->face(face)->set_boundary_id(50);
					

			}

	}
}