// we do not need to specialize this particular routine since we have already specialized develop_mesh routine
template<int dim>
void 
Base_MeshGenerator<dim>::set_bid_box_cylinder()
{
	AssertDimension(dim,3);
	const double left_plane = -2.0;		//x coordinate, yz plane 
	const double right_plane = 2.0;		// x coordinate, yz plane
	const double bottom_plane = -2.0; 	// y coordinate, xz plane
	const double top_plane = 2.0;		// y coordinate, xz plane
	const double front_plane = 1.0;		// z coordinate, xy plane
	const double back_plane = 0.0; 		// z coordinate, xy plane


	// we only iterate over the active cells
	typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(),endc = triangulation.end();


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
				{
					// inflow surface
					cell->face(face)->set_boundary_id(101);
					square_edge = true;
				}

				if (fabs(y_cord - bottom_plane) < 1e-16)
				{
					cell->face(face)->set_boundary_id(0);
					square_edge = true;
				}

				if (fabs(x_cord - right_plane) < 1e-16)
				{
					cell->face(face)->set_boundary_id(102);
					square_edge = true;
				}

				if (fabs(y_cord - top_plane) < 1e-16)
				{
					cell->face(face)->set_boundary_id(1);
					square_edge = true;
				}

				// specular walls on the front and on the back
				if (fabs(z_cord - back_plane) < 1e-16)
				{
					cell->face(face)->set_boundary_id(50);
					square_edge = true;
				}

				if (fabs(z_cord - front_plane) < 1e-16)
				{
					cell->face(face)->set_boundary_id(50);
					square_edge = true;
				}

				if (square_edge == false)
				{
					cell->face(face)->set_boundary_id(4);
					square_edge = true;
				}
			}

	}

	// now we set the spherical manifold
	triangulation.set_manifold(4,boundary);
}



template<int dim>
void
Base_MeshGenerator<dim>::mesh_internal_box_cylinder()
{
			AssertDimension(dim,3);
	        const double inner_radius = 0.5;
            double outer_radius = 2.0;
            const double length_in_z = 1.0;
            const int repetitions_in_z = 1;

            //The diagonal of the rectangle is the line joining p1 and p2
            GridGenerator::hyper_cube_with_cylindrical_hole (triangulation,
                                                            inner_radius,
                                                             outer_radius,
                                                            length_in_z, 
                                                            repetitions_in_z);


            // we only solve the inflow problem with the present geometry
            switch(problem_type)
            {
            	case inflow_outflow:
            	{
            		set_bid_box_cylinder();
            		break;
            	}

            	case heat_conduction:
            	case lid_driven_cavity:
            	case periodic:
            	default:
            	{
            		Assert(1 == 0 ,ExcMessage("Test Case not implemented for this boundary"));
            	}

            }

            //triangulation.refine_global(initial_refinement);
}
