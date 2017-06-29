namespace TestMeshGenerator
{
	using namespace dealii;

	TEST(MeshGenerator,HandlesMeshGeneration)
	{
		const unsigned int dim = 3;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		std::string output_file = "grid2_3D";


		MeshGenerator::Base_MeshGenerator<dim> 		mesh_gen(
							   							output_file,
							   							constants.constants_num.mesh_filename,
							   							constants.constants_num.mesh_type,
							   							constants.constants_num.problem_type,
							   							constants.constants_num.part_x,
							   							constants.constants_num.part_y,
							   							constants.constants_num.initial_refinement);



		//mesh_gen.triangulation.refine_global(constants.constants_num.initial_refinement);
		std::cout << "#Cells " << mesh_gen.triangulation.n_active_cells();
		mesh_gen.print_mesh_info();
		mesh_gen.print_grid(0);
		mesh_gen.print_grid_points();

	}

	TEST(DISABLED_MeshGenerator,HandlesMeshGeneration)
	{
		const unsigned int dim = 3;

		CylindricalManifold<dim> boundary_description(2);
		Triangulation<dim> triangulation;

		GridGenerator::hyper_cube_with_cylindrical_hole (triangulation,
                                                          0.5,
                                                           2.0,
                                                            1.0, 
                                                            1);

		typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(),
														 endc = triangulation.end();



		for (; cell != endc ; cell++)
			for (unsigned int face = 0 ; face < GeometryInfo<dim>::faces_per_cell ; face++)
				if (cell->face(face)->at_boundary())
					for (unsigned int v = 0 ; v < GeometryInfo<dim>::vertices_per_face ; v++)
					{
						Point<dim> projection(0,0,cell->face(face)->vertex(v)(2));

						if (fabs(projection.distance(cell->face(face)->vertex(v))-0.5) < 1e-16)
							cell->face(face)->set_boundary_id(1);
					}


	// cell = triangulation.begin_active();


	// for (; cell != endc ; cell++)
	// {
	// 	for (unsigned int face = 0 ; face < GeometryInfo<dim>::faces_per_cell ; face++)
	// 		if (cell->face(face)->at_boundary())
	// 		{ 
	// 			// x y and z coordinates of the face center
	// 			double x_cord = cell->face(face)->center()(0);
	// 			double y_cord = cell->face(face)->center()(1);
	// 			double z_cord = cell->face(face)->center()(2);


	// 			// specular walls on the front and on the back
	// 			if (fabs(z_cord - 0.0) < 1e-16)
	// 			{
	// 				cell->face(face)->set_boundary_id(50);
	// 			}

	// 			if (fabs(z_cord - 1.0) < 1e-16)
	// 			{
	// 				cell->face(face)->set_boundary_id(50);
	// 			}

	// 		}
	// }

	const double half_edge = 4.0;				// edge length of half of the cube
	const double left_plane = -half_edge;		//x coordinate, yz plane 
	const double right_plane = half_edge;		// x coordinate, yz plane
	const double bottom_plane = -half_edge; 	// y coordinate, xz plane
	const double top_plane = half_edge;			// y coordinate, xz plane
	const double front_plane = 1.0;				// z coordinate, xy plane
	const double back_plane = 0.0; 				// z coordinate, xy plane


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

			}

		}





		triangulation.set_manifold(1,boundary_description);

		triangulation.refine_global(2);
		std::string file_for_grid = "grid2_3D";
		std::ofstream out (file_for_grid.c_str());
		GridOut grid_out;
		grid_out.write_eps (triangulation, out);
	}
}