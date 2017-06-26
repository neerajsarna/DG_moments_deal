namespace TestMeshGenerator
{
	using namespace dealii;

	TEST(MeshGenerator,HandlesMeshGeneration)
	{
		const unsigned int dim = 3;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		std::string output_file = "grid_box";


		MeshGenerator::Base_MeshGenerator<dim> 		mesh_gen(
							   							output_file,
							   							constants.constants_num.mesh_filename,
							   							constants.constants_num.mesh_type,
							   							constants.constants_num.problem_type,
							   							constants.constants_num.part_x,
							   							constants.constants_num.part_y,
							   							constants.constants_num.initial_refinement);

		std::cout << "#Cells " << mesh_gen.triangulation.n_active_cells();
		mesh_gen.print_mesh_info();
		mesh_gen.print_grid(0);

		typename Triangulation<dim>::active_cell_iterator cell = mesh_gen.triangulation.begin_active(),
														endc = mesh_gen.triangulation.end();


	const double left_plane = -2.0;		//x coordinate, yz plane 
	const double right_plane = 2.0;		// x coordinate, yz plane
	const double bottom_plane = -2.0; 	// y coordinate, xz plane
	const double top_plane = 2.0;		// y coordinate, xz plane
	const double front_plane = 1.0;		// z coordinate, xy plane
	const double back_plane = 0.0; 		// z coordinate, xy plane



		for (; cell != endc ; cell++)
		{
			for (unsigned int face = 0 ; face < GeometryInfo<dim>::faces_per_cell ; face++)
				if (cell->face(face)->at_boundary())
				{
					const unsigned int b_id = cell->face(face)->boundary_id();
					bool square_edge = false;

				// x y and z coordinates of the face center
					double x_cord = cell->face(face)->center()(0);
					double y_cord = cell->face(face)->center()(1);
					double z_cord = cell->face(face)->center()(2);

					if (fabs(x_cord - left_plane) < 1e-16)
					{
						square_edge = true;
						EXPECT_EQ(b_id,101);
					}


					if (fabs(y_cord - bottom_plane) < 1e-16)
					{
						square_edge = true;
						EXPECT_EQ(b_id,0);
					}


					if (fabs(x_cord - right_plane) < 1e-16)
					{
						square_edge = true;
						EXPECT_EQ(b_id,102);
					}

					if (fabs(y_cord - top_plane) < 1e-16)
					{
						square_edge = true;
						EXPECT_EQ(b_id,1);
					}

					if (fabs(z_cord - back_plane) < 1e-16)
					{
						square_edge = true;
						EXPECT_EQ(b_id,2);
					}

					if (fabs(z_cord - front_plane) < 1e-16)
					{
						square_edge = true;
						EXPECT_EQ(b_id,3);
					}

					if (square_edge == false)
						EXPECT_EQ(b_id,4);
					
				}
		}

	}
}