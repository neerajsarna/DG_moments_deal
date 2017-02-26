namespace TestMeshGenerator
{
	using namespace dealii;

	// TEST(MeshGenerator,HandlesMeshGeneration)
	// {
	// 	const unsigned int dim = 2;
	// 	ASSERT_EQ(dim,2) << "3D not implemented" << std::endl;

	// 	std::string folder_name = "../system_matrices/";
	// 	Constants::Base_Constants constants(input_file);
	// 	std::string mesh_file = "mesh_file";
	// 	std::string output_file = "testing_grid";

	// 	MeshGenerator::Base_MeshGenerator<dim> mesh_gen(mesh_file,output_file,constants.constants);
	// 	std::cout << "#Cells " << mesh_gen.triangulation.n_active_cells();
	// 	mesh_gen.print_grid(0);

	// 	typename Triangulation<dim>::active_cell_iterator cell = mesh_gen.triangulation.begin_active(),
	// 													   endc = mesh_gen.triangulation.end();


	// 	for (; cell != endc ; cell++)
	// 		for (unsigned int face = 0 ; face < GeometryInfo<dim>::faces_per_cell ; face++)
	// 			if (cell->face(face)->at_boundary())
	// 			{
	// 				const double x_coord = cell->face(face)->center()(0);
	// 				const double y_coord = cell->face(face)->center()(1);
	// 				bool entered = false;

	// 				if (cell->face(face)->boundary_id() == 0)
	// 				{
	// 					entered = true;
	// 					EXPECT_NEAR(x_coord,constants.constants.xl,1e-5);
	// 				}

	// 				if (cell->face(face)->boundary_id() == 1)
	// 				{
	// 					entered = true;
	// 					EXPECT_NEAR(y_coord,constants.constants.yb,1e-5);
	// 				}

	// 				if (cell->face(face)->boundary_id() == 2)
	// 				{
	// 					entered = true;
	// 					EXPECT_NEAR(x_coord,constants.constants.xr,1e-5);
	// 				}

	// 				if (cell->face(face)->boundary_id() == 3)
	// 				{
	// 					entered = true;
	// 					EXPECT_NEAR(y_coord,constants.constants.yt,1e-5);
	// 				}

	// 				// if (!entered)
	// 				// 	std::cout << "Distance from the center: "<<cell->face(face)->center().norm() << std::endl;


	// 			}
	// }


	TEST(MeshGenerator,HandlesMeshGeneration)
	{
		const unsigned int dim = 2;
		ASSERT_EQ(dim,2) << "3D not implemented" << std::endl;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		std::string mesh_file = "mesh_file";
		std::string output_file = "Square_Cavity";

		MeshGenerator::Base_MeshGenerator<dim> mesh_gen(mesh_file,output_file,constants.constants);
		std::cout << "#Cells " << mesh_gen.triangulation.n_active_cells();
		mesh_gen.print_grid(0);

		// for (; cell != end_c ; cell++)
		// 	for (unsigned int face = 0 ; face < GeometryInfo<dim>::faces_per_cell ;face++)
		// 		if (cell->face(face)->at_boundary())
		// 			printf("%d\n",cell->face(face)->boundary_id()); 
	}
}