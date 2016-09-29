namespace Test_Solver
{
	using namespace dealii;

	TEST(MeshGeneration,HandlesMeshGeneration)
	{
		const unsigned int dim = 2;
		ASSERT_EQ(dim,2) << "3D not implemented" << std::endl;

		std::string input_file = "../test_input_files/input1.in";
		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		SystemA::SystemA<dim> systemA(constants.constants,folder_name);

		FEM_Solver::Base_Solver<dim> base_solver("mesh_file_name",
											 "grid",
											 constants.constants,
											 &systemA);


		base_solver.print_grid();
		EXPECT_NE(base_solver.triangulation.n_active_cells(),0);


	}
}