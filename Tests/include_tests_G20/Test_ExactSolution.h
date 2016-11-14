namespace Test_ExactSolution
{
	using namespace dealii;

	TEST(ExactSolutionPoissonHeatG26,HandlesExactSolutionG26)
	{
		const unsigned int dim = 2;
		std::string folder_name = "../system_matrices/";

		Constants::Base_Constants constants(input_file);
		G20::G20<dim> G20(constants.constants,folder_name);
		ExactSolution::G20_PoissonHeat<dim>  G20_PoissonHeat(constants.constants,G20.base_tensorinfo.S_half);

		Point<dim> p;
		Vector<double> value_manuel(constants.constants.nEqn);
		Vector<double> value(constants.constants.nEqn);

		const unsigned int ID_theta = constants.constants.variable_map.find("theta")->second;
		const unsigned int ID_stress = constants.constants.variable_map.find("sigmayy")->second;
		const unsigned int ID_heat = constants.constants.variable_map.find("qy")->second;

		// The value of the exact solution have been obtained from Mathematica

		p(0) = 0.6;
		p(1) = 0.6;

		if (constants.constants.mesh_type == periodic_square)
		{
			
			ASSERT_NEAR(G20.constants.alpha,0.816496580927726,1e-5);

			if (fabs(constants.constants.tau - 0.1) < 1e-5)
			{
				
				value_manuel[ID_theta] = -1.22478;
				value_manuel[ID_stress] = -0.00705621;
			}

			if (fabs(constants.constants.tau - 0.3) < 1e-5)
			{
				
				value_manuel[ID_theta] = -1.247;
				value_manuel[ID_stress] = 0.000254077;
				 
			}
		}



		// fix the solution for a symmetric system
		MatrixOpt::Base_MatrixOpt matrix_opt;
		value_manuel = matrix_opt.Sparse_matrix_dot_Vector(G20.base_tensorinfo.S_half,value_manuel);

		G20_PoissonHeat.vector_value(p,value);

		for (int i = 0 ; i < constants.constants.nEqn; i++)
			EXPECT_NEAR(value(i),value_manuel(i),1e-5) << "Failing at " << i;

	}
}