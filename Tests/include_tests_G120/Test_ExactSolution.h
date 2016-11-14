namespace Test_ExactSolution
{
	using namespace dealii;

	TEST(ExactSolutionPoissonHeatG120,HandlesExactSolutionG120)
	{
		const unsigned int dim = 2;
		std::string folder_name = "../system_matrices/";

		Constants::Base_Constants constants(input_file);
		G120::G120<dim> G120(constants.constants,folder_name);
		ExactSolution::G120_PoissonHeat<dim>  G120_PoissonHeat(constants.constants,G120.base_tensorinfo.S_half);

		Point<dim> p;
		Vector<double> value_manuel(constants.constants.nEqn);
		Vector<double> value(constants.constants.nEqn);

		const unsigned int ID_theta = constants.constants.variable_map.find("theta")->second;
		const unsigned int ID_stress = constants.constants.variable_map.find("sigmayy")->second;
		const unsigned int ID_heat = constants.constants.variable_map.find("qy")->second;

		// The value of the exact solution has been obtained from Mathematica

		p(0) = 0.6;
		p(1) = 0.6;

		if (constants.constants.mesh_type == periodic_square)
		{
			
			ASSERT_NEAR(G120.constants.alpha,0.816496580927726,1e-5);

			if (fabs(constants.constants.tau - 0.1) < 1e-5)
			{
				
				value_manuel[ID_theta] = -1.19819;	
				value_manuel[ID_stress] = 0.0121628;
				value_manuel[ID_heat] = -0.0455368;
			}

			if (fabs(constants.constants.tau - 0.3) < 1e-5)
			{
				
				value_manuel[ID_theta] = -1.24145;	
				value_manuel[ID_stress] = 0.00140134;
				value_manuel[ID_heat] = -0.0455368;
				 
			}
		}



		// fix the solution for a symmetric system
		MatrixOpt::Base_MatrixOpt matrix_opt;
		value_manuel = matrix_opt.Sparse_matrix_dot_Vector(G120.base_tensorinfo.S_half,value_manuel);

		G120_PoissonHeat.vector_value(p,value);

		for (int i = 0 ; i < constants.constants.nEqn; i++)
			EXPECT_NEAR(value(i),value_manuel(i),1e-5) << "Failing at " << i;

	}
}