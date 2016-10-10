namespace Test_ExactSolution
{
	using namespace dealii;

	TEST(ExactSolutionRingSystemA,HandlesExactSolutionRingSystemA)
	{
		const unsigned int dim = 2;
		std::string folder_name = "../system_matrices/";

		Constants::Base_Constants constants(input_file);
		SystemA::SystemA<dim> systemA(constants.constants,folder_name);
		ExactSolution::ExactSolution_SystemA_ring<dim>  exact_solution_systemA(constants.constants,systemA.base_tensorinfo.S_half);

		Point<dim> p;
		Vector<double> value_manuel(constants.constants.nEqn);
		Vector<double> value(constants.constants.nEqn);

		// The value of the exact solution have been obtained from Mathematica

		p(0) = 0.6;
		p(1) = 0.6;

		if (constants.constants.mesh_type == ring)
		{
			if (fabs(constants.constants.A0 - 0.0) < 1e-5 &&
				fabs(constants.constants.A1 - 0.2) < 1e-5 &&
				fabs(constants.constants.A2 - 0.1) < 1e-5 &&
				fabs(constants.constants.uW - constants.constants.tau) < 1e-5)
			{
				if (fabs(constants.constants.tau - 10.0) < 1e-5 )
				{
					value_manuel(0) = 0.0620844;
					value_manuel(1) = -1.07755;
					value_manuel(2) = 0.0492537;			
				}

				if (fabs(constants.constants.tau - 0.1) < 1e-5)
				{
					value_manuel(0) = 2.43003;
					value_manuel(1) = -0.0577167;
					value_manuel(2) = 0.00268024;
				}

				if (fabs(constants.constants.tau - 0.01) < 1e-5)
				{
					value_manuel(0) = 8.15317;
					value_manuel(1) = -0.0988336;
					value_manuel(2) = -0.0509269;
				}

				if (fabs(constants.constants.tau - 1.0) < 1e-5)
				{
					value_manuel(0) = 1.39944;
					value_manuel(1) = -0.0555114;
					value_manuel(2) = 0.115708;
				}


			}
			
		}



		MatrixOpt::Base_MatrixOpt matrix_opt;
		value_manuel = matrix_opt.Sparse_matrix_dot_Vector(systemA.base_tensorinfo.S_half,value_manuel);

		exact_solution_systemA.vector_value(p,value);

		for (int i = 0 ; i < constants.constants.nEqn; i++)
			EXPECT_NEAR(value(i),value_manuel(i),1e-5) << "Failing at " << i;

	}
}