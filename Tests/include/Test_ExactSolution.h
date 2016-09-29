namespace Test_ExactSolution
{
	using namespace dealii;

	TEST(ExactSolutionRingSystemA,HandlesExactSolutionRingSystemA)
	{
		const unsigned int dim = 2;
		std::string input_file = "../test_input_files/input1.in";
		std::string folder_name = "../system_matrices/";

		Constants::Base_Constants constants(input_file);
		ExactSolution::ExactSolution_SystemA_ring<dim>  exact_solution_systemA(constants.constants);

		Point<dim> p;
		Vector<double> value_manuel(constants.constants.nEqn);
		Vector<double> value(constants.constants.nEqn);

		p(0) = 0.6;
		p(1) = 0.6;

		if (fabs(constants.constants.tau - 10.0) < 1e-5)
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


		exact_solution_systemA.vector_value(p,value);

		for (int i = 0 ; i < constants.constants.nEqn; i++)
			EXPECT_NEAR(value(i),value_manuel(i),1e-5) << "Failing at " << i;

	}
}