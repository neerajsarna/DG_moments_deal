namespace Test_ForceType
{
	using namespace dealii;

	// we now test the force value arising from systemA
	TEST(ForceG56,HandlesForceG56)
	{
		const unsigned int dim = 2;
		std::string folder_name = "../system_matrices/";

		Constants::Base_Constants constants(input_file);
		G56::G56<dim> G56(constants.constants,folder_name); 

		// only do the computation if we have a particular force type and a particular value of alpha
		if (G56.constants.force_type == type3 )
		{
			std::vector<Point<dim>> p(1);

		// testing for forcetype1
			std::vector<Vector<double>> value(1,Vector<double>(G56.constants.nEqn));
			
		// x coordinate
			p[0][0] = 0.5;

		// y coordinate
			p[0][1] = 0.3;

			G56.source_term(p,value);
		// the exact value has been taken from mathematica file
			ASSERT_NEAR(value[0](dim + 1),-0.0734847,1e-5);			
		}

	}
}