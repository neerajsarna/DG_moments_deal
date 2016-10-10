namespace Test_ForceType
{
	using namespace dealii;

	// checks the forcing for un symmetric system
	TEST(ForceTypeUnSymm,HandlesForcingUnSymm)
	{
		// input file
		const unsigned int dim = 2;
		ASSERT_LT(dim,3)<< "3D not implemented " << std::endl;
		
		// first we need to read all the constants
		Constants::Base_Constants constants(input_file);

		ForceType::ForceType1<dim> force_type1(constants.constants);
		ForceType::ForceType2<dim> force_type2(constants.constants);
		ForceType::ForceType3<dim> force_type3(constants.constants);

		// We now consider a single point and compute the value of the forces at that point
		std::vector<Point<dim>> p(1);

		// testing for forcetype1
		std::vector<Vector<double>> value1(1,Vector<double>(force_type1.constants.nEqn));

		const double exact_value1 = 0.195;
	
		const double force_factor = 1;
		// x coordinate
		p[0][0] = 0.5;

		// y coordinate
		p[0][1] = 0.8;

		force_type1.source_term(p,value1,force_factor);

		// first we need to check whether the values of A0, A1 and A2 correspond to 
		// the exact solution or not 
		if (fabs(force_type1.constants.A0 - 0.0) < 1e-5 &&
			fabs(force_type1.constants.A1-0.2) < 1e-5 &&
			fabs(force_type1.constants.A2-0.1) < 1e-5)
		{

		// exact value obtained from mathematica
		EXPECT_LT(fabs(value1[0](0)-exact_value1),1e-5);

		// we also need to check for the other values of value
		for (unsigned int i =0 ; i < value1.size() ; i++)
			for (unsigned int j = 0 ; j < value1[i].size() ; i++)
				if (j != 0)
					EXPECT_FLOAT_EQ(value1[i](j),0);			
		}




		// testing for force type2
		std::vector<Vector<double>> value2(1,Vector<double>(force_type2.constants.nEqn));

		const double exact_value2 = -2.42555;

		force_type2.source_term(p,value2,force_factor);

		// first we need to check whether the values of A0, A1 and A2 correspond to 
		// the exact solution or not 
		if (fabs(force_type2.constants.A0 - 0.0) < 1e-5 &&
			fabs(force_type2.constants.A1-0.2) < 1e-5 &&
			fabs(force_type2.constants.A2-0.1) < 1e-5 &&
			fabs(force_type2.constants.tau - 0.1) < 1e-5)
		{
			// exact value obtained from mathematica
			ASSERT_LT(fabs(value2[0](0)-exact_value2),1e-5);


			// we also need to check for the other values of value
			for (unsigned int i =0 ; i < value2.size() ; i++)
				for (unsigned int j = 0 ; j < value2[i].size() ; i++)
					if (j != 0)
						EXPECT_FLOAT_EQ(value2[i](j),0);
		}




		// testing for the third force type
		std::vector<Vector<double>> value3(1,Vector<double>(force_type3.constants.nEqn));

		const double exact_value3 = -0.64;

		force_type3.source_term(p,value3,force_factor);

		// first we need to check whether the values of A0, A1 and A2 correspond to 
		// the exact solution or not 
		if (fabs(force_type3.constants.alpha-1.0) < 1e-5)
		{
				// exact value obtained from mathematica
				EXPECT_LT(fabs(value3[0](3)-exact_value3),1e-5);

				// we also need to check for the other values of value
				for (unsigned int i =0 ; i < value3.size() ; i++)
					for (unsigned int j = 0 ; j < value3[i].size() ; j++)
						if (j != 3)
							EXPECT_FLOAT_EQ(value3[i](j),0);
		}


	}

	// we now test the force value arising from systemA
	TEST(ForceSystemA,HandlesForceSystemA)
	{
		const unsigned int dim = 2;
		std::string folder_name = "../system_matrices/";

		Constants::Base_Constants constants(input_file);
		SystemA::SystemA<dim> systemA(constants.constants,folder_name); 

		if (systemA.constants.force_type == type1)
		{
			std::vector<Point<dim>> p(1);

		// testing for forcetype1
			std::vector<Vector<double>> value(1,Vector<double>(systemA.constants.nEqn));
			
		// x coordinate
			p[0][0] = 0.5;

		// y coordinate
			p[0][1] = 0.8;

			systemA.source_term(p,value);
		// the exact value has been taken from mathematica file
			ASSERT_NEAR(value[0](0),0.195,1e-5);			
		}

	}
}