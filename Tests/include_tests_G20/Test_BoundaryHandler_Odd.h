namespace Test_BoundaryHandler_Odd
{
	using namespace dealii;



	TEST(BCrhsOddG20Square,HandlingBCrhsOddG20Square)
	{
		const unsigned int dim = 2;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		G20::G20<dim> G20(constants.constants,folder_name);		

		const double theta0 = 1.0;
		const double theta1 = 1.0;
		const double theta2 = 2.0;
		const double theta3 = 2.0;

		Vector<double> bc_rhs(G20.constants.nEqn);
		Vector<double> bc_rhs_manuel(G20.constants.nEqn);

		Assert(constants.constants.bc_type == odd,ExcMessage("Incorrect boundary implementation for the current test"));

		if (constants.constants.bc_type == odd)		
		{
			// first we consider the left boundary
			Tensor<1,dim> p;
			Tensor<1,dim> normal_vector;

			p[0] = -0.5;
			p[1] = 0.0;


			bc_rhs_manuel(1) = 0;
			bc_rhs_manuel(5) = 0;
			bc_rhs_manuel(7) = 4 *  sqrt(2/(15 * M_PI)) * theta0 * sqrt(3.0/2.0);
			bc_rhs_manuel(9) = 0.106385 * theta0 * sqrt(3.0/2.0);
			bc_rhs_manuel(11) = -0.0531923 * theta0 * sqrt(3.0/2.0);


			G20.build_BCrhs(p,normal_vector,bc_rhs,0);

			for (unsigned int i = 0 ; i < bc_rhs.size() ; i++)
				EXPECT_NEAR(bc_rhs(i),bc_rhs_manuel(i),1e-5) << "value of i " << i ;

			// // bottom wall
			p[0] = 0;
			p[1] = -0.5;


			bc_rhs_manuel(1) = 0;
			bc_rhs_manuel(5) = 0;
			bc_rhs_manuel(7) = 4 *  sqrt(2/(15 * M_PI)) * theta1 * sqrt(3.0/2.0);
			bc_rhs_manuel(9) = 0.106385 * theta1 * sqrt(3.0/2.0);
			bc_rhs_manuel(11) = -0.0531923 * theta1 * sqrt(3.0/2.0);

			G20.build_BCrhs(p,normal_vector,bc_rhs,1);

			for (unsigned int i = 0 ; i < bc_rhs.size() ; i++)
				EXPECT_NEAR(bc_rhs(i),bc_rhs_manuel(i),1e-5) << "value of i " << i ;


			bc_rhs_manuel(1) = 0;
			bc_rhs_manuel(5) = 0;
			bc_rhs_manuel(7) = 4 *  sqrt(2/(15 * M_PI)) * theta2 * sqrt(3.0/2.0);
			bc_rhs_manuel(9) = 0.106385 * theta2 * sqrt(3.0/2.0);
			bc_rhs_manuel(11) = -0.0531923 * theta2 * sqrt(3.0/2.0);

			G20.build_BCrhs(p,normal_vector,bc_rhs,2);

			for (unsigned int i = 0 ; i < bc_rhs.size() ; i++)
				EXPECT_NEAR(bc_rhs(i),bc_rhs_manuel(i),1e-5) << "value of i " << i ;


			bc_rhs_manuel(1) = 0;
			bc_rhs_manuel(5) = 0;
			bc_rhs_manuel(7) = 4 *  sqrt(2/(15 * M_PI)) * theta3 * sqrt(3.0/2.0);
			bc_rhs_manuel(9) = 0.106385 * theta3 * sqrt(3.0/2.0);
			bc_rhs_manuel(11) = -0.0531923 * theta3 * sqrt(3.0/2.0);

			G20.build_BCrhs(p,normal_vector,bc_rhs,3);

			for (unsigned int i = 0 ; i < bc_rhs.size() ; i++)
				EXPECT_NEAR(bc_rhs(i),bc_rhs_manuel(i),1e-5) << "value of i " << i ;

		}
		

	}	
}
