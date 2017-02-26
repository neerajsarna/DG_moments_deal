namespace Test_BoundaryHandler_Odd
{
	using namespace dealii;

	TEST(BCrhsOddG20Square,HandlingBCrhsOddG20Square)
	{
		const unsigned int dim = 2;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		G20::G20<dim> G20(constants.constants,folder_name);		

		std::cout << "Tangential velocity " << constants.constants.vt4 << std::endl;

		const double theta0 = 1.0;
		const double theta1 = 1.0;
		const double theta2 = 2.0;
		const double theta3 = 2.0;
		const double theta4 = 1.5;
	
		const double vn0 = 0.0;
		const double vn1 = 0.0;
		const double vn2 = 0.0;
		const double vn3 = 0.0;
		const double vn4 = 0.0;

		const double vt0 = 0.0;
		const double vt1 = 0.0;
		const double vt2 = 0.0;
		const double vt3 = 0.0;
		const double vt4 = -1.0;

		Vector<double> bc_rhs(G20.constants.nBC);
		Vector<double> bc_rhs_manuel(G20.constants.nBC);


		
			// first we consider the left boundary
			Tensor<1,dim> p;
			Tensor<1,dim> normal_vector;

			p[0] = 0.5;
			p[1] = 0.6;


			bc_rhs_manuel(0) = vn0;
			bc_rhs_manuel(1) = -1/sqrt(M_PI) *vt0 ;
			bc_rhs_manuel(2) = 4 *  sqrt(2/(15 * M_PI)) * theta0 * sqrt(3.0/2.0);
			bc_rhs_manuel(3) = 0.106385 * theta0 * sqrt(3.0/2.0);
			bc_rhs_manuel(4) = -0.0531923 * theta0 * sqrt(3.0/2.0);


			G20.build_BCrhs(p,normal_vector,bc_rhs,0);

			for (unsigned int i = 0 ; i < bc_rhs.size() ; i++)
				EXPECT_NEAR(bc_rhs(i),bc_rhs_manuel(i),1e-5) << "value of i " << i ;

			// // bottom wall
			p[0] = 0;
			p[1] = -0.5;


			bc_rhs_manuel(0) = vn1;
			bc_rhs_manuel(1) = -1/sqrt(M_PI) *vt1 ;
			bc_rhs_manuel(2) = 4 *  sqrt(2/(15 * M_PI)) * theta1 * sqrt(3.0/2.0);
			bc_rhs_manuel(3) = 0.106385 * theta1 * sqrt(3.0/2.0);
			bc_rhs_manuel(4) = -0.0531923 * theta1 * sqrt(3.0/2.0);

			G20.build_BCrhs(p,normal_vector,bc_rhs,1);

			for (unsigned int i = 0 ; i < bc_rhs.size() ; i++)
				EXPECT_NEAR(bc_rhs(i),bc_rhs_manuel(i),1e-5) << "value of i " << i ;

			//RIGHT WALL

			bc_rhs_manuel(0) = vn2;
			bc_rhs_manuel(1) = -1/sqrt(M_PI) *vt2;
			bc_rhs_manuel(2) = 4 *  sqrt(2/(15 * M_PI)) * theta2 * sqrt(3.0/2.0);
			bc_rhs_manuel(3) = 0.106385 * theta2 * sqrt(3.0/2.0);
			bc_rhs_manuel(4) = -0.0531923 * theta2 * sqrt(3.0/2.0);

			G20.build_BCrhs(p,normal_vector,bc_rhs,2);

			for (unsigned int i = 0 ; i < bc_rhs.size() ; i++)
				EXPECT_NEAR(bc_rhs(i),bc_rhs_manuel(i),1e-5) << "value of i " << i ;


			// TOP WALL

			bc_rhs_manuel(0) = vn3;
			bc_rhs_manuel(1) = -1/sqrt(M_PI) *vt3;
			bc_rhs_manuel(2) = 4 *  sqrt(2/(15 * M_PI)) * theta3* sqrt(3.0/2.0);
			bc_rhs_manuel(3) = 0.106385 * theta3 * sqrt(3.0/2.0);
			bc_rhs_manuel(4) = -0.0531923 * theta3 * sqrt(3.0/2.0);

			G20.build_BCrhs(p,normal_vector,bc_rhs,3);

			for (unsigned int i = 0 ; i < bc_rhs.size() ; i++)
				EXPECT_NEAR(bc_rhs(i),bc_rhs_manuel(i),1e-5) << "value of i " << i ;		

			// The inner cylinder

			bc_rhs_manuel(0) = vn4;
			bc_rhs_manuel(1) = -1/sqrt(M_PI) *vt4;
			bc_rhs_manuel(2) = 4 *  sqrt(2/(15 * M_PI)) * theta4 * sqrt(3.0/2.0);
			bc_rhs_manuel(3) = 0.106385 * theta4 * sqrt(3.0/2.0);
			bc_rhs_manuel(4) = -0.0531923 * theta4 * sqrt(3.0/2.0);

			G20.build_BCrhs(p,normal_vector,bc_rhs,4);

			for (unsigned int i = 0 ; i < bc_rhs.size() ; i++)
				EXPECT_NEAR(bc_rhs(i),bc_rhs_manuel(i),1e-5) << "value of i " << i ;
	}	


	TEST(BCrhsInflowG20,HandlingBCrhsInflowG20)
	{
		const unsigned int dim = 2;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		G20::G20<dim> G20(constants.constants,folder_name);		

		std::cout << "Tangential velocity " << constants.constants.vt4 << std::endl;

 		double rho101 = 2.0;
	 	double theta101 = 1.0;
 		double vn101 = -1;
 		double vt101 = 1;

 
 		double rho102 = 1.0;
 		double theta102 = 1.0;
 		double vn102 = 1;
 		double vt102 = 1;
		
		// conversion to our variables
		double w_theta101 = -sqrt(3.0/2.0)*theta101;
		double w_theta102 = -sqrt(3.0/2.0)*theta102;

		Vector<double> bc_rhs(G20.constants.nBC);
		Vector<double> bc_rhs_manuel(G20.constants.nBC);


			// first we consider the boundary with bid=101
			Tensor<1,dim> p;
			Tensor<1,dim> normal_vector;

			p[0] = 0.5;
			p[1] = 0.6;

			bc_rhs_manuel(0) = vn101 + 1/sqrt(3*M_PI) * w_theta101-sqrt(2/M_PI) * rho101;
			bc_rhs_manuel(1) = -1/sqrt(M_PI) * vt101 ;
			bc_rhs_manuel(2) = -7/sqrt(30*M_PI) * w_theta101-(1/sqrt(5*M_PI)) * rho101;
			bc_rhs_manuel(3) = -sqrt(2/M_PI)/5.0 * w_theta101+2/(5.*sqrt(3*M_PI)) * rho101 ;
			bc_rhs_manuel(4) = 1/(5.*sqrt(2*M_PI)) * w_theta101-1/(5.*sqrt(3*M_PI)) * rho101;


			G20.bcrhs_inflow.BCrhs(p,normal_vector,bc_rhs,101);

			for (unsigned int i = 0 ; i < bc_rhs.size() ; i++)
				EXPECT_NEAR(bc_rhs(i),bc_rhs_manuel(i),1e-5) << "value of i " << i ;

			// consider bid = 102
			p[0] = 0;
			p[1] = -0.5;


			bc_rhs_manuel(0) = vn102 + 1/sqrt(3*M_PI) * w_theta102 - sqrt(2/M_PI) * rho102;
			bc_rhs_manuel(1) = -1/sqrt(M_PI) * vt102 ;
			bc_rhs_manuel(2) = -7/sqrt(30*M_PI) * w_theta102-(1/sqrt(5*M_PI)) * rho102;
			bc_rhs_manuel(3) = -sqrt(2/M_PI)/5.0 * w_theta102+2/(5.*sqrt(3*M_PI)) * rho102;
			bc_rhs_manuel(4) = 1/(5.*sqrt(2*M_PI)) * w_theta102-1/(5.*sqrt(3*M_PI)) * rho102;

			G20.bcrhs_inflow.BCrhs(p,normal_vector,bc_rhs,102);

			for (unsigned int i = 0 ; i < bc_rhs.size() ; i++)
				EXPECT_NEAR(bc_rhs(i),bc_rhs_manuel(i),1e-5) << "value of i " << i ;
	}	
}
