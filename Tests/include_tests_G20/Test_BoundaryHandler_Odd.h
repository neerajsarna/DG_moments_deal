namespace Test_BoundaryHandler_Odd
{
	using namespace dealii;

	TEST(BCrhsWallG20Square,HandlingBCrhsWallG20Square)
	{
		const unsigned int dim = 1;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		G20::G20<dim> G20(constants.constants,folder_name);		


		if (dim == 2)
		{
		const double  theta0 = 1.0;
		const double  theta1 = 1.0;
		const double  theta2 = 2.0;
		const double  theta3 = 2.0;
		const double  theta4 = 0.0;


		const double  vx0 = 0.0;
		const double  vx1 = 0.0;
		const double  vx2 = 0.0;
		const double  vx3 = 0.0;
		const double  vx4 = 0.0;


		const double  vy0 = 0.0;
		const double  vy1 = 0.0;
		const double  vy2 = 0.0;
		const double  vy3 = 0.0;
		const double  vy4 = 0.0;

		// first we consider the left boundary
		Tensor<1,dim> p;
		Tensor<1,dim> normal_vector;

		p[0] = 0.5;
		p[1] = 0.6;

		normal_vector[0] = -1/sqrt(2);
		normal_vector[1] = 1/sqrt(2);

		Vector<double> bc_rhs(G20.constants.nBC);
		Vector<double> bc_rhs_manuel(G20.constants.nBC);
		FullMatrix<double> projector(2,2);

		projector(0,0) = normal_vector[0];
		projector(0,1) = normal_vector[1];
		projector(1,0) = -normal_vector[1];
		projector(1,1) = normal_vector[0];

		Vector<double> wall_velocity = G20.bcrhs_wall.local_wall_velocity(vx0,vy0,projector);



		bc_rhs_manuel(0) = wall_velocity(0);
		bc_rhs_manuel(1) = -1/sqrt(M_PI) * wall_velocity(1) ;
		bc_rhs_manuel(2) = 4 *  sqrt(2/(15 * M_PI)) * theta0 * sqrt(3.0/2.0);
		bc_rhs_manuel(3) = 0.106385 * theta0 * sqrt(3.0/2.0);
		bc_rhs_manuel(4) = -0.0531923 * theta0 * sqrt(3.0/2.0);


		

		G20.build_BCrhs(p,normal_vector,bc_rhs,0);

		for (unsigned int i = 0 ; i < bc_rhs.size() ; i++)
			EXPECT_NEAR(bc_rhs(i),bc_rhs_manuel(i),1e-5) << "value of i " << i ;

		wall_velocity = G20.bcrhs_wall.local_wall_velocity(vx1,vy1,projector);

		bc_rhs_manuel(0) = wall_velocity(0);
		bc_rhs_manuel(1) = -1/sqrt(M_PI) * wall_velocity(1) ;
		bc_rhs_manuel(2) = 4 *  sqrt(2/(15 * M_PI)) * theta1 * sqrt(3.0/2.0);
		bc_rhs_manuel(3) = 0.106385 * theta1 * sqrt(3.0/2.0);
		bc_rhs_manuel(4) = -0.0531923 * theta1 * sqrt(3.0/2.0);

		G20.build_BCrhs(p,normal_vector,bc_rhs,1);

		for (unsigned int i = 0 ; i < bc_rhs.size() ; i++)
			EXPECT_NEAR(bc_rhs(i),bc_rhs_manuel(i),1e-5) << "value of i " << i ;

			//RIGHT WALL
		wall_velocity = G20.bcrhs_wall.local_wall_velocity(vx2,vy2,projector);

		bc_rhs_manuel(0) = wall_velocity(0);
		bc_rhs_manuel(1) = -1/sqrt(M_PI) * wall_velocity(1);
		bc_rhs_manuel(2) = 4 *  sqrt(2/(15 * M_PI)) * theta2 * sqrt(3.0/2.0);
		bc_rhs_manuel(3) = 0.106385 * theta2 * sqrt(3.0/2.0);
		bc_rhs_manuel(4) = -0.0531923 * theta2 * sqrt(3.0/2.0);

		G20.build_BCrhs(p,normal_vector,bc_rhs,2);

		for (unsigned int i = 0 ; i < bc_rhs.size() ; i++)
			EXPECT_NEAR(bc_rhs(i),bc_rhs_manuel(i),1e-5) << "value of i " << i ;


		wall_velocity = G20.bcrhs_wall.local_wall_velocity(vx3,vy3,projector);
			// TOP WALL

		bc_rhs_manuel(0) = wall_velocity(0);
		bc_rhs_manuel(1) = -1/sqrt(M_PI) * wall_velocity(1);
		bc_rhs_manuel(2) = 4 *  sqrt(2/(15 * M_PI)) * theta3* sqrt(3.0/2.0);
		bc_rhs_manuel(3) = 0.106385 * theta3 * sqrt(3.0/2.0);
		bc_rhs_manuel(4) = -0.0531923 * theta3 * sqrt(3.0/2.0);

		G20.build_BCrhs(p,normal_vector,bc_rhs,3);

		for (unsigned int i = 0 ; i < bc_rhs.size() ; i++)
			EXPECT_NEAR(bc_rhs(i),bc_rhs_manuel(i),1e-5) << "value of i " << i ;					
		}

		if (dim == 1)
		{
			const double  theta0 = 1.0;
			const double  theta1 = 2.0;
			const double  theta2 = 3.0;
			const double  theta3 = 4.0;
			const double  theta4 = 4.0;


			double vx0 = 1.0;
			double vx1 = 2.0;
			double vx2 = 4.0;
			double vx3 = 5.0;
			double vx4 = 10.0;


			double vy0 = 5.0;
			double vy1 = 6.0;
			double vy2 = 7.0;
			double vy3 = 8.0;
			double vy4 = 9.0;

		// first we consider the left boundary
		Tensor<1,dim> p;
		Tensor<1,dim> normal_vector;

		p[0] = 0.5;

		normal_vector[0] = -1;
		

		Vector<double> bc_rhs(G20.constants.nBC);
		Vector<double> bc_rhs_manuel(G20.constants.nBC);
		
		double wall_velocity = -vx0;


		bc_rhs_manuel(0) = wall_velocity;
		bc_rhs_manuel(1) = -0.921318 * theta0 /sqrt(2) ;
		bc_rhs_manuel(2) = 0.412026 * theta0 * sqrt(1.0/2.0);


		G20.build_BCrhs(p,normal_vector,bc_rhs,0);

		for (unsigned int i = 0 ; i < bc_rhs.size() ; i++)
			EXPECT_NEAR(bc_rhs(i),bc_rhs_manuel(i),1e-5) << "value of i " << i ;

		normal_vector[0] = 1;
		wall_velocity = vx1;

		bc_rhs_manuel(0) = wall_velocity;
		bc_rhs_manuel(1) = -0.921318 * theta1 /sqrt(2) ;
		bc_rhs_manuel(2) = 0.412026 * theta1 * sqrt(1.0/2.0);

		G20.build_BCrhs(p,normal_vector,bc_rhs,1);

		for (unsigned int i = 0 ; i < bc_rhs.size() ; i++)
			EXPECT_NEAR(bc_rhs(i),bc_rhs_manuel(i),1e-5) << "value of i " << i ;
	
		}


		
	}	


	TEST(BCrhsInflowG20,HandlingBCrhsInflowG20)
	{
		const unsigned int dim = 1;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		G20::G20<dim> G20(constants.constants,folder_name);		

 		double rho101 = 2.0;
	 	double theta101 = 1.0;
 		double vx101 = -1;
 		double vy101 = 0;


 		double rho102 = 1.0;
 		double theta102 = 1.0;
 		double vx102 = 1;
 		double vy102 = 0;

 		if (dim == 2)
 		{
		// conversion to our variables
		double w_theta101 = -sqrt(3.0/2.0)*theta101;
		double w_theta102 = -sqrt(3.0/2.0)*theta102;

		Vector<double> bc_rhs(G20.constants.nBC);
		Vector<double> bc_rhs_manuel(G20.constants.nBC);


		// first we consider the left boundary
		Tensor<1,dim> p;
		Tensor<1,dim> normal_vector;

		p[0] = 0.5;
		p[1] = 0.6;

		normal_vector[0] = -1/sqrt(2);
		normal_vector[1] = 1/sqrt(2);

		FullMatrix<double> projector(2,2);

		projector(0,0) = normal_vector[0];
		projector(0,1) = normal_vector[1];
		projector(1,0) = -normal_vector[1];
		projector(1,1) = normal_vector[0];

		Vector<double> wall_velocity = G20.bcrhs_wall.local_wall_velocity(vx101,vy101,projector);

		bc_rhs_manuel(0) = wall_velocity(0) + 1/sqrt(3*M_PI) * w_theta101-sqrt(2/M_PI) * rho101;
		bc_rhs_manuel(1) = -1/sqrt(M_PI) * wall_velocity(1) ;
		bc_rhs_manuel(2) = -7/sqrt(30*M_PI) * w_theta101-(1/sqrt(5*M_PI)) * rho101;
		bc_rhs_manuel(3) = -sqrt(2/M_PI)/5.0 * w_theta101+2/(5.*sqrt(3*M_PI)) * rho101 ;
		bc_rhs_manuel(4) = 1/(5.*sqrt(2*M_PI)) * w_theta101-1/(5.*sqrt(3*M_PI)) * rho101;


			G20.bcrhs_inflow.BCrhs(p,normal_vector,bc_rhs,101);

			for (unsigned int i = 0 ; i < bc_rhs.size() ; i++)
				EXPECT_NEAR(bc_rhs(i),bc_rhs_manuel(i),1e-5) << "value of i " << i ;

		wall_velocity = G20.bcrhs_wall.local_wall_velocity(vx102,vy102,projector);


			bc_rhs_manuel(0) = wall_velocity(0) + 1/sqrt(3*M_PI) * w_theta102 - sqrt(2/M_PI) * rho102;
			bc_rhs_manuel(1) = -1/sqrt(M_PI) * wall_velocity(1);
			bc_rhs_manuel(2) = -7/sqrt(30*M_PI) * w_theta102-(1/sqrt(5*M_PI)) * rho102;
			bc_rhs_manuel(3) = -sqrt(2/M_PI)/5.0 * w_theta102+2/(5.*sqrt(3*M_PI)) * rho102;
			bc_rhs_manuel(4) = 1/(5.*sqrt(2*M_PI)) * w_theta102-1/(5.*sqrt(3*M_PI)) * rho102;

			G20.bcrhs_inflow.BCrhs(p,normal_vector,bc_rhs,102);

			for (unsigned int i = 0 ; i < bc_rhs.size() ; i++)
				EXPECT_NEAR(bc_rhs(i),bc_rhs_manuel(i),1e-5) << "value of i " << i ; 			
 		}

 		if (dim == 1)
 		{

		Vector<double> bc_rhs(G20.constants.nBC);
		Vector<double> bc_rhs_manuel(G20.constants.nBC);

		double w_theta101 = theta101/sqrt(2);
		double w_theta102 = theta102/sqrt(2);

		// first we consider the left boundary
		Tensor<1,dim> p;
		Tensor<1,dim> normal_vector;

		p[0] = 0.5;

		normal_vector[0] = -1;

		double wall_velocity = -vx101;

		bc_rhs_manuel(0) = wall_velocity + -0.56419 * w_theta101 -0.797885 * rho101;
		bc_rhs_manuel(1) = -0.690988 * w_theta101 + 0.325735 * rho101;
		bc_rhs_manuel(2) = 0.257516 * w_theta101-0.21851 * rho101;



			G20.bcrhs_inflow.BCrhs(p,normal_vector,bc_rhs,101);

			for (unsigned int i = 0 ; i < bc_rhs.size() ; i++)
				EXPECT_NEAR(bc_rhs(i),bc_rhs_manuel(i),1e-5) << "value of i " << i ;

		// wall_velocity = G20.bcrhs_wall.local_wall_velocity(vx102,vy102,projector);


		// 	bc_rhs_manuel(0) = wall_velocity(0) + 1/sqrt(3*M_PI) * w_theta102 - sqrt(2/M_PI) * rho102;
		// 	bc_rhs_manuel(1) = -1/sqrt(M_PI) * wall_velocity(1);
		// 	bc_rhs_manuel(2) = -7/sqrt(30*M_PI) * w_theta102-(1/sqrt(5*M_PI)) * rho102;
		// 	bc_rhs_manuel(3) = -sqrt(2/M_PI)/5.0 * w_theta102+2/(5.*sqrt(3*M_PI)) * rho102;
		// 	bc_rhs_manuel(4) = 1/(5.*sqrt(2*M_PI)) * w_theta102-1/(5.*sqrt(3*M_PI)) * rho102;

		// 	G20.bcrhs_inflow.BCrhs(p,normal_vector,bc_rhs,102);

		// 	for (unsigned int i = 0 ; i < bc_rhs.size() ; i++)
		// 		EXPECT_NEAR(bc_rhs(i),bc_rhs_manuel(i),1e-5) << "value of i " << i ; 				
 		}

	}	

	// TEST(BspecularG20,HandlesBspecularG20)
	// {
	// 	const unsigned int dim = 2;

	// 	std::string folder_name = "../system_matrices/";
	// 	Constants::Base_Constants constants(input_file);
	// 	G20::G20<dim> G20(constants.constants,folder_name);		

	// 	for (unsigned int i = 0 ; i < constants.constants.nBC; i++)	
	// 	{
	// 		for (unsigned int j = 0 ; j < constants.constants.nEqn ; j++)
	// 			std::cout << " " << G20.system_data.B_specular.matrix.coeffRef(i,j);

	// 		printf("\n");
	// 	}	


	// }
}
