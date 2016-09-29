namespace Test_BoundaryHandler_Odd
{
	using namespace dealii;

	Sparse_matrix develop_BC_systemA(const unsigned int nEqn)
	{
		Sparse_matrix BC;
		BC.resize(nEqn,nEqn);

		BC.coeffRef(0, 0) = 1.0;
		BC.coeffRef(1, 0) = 1.0;
		BC.coeffRef(1, 3) = 1.0;
		BC.coeffRef(2, 2) = 1.0;
		BC.coeffRef(3, 3) = 1.0;
		BC.coeffRef(4, 2) = 1.0;
		BC.coeffRef(5, 5) = 1.0;

		return(BC);
	}

	TEST(BCOddSystemA,HandlingBCOddSystemA)
	{
		const unsigned int dim = 2;
		std::string input_file = "../test_input_files/input1.in";
		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		SystemA::SystemA<dim> systemA(constants.constants,folder_name);

		Sparse_matrix BC = develop_BC_systemA(systemA.constants.nEqn);

		Compare_Float_Mat(BC,systemA.BC);
	}

	TEST(BCrhsOddSystemA,HandlingBCrhsOddSystemA)
	{
		const unsigned int dim = 2;
		std::string input_file = "../test_input_files/input1.in";
		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		SystemA::SystemA<dim> systemA(constants.constants,folder_name);		

		const double theta1 = 1.0;
		const double theta0 = 2.0;
		const double uW = 0.1;		
		
		Tensor<1,dim> p;
		Tensor<1,dim> normal_vector;
		Vector<double> bc_rhs_manuel(systemA.constants.nEqn);

		//consider a point on the inner ring
		p[0]= sqrt(pow(0.5,2)/2);
		p[1] = sqrt(pow(0.5,2)/2);

		normal_vector[0]= 1/sqrt(2.0);
		normal_vector[1] = 1/sqrt(2.0);

		// boundary condition for qx and Rxy
		bc_rhs_manuel(1) = -systemA.constants.chi * theta0;
		bc_rhs_manuel(4) = -uW * normal_vector[1];

		Vector<double> bc_rhs(systemA.constants.nEqn);
		systemA.build_BCrhs(p,normal_vector,bc_rhs);

		for (unsigned int i = 0 ; i < 2 ; i++)
			EXPECT_NEAR(bc_rhs_manuel(i),bc_rhs(i),1e-5);

		// consider a point on the outer ring now 
		bc_rhs = 0;
		bc_rhs_manuel = 0;

		p[0]= 1.0;
		p[1] = 0.0;

		normal_vector[0]= 1/sqrt(2.0);
		normal_vector[1] = 1/sqrt(2.0);


		// boundary condition for qx
		bc_rhs_manuel(1) = -systemA.constants.chi * theta1;

		systemA.build_BCrhs(p,normal_vector,bc_rhs);

		for (int i = 0 ; i < systemA.constants.nEqn ; i++)
			EXPECT_NEAR(bc_rhs_manuel(i),bc_rhs(i),1e-5);
	}

}