namespace Test_BoundaryHandler_Odd
{
	using namespace dealii;

	Full_matrix develop_BC(const unsigned int nEqn)
	{
		Full_matrix BC;
		BC.resize(nEqn,nEqn);

		BC << 1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.00001,0.,0.,
   	-8.164965809277262e-6,0.000014142135623730951,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
   0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,
   0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
   0.5641895835477563,0.,0.,0.,0.,0.,-0.1784124116152771,0.,1.1845513696160073,0.,0.,
   0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
   0.8240516309828045,-0.3568248232305542,0.,0.,0.,0.,0.,0.,0.,0.,-0.737054185538849,
   0.9536544540177923,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
   0.,0.10638460810704872,0.6449224123464928,0.,0.,0.,0.,0.,0.,0.,0.,
   -0.09515328619481445,-0.24623252122982905,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,
   0.,0.,0.,0.,0.,0.,0.,0.,-0.05319230405352436,-0.09213177319235613,0.,
   0.46065886596178063,0.,0.,0.,0.,0.,0.,0.04757664309740722,0.,0.,
   -0.24623252122982908,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,
   0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,
   0.,0.,0.,0.,0.15078600877302686,0.,0.,0.,0.,0.,0.5245099497097857,0.,
   -0.3165846701526373,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
   0.,1.;

		return(BC);
	}

	TEST(BCOddG26,HandlingBCOddG26)
	{
		const unsigned int dim = 2;
		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);

		G26::G26<dim> G26(constants.constants,folder_name);

		if (constants.constants.bc_type == odd)
		{
			Full_matrix BC = develop_BC(G26.constants.nEqn);
			Compare_Float_Mat(BC,G26.BC);			
		}

	}

	TEST(BCrhsOddG26,HandlingBCrhsOddG26)
	{
		const unsigned int dim = 2;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		G26::G26<dim> G26(constants.constants,folder_name);		

		const double  theta0 = 1.0;
		const double  theta1 = 2.0;
		const double  theta2 = 3.0;
		const double  theta3 = 4.0;

 
		const double  vn0 = 1.0;
		const double  vn1 = 2.0;
		const double  vn2 = 3.0;
		const double  vn3 = 4.0;

 
		const double  vt0 = 1.0;
		const double  vt1 = 2.0;
		const double  vt2 = 3.0;
		const double  vt3 = 4.0;
		
		if (constants.constants.bc_type == odd)
		{
			// for this particular system, p and normal_vector do not matter. Everything regulated by b_id
			// consider the wall with b_id == 0
			Tensor<1,dim> p;
			Tensor<1,dim> normal_vector;
			unsigned int b_id = 0;
			Vector<double> bc_rhs_manuel_char(G26.constants.nBC);
			Vector<double> bc_rhs_manuel(G26.constants.nEqn);
			Vector<double> bc_rhs(G26.constants.nEqn);
			bc_rhs_manuel = 0;
			std::vector<double> odd_coeff(G26.system_data.odd_ID.size());

			for (int i = 0 ; i < G26.system_data.odd_ID.size() ; i++)
				odd_coeff[i] = G26.system_data.B.matrix.coeffRef(i,G26.system_data.odd_ID.coeffRef(i,0));

			// we add just for precaution
			bc_rhs_manuel_char(0) += vn0/odd_coeff[0];
			bc_rhs_manuel_char(1) += -0.5 * vt0/odd_coeff[1];
			bc_rhs_manuel_char(2) += -2*sqrt(2/15.0) * (-sqrt(3/2.0)) * theta0/odd_coeff[2];
			bc_rhs_manuel_char(3) += -sqrt(2)/15.0 * (-sqrt(3/2.0)) * theta0/odd_coeff[3];
			bc_rhs_manuel_char(4) += 1/(15.0 * sqrt(2)) * (-sqrt(3.0/2.0)) * theta0/odd_coeff[4];
			bc_rhs_manuel_char(5) += -1/(2*sqrt(14)) * vt0/odd_coeff[5];

			bc_rhs_manuel = 0;
			bc_rhs = 0;

			for (int i = 0 ; i < G26.system_data.odd_ID.rows() ; i++)
				bc_rhs_manuel(G26.system_data.odd_ID(i,0)) = bc_rhs_manuel_char(i);
			


			G26.build_BCrhs(p,normal_vector,bc_rhs,b_id);

			for (int i = 0 ; i < G26.constants.nEqn ; i++)
				EXPECT_NEAR(bc_rhs(i),bc_rhs_manuel(i),1e-5) << "Failed at " << i;	
		}

	}


}