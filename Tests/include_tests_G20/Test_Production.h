namespace TestProduction
{

	// production term for the Maxwell Molecules
	MatrixXd build_P_MM(const unsigned int nEqn,const double tau)
	{
		MatrixXd P(nEqn,nEqn);

		P <<0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.,0,0,0,0,0,0,0,0,0,0,0,0,0,1.0000000000000002,0,0,0,
   0,0,0,0,0,0,0,0,0,0,1.,0,0,0,0,0,0,0,0,0,0,0,0,0,0.6666666666666666,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0.6666666666666666,0,0,0,0,0,0,0,0,0,0,0,0,0,1.4999999999999998,0,0,0,0,0,
   0,0,0,0,0,0,0,0,1.4999999999999996,0,0,0,0,0,0,0,0,0,0,0,0,0,1.5000000000000002,0,0,
   0,0,0,0,0,0,0,0,0,0,0,1.5000000000000004;

		return(P/tau);
	}

	TEST(ProductionMM,HandlesProductionTermMM)
	{
		const unsigned int dim = 2;
		ASSERT_EQ(dim,2) << "3D not implemented" << std::endl;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		G20::G20<dim> G20(constants.constants,folder_name);

		MatrixXd P_manuel = build_P_MM(constants.constants.nEqn,constants.constants.tau);

		for (unsigned int i = 0 ; i < constants.constants.nEqn ; i++)
				EXPECT_NEAR(P_manuel(i,i),G20.system_data.P.matrix.coeffRef(i,i),1e-5);
	}
}