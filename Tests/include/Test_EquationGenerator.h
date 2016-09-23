// in the following namespace we test the Equation Generator with the help of A matrix
namespace Test_EquationGenerator
{
	using namespace dealii;

	void develop_A1(Sparse_matrix &A1, const unsigned int nEqn)
	{
		A1.resize(nEqn,nEqn);
		Assert(nEqn == 6, ExcNotImplemented());
		switch(nEqn)
		{
			case 6:
			{
				Assert(A1.rows() !=0 ,ExcNotInitialized());
				A1.coeffRef(0, 1)  = 1.;
				A1.coeffRef(1, 0)  = 1.;
				A1.coeffRef(1, 3)   = 1.;
				A1.coeffRef(2, 4)  = 1.;
				A1.coeffRef(3, 1) =  0.66666666666667;
				A1.coeffRef(4, 2) = 0.50000000000000;
				A1.coeffRef(5, 1) =  -0.33333333333333;

				A1.makeCompressed();
				break;
			}
			default:
			{
				Assert(1 == 0,ExcMessage("Should not have reached here"));
			}
		}

	}

	void develop_A2(Sparse_matrix &A2,const unsigned int nEqn)
	{
		A2.resize(nEqn,nEqn);
		Assert(nEqn == 6, ExcNotImplemented());
		switch(nEqn)
		{
			case 6:
			{
				Assert(A2.rows() !=0,ExcNotInitialized());
				A2.coeffRef(0, 2) =  1.;
				A2.coeffRef(1, 4) =  1.;
				A2.coeffRef(2, 0) =  1.;
				A2.coeffRef(2, 5) =  1.;
				A2.coeffRef(3, 2) =-0.33333333333333;
				A2.coeffRef(4, 1) =  0.50000000000000;
				A2.coeffRef(5, 2) =  0.66666666666667;

				A2.makeCompressed();
				break;
			}
			default:
			{
				Assert(1 == 0,ExcMessage("Should not have reached here"));
			}
		}		
	}

	void develop_B(Sparse_matrix &B,const unsigned int nEqn, const unsigned int nBC)
	{
		B.resize(nBC,nEqn);
		Assert(nEqn == 6, ExcNotImplemented());
		switch(nEqn)
		{
			case 6:
			{
				Assert(B.rows() !=0,ExcNotInitialized());
				B.coeffRef(0, 0) =  -1.0000000000000;
				B.coeffRef(0, 1) =  1.0000000000000;
				B.coeffRef(0, 3) =  -1.0000000000000;
				B.coeffRef(1, 2) =  -1.0000000000000;
				B.coeffRef(1, 4) =   1.0000000000000;

				B.makeCompressed();
				break;
			}
			default:
			{
				Assert(1 == 0,ExcMessage("Should not have reached here"));
			}
		}
	}

	void develop_S_half(Sparse_matrix &S_half,const unsigned int nEqn)
	{
		S_half.resize(nEqn,nEqn);
		Assert(nEqn == 6, ExcNotImplemented());
		switch(nEqn)
		{
			case 6:
			{
				Assert(S_half.rows() != 0,ExcNotInitialized());
				S_half.coeffRef(0, 0) =  1.41421356237310;
				S_half.coeffRef(1, 1) =  1.41421356237310;
				S_half.coeffRef(2, 2) =  1.41421356237310;
				S_half.coeffRef(3, 3) =  1.93185165257814;
				S_half.coeffRef(3, 5) =  0.517638090205042;
				S_half.coeffRef(4, 4) =  2.00000000000000;
				S_half.coeffRef(5, 3) =  0.517638090205042;
				S_half.coeffRef(5, 5) =  1.93185165257814;

				S_half.makeCompressed();
				break;
			}
			default:
			{
				Assert(1 == 0,ExcMessage("Should not have reached here"));
			}
		}
	}

	void develop_S_half_inv(Sparse_matrix &S_half_inv,const unsigned int nEqn)
	{
		S_half_inv.resize(nEqn,nEqn);
		Assert(nEqn == 6, ExcNotImplemented());
		switch(nEqn)
		{
			case 6:
			{
				Assert(S_half_inv.rows() !=0,ExcNotInitialized());
				S_half_inv.coeffRef(0, 0) =  0.70710678118655;
				S_half_inv.coeffRef(1, 1) =  0.70710678118655;
				S_half_inv.coeffRef(2, 2) = 0.70710678118655;
				S_half_inv.coeffRef(3, 3) = 0.55767753582521;
				S_half_inv.coeffRef(3, 5) = -0.14942924536134;
				S_half_inv.coeffRef(4, 4) = 0.50000000000000;
				S_half_inv.coeffRef(5, 3) = -0.14942924536134;
				S_half_inv.coeffRef(5, 5) = 0.55767753582521;

				S_half_inv.makeCompressed();
				break;
			}
			default:
			{
				Assert(1 == 0,ExcMessage("Should not have reached here"));
			}
		}
	}


	void develop_P(Sparse_matrix &P,const unsigned int nEqn,const double tau)
	{
		P.resize(nEqn,nEqn);
		switch(nEqn)
		{
			case 6:
			{
				Assert(P.rows() != 0,ExcNotInitialized());

				// only the first one is the conservation law
				for (unsigned int i = 1 ; i < nEqn ;  i++)
					P.coeffRef(i,i) = 1/tau;

				break;


			}
		}
	}

	TEST(DevelopingSystem,HandlingDevelopingSystems)
	{
		const unsigned int dim = 2;
		ASSERT_EQ(dim,2) << "3D not implemented" << std::endl;

		std::string input_file = "../test_input_files/input1.in";
		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		SystemA::Base_SystemA<dim> systemA(constants.constants,folder_name);

		// no we will manually enter the matrices so as to compare with the read values
		// we will now input all the matrices manually so as to compare with the results
		Sparse_matrix A1;
		Sparse_matrix A2;
		Sparse_matrix B;
		Sparse_matrix S_half;
		Sparse_matrix S_half_inv;
		Sparse_matrix P;

		

		// develop A1
		develop_A1(A1,systemA.constants.nEqn);
		develop_A2(A2,systemA.constants.nEqn);
		develop_B(B,systemA.constants.nEqn,systemA.constants.nBC);
		develop_S_half(S_half,systemA.constants.nEqn);
		develop_S_half_inv(S_half_inv,systemA.constants.nEqn);
		develop_P(P,systemA.constants.nEqn,systemA.constants.tau);

		// now we will compare with the values that have been read
		for (int i = 0 ; i < systemA.constants.nEqn ; i ++)
			for (int j = 0 ; j < systemA.constants.nEqn ; j++)
			{
				EXPECT_LT(fabs(systemA.systemA_data.A[0].matrix.coeffRef(i,j)-A1.coeffRef(i,j)),1e-5);
				EXPECT_LT(fabs(systemA.systemA_data.A[1].matrix.coeffRef(i,j)-A2.coeffRef(i,j)),1e-5);

				EXPECT_LT(fabs(systemA.systemA_data.S_half.matrix.coeffRef(i,j)-S_half.coeffRef(i,j)),1e-5);
				EXPECT_LT(fabs(systemA.systemA_data.S_half_inv.matrix.coeffRef(i,j)-S_half_inv.coeffRef(i,j)),1e-5);
				EXPECT_LT(fabs(systemA.systemA_data.P.matrix.coeffRef(i,j)-P.coeffRef(i,j)),1e-5);
			}


		for (int i = 0 ; i < systemA.constants.nBC ; i ++)
				for (int j = 0 ; j < systemA.constants.nEqn; j ++)
					EXPECT_LT(fabs(systemA.systemA_data.B.matrix.coeffRef(i,j) - B.coeffRef(i,j)),1e-5);
		
	}
}