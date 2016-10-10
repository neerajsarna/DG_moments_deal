// In the following routines we test the boundary handler for characteristic boundary conditions
namespace Test_BoundaryHanlder_Char
{
	using namespace dealii;
	#include "basic_routines.h"

	// we first test the routine for comparing full matrices
	TEST(CompareFloatMat,HandlesCompareFloatMat)
	{
		Full_matrix A;
		Full_matrix B;

		A.resize(2,2);
		B.resize(2,2);

		// we provide unequal matrices and see whether whether the code shouts or not
		A << 1,1,1,1;
		B << 1,1,1,1;

		Compare_Float_Mat(A,B);
	}

	// the systemA does not consist of any velocity therefore the matrix B is not needed to be fixed by 
	// introducing epsilons. So we check this 
	TEST(BMatSystemA,HandlesBMatSystemA)
	{
		const unsigned int dim = 2;
		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		SystemA::SystemA<dim> systemA(constants.constants,folder_name);

		Sparse_matrix B(systemA.constants.nBC,systemA.constants.nEqn);

		B.coeffRef(0, 0) =  -1.0000000000000;
		B.coeffRef(0, 1) =  1.0000000000000;
		B.coeffRef(0, 3) = -1.0000000000000;
		B.coeffRef(1, 2) =  -1.0000000000000;
		B.coeffRef(1, 4) =  1.0000000000000;

		Compare_Float_Mat(B,systemA.system_data.B.matrix);
	}	

	//The matrix X_minus will not necessarily be equal to the one obtained from Mathematica since the eigenvectors
	//are unique only upto a multiplication with a scalar. So we check the matrix Xminus.(B.Xminus).inverse
	TEST(XminusBXminusInv,HandlesXminusBXminusInv)
	{
		const unsigned int dim = 2;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		SystemA::SystemA<dim> systemA(constants.constants,folder_name);

		if (constants.constants.bc_type == characteristic)
		{
			MatrixOpt::Base_MatrixOpt matrix_opt;
			Full_matrix Xminus = matrix_opt.compute_Xminus(systemA.system_data.Ax.matrix,
													  systemA.constants.nBC);

			Full_matrix XminusBXminusInv = Xminus * systemA.B_tilde_inv;

			// result from Mathematica
			Full_matrix XminusBXminusInv_result;
			XminusBXminusInv_result.resize(systemA.constants.nEqn,systemA.constants.nBC);

			XminusBXminusInv_result << -0.338105, 0., 0.436492, 0., 0., -0.585786, -0.225403, 0., 0., 
									0.414214, 0.112702, 0.;

			Compare_Float_Mat(XminusBXminusInv_result,XminusBXminusInv);			
		}

	}

	TEST(BhatCheck,HandlesBhat)
	{
		const unsigned int dim = 2;
		
		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		SystemA::SystemA<dim> systemA(constants.constants,folder_name);

		if (constants.constants.bc_type == characteristic)
		{
			Full_matrix B_hat;
			const unsigned int nEqn = systemA.constants.nEqn;
			B_hat.resize(nEqn,nEqn);

			B_hat << 0.338105, -0.338105, 0., 0.338105, 0., 0., -0.436492, 0.436492, 0., 
			-0.436492, 0., 0., 0., 0., 0.585786, 0., -0.585786, 0., 0.225403, 
			-0.225403, 0., 0.225403, 0., 0., 0., 0., -0.414214, 0., 0.414214, 0., 
			-0.112702, 0.112702, 0., -0.112702, 0., 0.;


			Compare_Float_Mat(B_hat,systemA.B_hat);			
		}
	}

	TEST(BCrhs,HandlesBCrhsSystemA)
	{
		const unsigned int dim = 2;
		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		SystemA::SystemA<dim> systemA(constants.constants,folder_name);

		const double theta1 = 1.0;
		const double theta0 = 2.0;
		const double uW = 0.1;		
		
		if (constants.constants.bc_type = characteristic)
		{
		Tensor<1,dim> p;
		Tensor<1,dim> normal_vector;
		Vector<double> bc_rhs_manuel(systemA.constants.nBC);

		//consider a point on the inner ring
		p[0]= sqrt(pow(0.5,2)/2);
		p[1] = sqrt(pow(0.5,2)/2);

		normal_vector[0]= 1/sqrt(2.0);
		normal_vector[1] = 1/sqrt(2.0);

		bc_rhs_manuel(0) = -systemA.constants.chi * theta0;
		bc_rhs_manuel(1) = -uW * normal_vector[1];

		Vector<double> bc_rhs(2);
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


		bc_rhs_manuel(0) = -systemA.constants.chi * theta1;

		systemA.build_BCrhs(p,normal_vector,bc_rhs);

		for (unsigned int i = 0 ; i < 2 ; i++)
			EXPECT_NEAR(bc_rhs_manuel(i),bc_rhs(i),1e-5);
	
		}
		
	}
}