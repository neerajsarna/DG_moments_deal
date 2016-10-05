//Testing of the basic matrix routines
namespace TestMatrixOpt
{
	using namespace dealii;

	// test the ordering function
	TEST(TestingOrdering,HandlesOrdering)
	{
		VectorXd array;
		const unsigned int num_entries = 5;
		std::vector<int> result(num_entries);
		array.resize(num_entries);

		array << 7.8, 2.5, 9.8,100.2 , -0.5;

		// location of the first lowest number
		result[0] = 4;

		// location of the second lowest number
		result[1] = 1;

		result[2] = 0;

		result[3] = 2;

		result[4] = 3;

		MatrixOpt::Base_MatrixOpt::Ordering ordering(array);

		for (unsigned int i = 0 ; i < result.size() ; i ++)
			EXPECT_EQ(result[i],ordering.index[i]);

	}

	// we will first consider a small matrix
	TEST(TestXminusDummyMat,HandlingXminus)
	{
		// we will first consider a simple matrix
		Sparse_matrix A;
		A.resize(5,5);

		A.coeffRef(0,0) = -1.0;
		A.coeffRef(1,1) = -1.0;

		// number of negative eigenvalues in A
		const unsigned int num_neg = 2;

		MatrixOpt::Base_MatrixOpt matrix_opt;
		Full_matrix Xminus1_result;

		Xminus1_result.resize(5,num_neg);
		Xminus1_result << 1,0,0,1,0,0,0,0,0,0;

		Full_matrix Xminus1 = matrix_opt.compute_Xminus(A,num_neg);

		for(unsigned int i = 0 ; i < Xminus1.rows() ; i++)
			for (unsigned int j = 0 ; j < Xminus1.cols() ; j++)
				EXPECT_FLOAT_EQ(Xminus1(i,j),Xminus1_result(i,j)) << "Incorrect negative eigenvectors";


		// now we compute the negative eigenvalues
		VectorXd neg_vals = matrix_opt.compute_Lambda_minus(A,num_neg);
		for (unsigned int i = 0 ; i < num_neg ; i ++)
			EXPECT_FLOAT_EQ(neg_vals(i),-1.0) << "Incorrect negative eigenvalues " ;

	}

	// we will now consider the matrix from system A
	TEST(TestXminusSystemA,HandlingXminusSystemA)
	{
		
		std::string folder_name = "../system_matrices/";

		Constants::Base_Constants constants(input_file);
		const unsigned int dim = 2;
		SystemA::SystemA<dim> systemA(constants.constants,folder_name);

		MatrixOpt::Base_MatrixOpt matrix_opt;
		Full_matrix Xminus = matrix_opt.compute_Xminus(systemA.system_data.A[0].matrix,2);
		VectorXd neg_vals = matrix_opt.compute_Lambda_minus(systemA.system_data.A[0].matrix,2);

		// now we check the computations
		Full_matrix A_Xminus = systemA.system_data.A[0].matrix * Xminus;
		Full_matrix Xminus_negvals;

		Xminus_negvals.resize(A_Xminus.rows(),A_Xminus.cols());

		// we now multiply Xminus from the right by the negative eigenvalues
		for (unsigned int i = 0 ; i < Xminus_negvals.cols(); i ++)
				Xminus_negvals.col(i) = Xminus.col(i)  * neg_vals(i);

		// we now compute the difference
		Full_matrix Diff = A_Xminus-Xminus_negvals;

		for (unsigned int i = 0 ; i < Diff.rows(); i ++)
			for (unsigned int j = 0 ; j < Diff.cols() ; j ++)
				EXPECT_NEAR(A_Xminus(i,j),Xminus_negvals(i,j),1e-10);

	}

	TEST(TestLambdaXcomputationSystemA,HandlingEigenDecomp)
	{
		std::string folder_name = "../system_matrices/";

		Constants::Base_Constants constants(input_file);
		const unsigned int dim = 2;
		SystemA::SystemA<dim> systemA(constants.constants,folder_name);

		MatrixOpt::Base_MatrixOpt matrix_opt;
		Sparse_matrix A = systemA.system_data.A[0].matrix;

		// first we compute the eigenvectors
		Full_matrix vecs = matrix_opt.compute_X(A);

		// next we compute the eigenvalues
		VectorXd vals = matrix_opt.compute_Lambda(A);

		// multiply the matrix from the right by the eigenvectors
		Full_matrix A_vecs = A * vecs;

		// we multiply the vecs by the eigenvalues from the right
		for (unsigned int i = 0 ; i < vals.size() ; i ++)
			vecs.col(i) *= vals(i);

		// compute the difference in the computation
		Full_matrix Diff = A_vecs-vecs;

		for (unsigned int i = 0 ; i < Diff.rows() ; i++)
			for (unsigned int j = 0 ; j < Diff.cols(); j ++)
				EXPECT_NEAR(Diff(i,j),0,1e-7);

	}

	TEST(AmodComputationSystemA,HandlingAmod)
	{
		std::string folder_name = "../system_matrices/";

		Constants::Base_Constants constants(input_file);
		const unsigned int dim = 2;
		SystemA::SystemA<dim> systemA(constants.constants,folder_name);	

		MatrixOpt::Base_MatrixOpt matrix_opt;
		Sparse_matrix A = systemA.system_data.Ax.matrix;
		Full_matrix Amod = matrix_opt.compute_Amod(A);

		Full_matrix Amod_result;
		Amod_result.resize(A.rows(),A.cols());

		// value taken from mathematica
		Amod_result << 0.774597, 0, 0, 0.774597, 0, 0, 0, 1.29099, 0, 0, 0, 0, 0, 0, 
						0.707107, 0, 0, 0, 0.516398, 0, 0, 0.516398, 0, 0, 0, 0, 0, 0, 
						0.707107, 0, -0.258199, 0, 0, -0.258199, 0, 0;

		Compare_Float_Mat(Amod_result,Amod);
	}

	TEST(AminusComputationSystemA,HandlingAminus)
	{
		
		std::string folder_name = "../system_matrices/";

		Constants::Base_Constants constants(input_file);
		const unsigned int dim = 2;
		SystemA::SystemA<dim> systemA(constants.constants,folder_name);	

		MatrixOpt::Base_MatrixOpt matrix_opt;
		Sparse_matrix A = systemA.system_data.Ax.matrix;
		Full_matrix Aminus = matrix_opt.compute_Aminus(A);

		Full_matrix Aminus_result;
		Aminus_result.resize(A.rows(),A.cols());

		Aminus_result << 0.774597, -1., 0, 0.774597, 0, 0, -1., 1.29099, 0, -1., 0, 0, 0, 0, 
						0.707107, 0, -1., 0, 0.516398, -0.666667, 0, 0.516398, 0, 0, 0, 0, 
						-0.5, 0, 0.707107, 0, -0.258199, 0.333333, 0, -0.258199, 0, 0;

		Compare_Float_Mat(Aminus,Aminus_result);
	}
}