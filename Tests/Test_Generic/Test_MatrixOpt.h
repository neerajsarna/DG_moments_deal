//Testing of the basic matrix routines
namespace TestMatrixOpt
{
	using namespace dealii;

	// // test the ordering function
	// TEST(TestingOrdering,HandlesOrdering)
	// {
	// 	VectorXd array;
	// 	const unsigned int num_entries = 5;
	// 	std::vector<int> result(num_entries);
	// 	array.resize(num_entries);

	// 	array << 7.8, 2.5, 9.8,100.2 , -0.5;

	// 	// location of the first lowest number
	// 	result[0] = 4;

	// 	// location of the second lowest number
	// 	result[1] = 1;

	// 	result[2] = 0;

	// 	result[3] = 2;

	// 	result[4] = 3;

	// 	MatrixOpt::Base_MatrixOpt::Ordering ordering(array);

	// 	for (unsigned int i = 0 ; i < result.size() ; i ++)
	// 		EXPECT_EQ(result[i],ordering.index[i]);

	// }

	// // we will first consider a small matrix
	// TEST(TestXminusDummyMat,HandlingXminus)
	// {
	// 	// we will first consider a simple matrix
	// 	Sparse_matrix A;
	// 	A.resize(5,5);

	// 	A.coeffRef(0,0) = -1.0;
	// 	A.coeffRef(1,1) = -1.0;

	// 	// number of negative eigenvalues in A
	// 	const unsigned int num_neg = 2;

	// 	MatrixOpt::Base_MatrixOpt matrix_opt;
	// 	Full_matrix Xminus1_result;

	// 	Xminus1_result.resize(5,num_neg);
	// 	Xminus1_result << 1,0,0,1,0,0,0,0,0,0;

	// 	Full_matrix Xminus1 = matrix_opt.compute_Xminus(A,num_neg);

	// 	for(unsigned int i = 0 ; i < Xminus1.rows() ; i++)
	// 		for (unsigned int j = 0 ; j < Xminus1.cols() ; j++)
	// 			EXPECT_FLOAT_EQ(Xminus1(i,j),Xminus1_result(i,j)) << "Incorrect negative eigenvectors";


	// 	// now we compute the negative eigenvalues
	// 	VectorXd neg_vals = matrix_opt.compute_Lambda_minus(A,num_neg);
	// 	for (unsigned int i = 0 ; i < num_neg ; i ++)
	// 		EXPECT_FLOAT_EQ(neg_vals(i),-1.0) << "Incorrect negative eigenvalues " ;

	// }

	// // we will now consider the matrix from system A
	// TEST(TestXminusG26,HandlingXminusg26)
	// {
		
	// 	std::string folder_name = "../system_matrices/";

	// 	Constants::Base_Constants constants(input_file);
	// 	const unsigned int dim = 2;
	// 	G26::G26<dim> G26(constants.constants,folder_name);

	// 	MatrixOpt::Base_MatrixOpt matrix_opt;
	// 	Full_matrix Xminus = matrix_opt.compute_Xminus(G26.system_data.A[0].matrix,2);
	// 	VectorXd neg_vals = matrix_opt.compute_Lambda_minus(G26.system_data.A[0].matrix,2);

	// 	// now we check the computations
	// 	Full_matrix A_Xminus = G26.system_data.A[0].matrix * Xminus;
	// 	Full_matrix Xminus_negvals;

	// 	Xminus_negvals.resize(A_Xminus.rows(),A_Xminus.cols());

	// 	// we now multiply Xminus from the right by the negative eigenvalues
	// 	for (unsigned int i = 0 ; i < Xminus_negvals.cols(); i ++)
	// 			Xminus_negvals.col(i) = Xminus.col(i)  * neg_vals(i);

	// 	// we now compute the difference
	// 	Full_matrix Diff = A_Xminus-Xminus_negvals;

	// 	for (unsigned int i = 0 ; i < Diff.rows(); i ++)
	// 		for (unsigned int j = 0 ; j < Diff.cols() ; j ++)
	// 			EXPECT_NEAR(A_Xminus(i,j),Xminus_negvals(i,j),1e-10);

	// }

	// TEST(TestLambdaXcomputationG26,HandlingEigenDecomp)
	// {
	// 	std::string folder_name = "../system_matrices/";

	// 	Constants::Base_Constants constants(input_file);
	// 	const unsigned int dim = 2;
	// 	G26::G26<dim> G26(constants.constants,folder_name);

	// 	MatrixOpt::Base_MatrixOpt matrix_opt;
	// 	Sparse_matrix A = G26.system_data.A[0].matrix;

	// 	// first we compute the eigenvectors
	// 	Full_matrix vecs = matrix_opt.compute_X(A);

	// 	// next we compute the eigenvalues
	// 	VectorXd vals = matrix_opt.compute_Lambda(A);

	// 	// multiply the matrix from the right by the eigenvectors
	// 	Full_matrix A_vecs = A * vecs;

	// 	// we multiply the vecs by the eigenvalues from the right
	// 	for (unsigned int i = 0 ; i < vals.size() ; i ++)
	// 		vecs.col(i) *= vals(i);


	// 	Compare_Float_Mat(A_vecs,vecs);

	// }

	// TEST(AmodComputationG26,HandlingAmod)
	// {
	// 	std::string folder_name = "../system_matrices/";

	// 	Constants::Base_Constants constants(input_file);
	// 	const unsigned int dim = 2;
	// 	G26::G26<dim> G26(constants.constants,folder_name);	

	// 	MatrixOpt::Base_MatrixOpt matrix_opt;
	// 	Sparse_matrix A = G26.system_data.Ax.matrix;
	// 	Full_matrix Amod = matrix_opt.compute_Amod(A);

	// 	Full_matrix Amod_result;
	// 	Amod_result.resize(A.rows(),A.cols());

	// 	// value taken from mathematica
	// 	Amod_result << 0.6983500078748186,0,0,-0.4489457046876043,0.5838140478417588,0,0,0,0,0,0,0,0,
 //   -0.10845349090738622,0.3281012787590295,0,0,0,1.6153378036822688,0,0,0,0,0,
 //   -0.4364622479906845,0,0.7074328343981866,0,0,0,0,0,0,0,0,0,0.7183599712790059,0,0,0,
 //   0,0,-0.07458187476572165,0,1.3393979463089787,0,0,0,0,0,0,-0.4489457046876042,0,0,
 //   1.1002483706152166,-0.5675522528130044,0,0,0,0,0,0,0,0,-0.6562284734616644,
 //   0.6432593894035369,0,0,0.3892093652278391,0,0,-0.37836816854200267,
 //   1.670740290296107,0,0,0,0,0,0,0,0,0.054184446249400345,-0.6549198548712548,0,0,0,0,
 //   0,0,0,1.6404207654943488,0,0,0,0,0,0,0,0,0,-0.5558954147453772,0,
 //   -0.1946046826139195,0,0,0.18918408427100128,-0.39441159330395537,0,
 //   0.8819171036881968,0,0,0,0,0,0,-0.027092223124700186,0.0917576670401116,0,
 //   -0.4714045207910316,0,-0.43646224799068517,0,0,0,0,0,2.2291347948801463,0,
 //   -0.7753366564492661,0,0,0,0,0,0,0,0,0,-0.07458187476572192,0,0,0,0,0,
 //   1.1329317688178246,0,-1.3841556331746234,0,0,0,0,0,0,0,0.28297313375927463,0,0,0,0,
 //   0,-0.31013466257970634,0,1.3688105989793777,0,0,0,0,0,0,0,0,0,0.3571727856823945,0,
 //   0,0,0,0,-0.3691081688465664,0,1.0333649937764755,0,0,0,0,0,0,0,-0.14148656687963732,
 //   0,0,0,0,0,0.15506733128985317,0,-0.11745858997584796,0,1.1338934190276817,0,0,0,0,0,
 //   0,0,-0.26787958926179584,0,0,0,0,0,0.2768311266349249,0,-0.7750237453323567,0,0,0,0,
 //   0,0,-0.10845349090738618,0,0,-0.6562284734616646,0.08127666937410112,0,0,0,0,0,0,0,
 //   0,0.666151821386022,-0.8149599841120081,0,0,0.218734185839353,0,0,
 //   0.42883959293569096,-0.6549198548712545,0,0,0,0,0,0,0,0,-0.5433066560746718,
 //   1.077792707369628,0,0,0,0,0,0,0,-0.555895414745377,0,0,0,0,0,0,0,0,0,
 //   1.2442359683789577,0,-0.10936709291967656,0,0,-0.21441979646784548,
 //   0.09175766704011143,0,-0.4714045207910316,0,0,0,0,0,0,0.2716533280373359,
 //   -0.4129081960150717,0,0.2519763153394847;

	// 	Compare_Float_Mat(Amod_result,Amod);
	// }

	// TEST(AminusComputationG26,HandlingAminus)
	// {
		
	// 	std::string folder_name = "../system_matrices/";

	// 	Constants::Base_Constants constants(input_file);
	// 	const unsigned int dim = 2;
	// 	G26::G26<dim> G26(constants.constants,folder_name);	

	// 	MatrixOpt::Base_MatrixOpt matrix_opt;
	// 	Sparse_matrix A = G26.system_data.Ax.matrix;
	// 	Full_matrix Aminus = matrix_opt.compute_Aminus(A);

	// 	Full_matrix Aminus_result;
	// 	Aminus_result.resize(A.rows(),A.cols());

	// 	Aminus_result << 0.6983500078748186,-1.,0,-0.4489457046876043,0.5838140478417588,0,0,0,0,0,0,0,0,
 //   -0.10845349090738622,0.3281012787590295,0,0,-1.,1.6153378036822688,0,
 //   0.816496580927726,-1.4142135623730951,0,0,-0.4364622479906845,0,0.7074328343981866,
 //   0,0,0,0,0,0,0,0,0,0.7183599712790059,0,0,-1.4142135623730951,0,0,
 //   -0.07458187476572165,0,1.3393979463089787,0,0,0,0,0,0,-0.4489457046876042,
 //   0.816496580927726,0,1.1002483706152166,-0.5675522528130044,0,0,-1.2909944487358056,
 //   0,0,0,0,0,-0.6562284734616644,0.6432593894035369,0,0,0.3892093652278391,
 //   -0.9428090415820634,0,-0.37836816854200267,1.670740290296107,0,0,0.5962847939999438,
 //   0,-1.7320508075688772,0,0,0,0.054184446249400345,-0.6549198548712548,0,0,0,0,
 //   -0.7071067811865475,0,0,1.6404207654943488,0,0,0.4472135954999579,0,
 //   -1.7320508075688772,0,0,0,0,-0.5558954147453772,0,-0.1946046826139195,
 //   0.4714045207910317,0,0.18918408427100128,-0.39441159330395537,0,0.8819171036881968,
 //   -0.2981423969999719,0,0,0,-1.7320508075688772,0,-0.027092223124700186,
 //   0.0917576670401116,0,-0.4714045207910316,0,-0.43646224799068517,0,
 //   -1.2909944487358056,0.8944271909999159,0,0,2.2291347948801463,0,-0.7753366564492661,
 //   0,0,0,1.1547005383792517,-1.6733200530681511,0,0,0,0,-0.07458187476572192,0,0,
 //   0.8944271909999159,0,0,1.1329317688178246,0,-1.3841556331746234,0,0,0,0,
 //   -1.6733200530681511,0,0,0.28297313375927463,0,0,-1.0392304845413263,0,0,
 //   -0.31013466257970634,0,1.3688105989793777,0,0,0,0,0.5554920598635308,0,0,0,0,
 //   0.3571727856823945,0,0,-0.9237604307034014,0,0,-0.3691081688465664,0,
 //   1.0333649937764755,0,0,0,0,0.49377071987869414,0,0,-0.14148656687963732,0,0,
 //   0.23094010767585035,0,-0.5773502691896258,0.15506733128985317,0,
 //   -0.11745858997584796,0,1.1338934190276817,0,0,-0.12344267996967354,0,
 //   0.3086066999241838,0,0,-0.26787958926179584,0,0,0.6928203230275509,0,0,
 //   0.2768311266349249,0,-0.7750237453323567,0,0,0,0,-0.3703280399090206,0,
 //   -0.10845349090738618,0,0,-0.6562284734616646,0.08127666937410112,0,0,
 //   1.1547005383792517,0,0,0,0,0,0.666151821386022,-0.8149599841120081,0,0,
 //   0.218734185839353,0,0,0.42883959293569096,-0.6549198548712545,0,0,
 //   -1.115546702045434,0,0.9258200997725514,0,0,0,-0.5433066560746718,1.077792707369628,
 //   0,0,0,0,0,0,0,-0.555895414745377,0,0,-0.8366600265340756,0,0.9258200997725514,0,0,0,
 //   0,1.2442359683789577,0,-0.10936709291967656,0,0,-0.21441979646784548,
 //   0.09175766704011143,0,-0.4714045207910316,0.557773351022717,0,0,0,
 //   0.9258200997725514,0,0.2716533280373359,-0.4129081960150717,0,0.2519763153394847;

	// 	Compare_Float_Mat(Aminus,Aminus_result);
	// }

	TEST(SquareOuterProductAB,HandlesSquareOuterProductAB)
	{
		Full_matrix A(2,2);
		FullMatrix<double> B(2,2);
		FullMatrix<double> result(4,4);
		FullMatrix<double> result_manuel(4,4);

		for (unsigned int i = 0 ; i < 2 ; i ++)
			for (unsigned int j = 0 ; j < 2 ; j ++)
			{
				A(i,j) = pow(i + j,3);
				B(i,j) = i + j + 2;
			}

		MatrixOpt::Base_MatrixOpt matrix_opt;
		result = matrix_opt.compute_A_outer_B(A,B);

		result_manuel = 0;

		result_manuel.fill(matrix_opt.multiply_scalar(A(0,1),B),0,2);
		result_manuel.fill(matrix_opt.multiply_scalar(A(1,0),B),2,0);
		result_manuel.fill(matrix_opt.multiply_scalar(A(1,1),B),2,2);

		for (unsigned int i = 0 ; i < 4 ; i ++)
			for (unsigned int j = 0 ; j < 4 ; j ++)
				EXPECT_NEAR(result(i,j),result_manuel(i,j),1e-5);

	}


	TEST(RectangeOuterProductSparseAB,HandlesRectangleOuterProductSparseAB)
	{
		Sparse_matrix A(2,2);
		FullMatrix<double> B(2,3);
		FullMatrix<double> result(4,6);
		FullMatrix<double> result_manuel(4,6);

		for (unsigned int i = 0 ; i < 2 ; i ++)
			for (unsigned int j = 0 ; j < 3 ; j ++)
			{
				if (j < 2)
					A.coeffRef(i,j) = pow(i + j , 3);

				B(i,j) = i + j + 2;
			}

		A.makeCompressed();
		MatrixOpt::Base_MatrixOpt matrix_opt;
		result = matrix_opt.compute_A_outer_B(A,B);

		// manual result
		result_manuel.fill(matrix_opt.multiply_scalar(A.coeffRef(0,0),B),0,0);
		result_manuel.fill(matrix_opt.multiply_scalar(A.coeffRef(0,1),B),0,3);
		result_manuel.fill(matrix_opt.multiply_scalar(A.coeffRef(1,0),B),2,0);
		result_manuel.fill(matrix_opt.multiply_scalar(A.coeffRef(1,1),B),2,3);

		for (unsigned int i = 0 ; i < result_manuel.m() ; i++)
			for (unsigned int j = 0 ; j < result_manuel.n() ; j++)
				EXPECT_NEAR(result(i,j),result_manuel(i,j),1e-5);
	}

}