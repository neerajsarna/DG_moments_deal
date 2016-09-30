// in this class we define certain basi matrix operatorions
namespace MatrixOpt
{
	using namespace dealii;
	using namespace std;

	class
	Base_MatrixOpt
	{
		public:
			Base_MatrixOpt(){;};

			struct Ordering
			{
				Ordering( VectorXd &Array );
				bool operator() (int i,int j);    
				VectorXi index;  
				VectorXd array;
			};

			// given as matrix, computes the eigenvectors corresonding to the negative eigenvalues
			// A is the sparse matrix and num_neg is the number of negative eigenvalues
			Full_matrix compute_Xminus(const Sparse_matrix &A,const unsigned int num_neg);

			// same as above but for full matrices
			Full_matrix compute_Xminus(const Full_matrix &A,const unsigned int num_neg);

			// Return a vector containing the negative eigenvalues of the input matrix
			VectorXd compute_Lambda_minus(Sparse_matrix &A,const unsigned int num_neg);

			// same as above but for full matrices
			VectorXd compute_Lambda_minus(Full_matrix &A,const unsigned int num_neg);


			// computes all the eigenvectors of the matrix A
			Full_matrix compute_X(Sparse_matrix &A);
			Full_matrix compute_X(Full_matrix &A);


			// computes all the eigenvalues of the matrix A
			VectorXd compute_Lambda(Sparse_matrix &A);
			VectorXd compute_Lambda(Full_matrix &A);


			// computes the modulus of a matrix
			Full_matrix compute_Amod(Sparse_matrix &A);
			Full_matrix compute_Amod(Full_matrix &A);


			// computes the flux matrix for a sparse matrix
			Full_matrix compute_Aminus(Sparse_matrix &A);

			// computes the difference between a full matrix and a sparse matrix

			// prints a full matrix
			void print_eigen_mat(Full_matrix &full_matrix,
								std::string matrix_name);

			// same as above but for sparse matrix
			void print_eigen_mat(Sparse_matrix &full_matrix,
								std::string matrix_name);
	};

	// constructor of the structure
	Base_MatrixOpt::Ordering::Ordering(VectorXd &Array)
	{
		int n = Array.size();
		array = Array;

		vector<int> sorted(n);
		for( int i = 0; i<n; i++ ) sorted[i] = i;
			sort(sorted.begin(),sorted.end(),*this);
		index = Map<VectorXi>(sorted.data(),n);
	}

	bool Base_MatrixOpt::Ordering::operator()(int i,int j)
	{
		  return( array[i]<array[j] );
	}

	Full_matrix Base_MatrixOpt::compute_Xminus(const Sparse_matrix &A,const unsigned int num_neg)
	{
		Full_matrix Xminus;
		Xminus.resize(A.rows(),num_neg);

		EigenSolver<MatrixXd> ES(A);
		MatrixXd vecs = ES.pseudoEigenvectors();
		VectorXd vals = ES.pseudoEigenvalueMatrix().diagonal();

		Ordering order(vals);
		

		Assert(Xminus.cols() != 0 && Xminus.rows() !=0 ,ExcNotInitialized());
		for (unsigned int i = 0 ; i < num_neg; i ++)
			Xminus.col(i) = vecs.col(order.index(i));

		return (Xminus);
	}

	Full_matrix Base_MatrixOpt::compute_Xminus(const Full_matrix &A,const unsigned int num_neg)
	{
		Full_matrix Xminus;
		Xminus.resize(A.rows(),num_neg);

		EigenSolver<MatrixXd> ES(A);
		MatrixXd vecs = ES.pseudoEigenvectors();
		VectorXd vals = ES.pseudoEigenvalueMatrix().diagonal();

		Ordering order(vals);
		

		Assert(Xminus.cols() != 0 && Xminus.rows() !=0 ,ExcNotInitialized());
		for (unsigned int i = 0 ; i < num_neg; i ++)
			Xminus.col(i) = vecs.col(order.index(i));

		return (Xminus);
	}

	VectorXd Base_MatrixOpt::compute_Lambda_minus(Sparse_matrix &A,const unsigned int num_neg)
	{
		// Create the eigensolver object for A
		EigenSolver<MatrixXd> ES(A);

		// compute the eigenvectors of the matrix A
		MatrixXd vecs = ES.pseudoEigenvectors();

		// compute the eigen values and the eigenvectors of the matrix
		VectorXd vals = ES.pseudoEigenvalueMatrix().diagonal();
		VectorXd neg_vals;

		neg_vals.resize(num_neg);
		Ordering order(vals);

		for (unsigned int i = 0 ; i < num_neg ; i ++)
			neg_vals(i) = vals(order.index[i]);

		return(neg_vals);
	}


	VectorXd Base_MatrixOpt::compute_Lambda_minus(Full_matrix &A,const unsigned int num_neg)
	{
		// Create the eigensolver object for A
		EigenSolver<MatrixXd> ES(A);

		// compute the eigenvectors of the matrix A
		MatrixXd vecs = ES.pseudoEigenvectors();

		// compute the eigen values and the eigenvectors of the matrix
		VectorXd vals = ES.pseudoEigenvalueMatrix().diagonal();
		VectorXd neg_vals;

		neg_vals.resize(num_neg);
		Ordering order(vals);

		for (unsigned int i = 0 ; i < num_neg ; i ++)
			neg_vals(i) = vals(order.index[i]);

		return(neg_vals);
	}

	void Base_MatrixOpt::print_eigen_mat(Full_matrix &full_matrix,
											std::string matrix_name)
	{
		const unsigned int n_rows = full_matrix.rows();
		const unsigned int n_cols = full_matrix.cols();

		std::cout << matrix_name << std::endl;
		fflush(stdout);

		for (unsigned int i = 0 ; i < n_rows; i ++)
		{
			for (unsigned int j = 0 ; j < n_cols ; j++)
				std::cout << " " << full_matrix.coeffRef(i,j) ;

			std::cout << std::endl;

		}

	}


	void Base_MatrixOpt::print_eigen_mat(Sparse_matrix &full_matrix,
											std::string matrix_name)
	{
		const unsigned int n_rows = full_matrix.rows();
		const unsigned int n_cols = full_matrix.cols();

		std::cout << matrix_name << std::endl;
		fflush(stdout);
		
		for (unsigned int i = 0 ; i < n_rows; i ++)
		{
			for (unsigned int j = 0 ; j < n_cols ; j++)
				std::cout << " " << full_matrix.coeffRef(i,j) ;

			std::cout << std::endl;

		}

	}


	Full_matrix Base_MatrixOpt::compute_X(Sparse_matrix &A)
	{
		EigenSolver<MatrixXd> ES(A);

		// compute the eigenvectors of the matrix A
		Full_matrix vecs = ES.pseudoEigenvectors();
		return(vecs);
	}

	Full_matrix Base_MatrixOpt::compute_X(Full_matrix &A)
	{
		EigenSolver<MatrixXd> ES(A);

		// compute the eigenvectors of the matrix A
		Full_matrix vecs = ES.pseudoEigenvectors();
		return(vecs);
	}


	VectorXd Base_MatrixOpt::compute_Lambda(Sparse_matrix &A)
	{
		EigenSolver<MatrixXd> ES(A);

		// compute the eigenvectors of the matrix A
		VectorXd vals = ES.pseudoEigenvalueMatrix().diagonal();

		return(vals);
	}

	VectorXd Base_MatrixOpt::compute_Lambda(Full_matrix &A)
	{
		EigenSolver<MatrixXd> ES(A);

		// compute the eigenvectors of the matrix A
		VectorXd vals = ES.pseudoEigenvalueMatrix().diagonal();

		return(vals);
	}

	Full_matrix Base_MatrixOpt::compute_Amod(Sparse_matrix &A)
	{
		Full_matrix vecs = compute_X(A);
		VectorXd vals = compute_Lambda(A);
		Full_matrix Amod;
		Amod = vecs*vals.cwiseAbs().asDiagonal()*vecs.inverse();

		return(Amod);
	}

	// same as above but for full matrix
	Full_matrix Base_MatrixOpt::compute_Amod(Full_matrix &A)
	{
		Full_matrix vecs = compute_X(A);
		VectorXd vals = compute_Lambda(A);
		Full_matrix Amod;
		Amod = vecs*vals.cwiseAbs().asDiagonal()*vecs.inverse();

		return(Amod);
	}

	// compute the Upwind Matrix
	Full_matrix Base_MatrixOpt::compute_Aminus(Sparse_matrix &A)
	{
		Full_matrix X = compute_X(A);
		VectorXd Lambda = compute_Lambda(A);


		return(X * (Lambda.cwiseAbs() - Lambda).asDiagonal() * X.inverse());
	}
}