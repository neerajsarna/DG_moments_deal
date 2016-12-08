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

			//compute the dydadic product of A and B i.e 
			Full_matrix compute_A_outer_B(Full_matrix &A,Full_matrix &B);

			Sparse_matrix compute_A_outer_B(Full_matrix &A,Sparse_matrix &B);

			Sparse_matrix compute_A_outer_B(Sparse_matrix &A,Full_matrix &B);

			Sparse_matrix compute_A_outer_B(Sparse_matrix &A,Sparse_matrix &B);


			// prints a full matrix
			void print_eigen_mat(Full_matrix &full_matrix,
								std::string matrix_name);

			// same as above but for sparse matrix
			void print_eigen_mat(Sparse_matrix &full_matrix,
								std::string matrix_name);

			// print a matrix to a file
			void print_eigen_mat_to_file(Sparse_matrix &full_matrix,
										 std::string matrix_name);

			void print_eigen_mat_to_file(Full_matrix &full_matrix,
										 std::string matrix_name);

			// print dealii sparse matrix to a file
			void print_dealii_sparse(const TrilinosWrappers::SparseMatrix &matrix,
									std::string matrix_name);

			void print_dealii_vector(const Vector<double> &vec, 
									std::string vector_name);

			// a dot product between sparse matrix and a vector
			Vector<double> Sparse_matrix_dot_Vector(const Sparse_matrix &matrix,const Vector<double> &vec);
			void COO_to_CSR(const TrilinosWrappers::SparseMatrix &matrix, MKL_INT *IA,MKL_INT *JA,double *V);

			Vector<double> Sparse_matrix_dot_Vector(const MKL_INT *ia,const MKL_INT *ja,
                                                                 const double *values,const Vector<double> &vec,
                                                                const double size);
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

	// print the full matrix to a file
	void Base_MatrixOpt::print_eigen_mat_to_file(Full_matrix &full_matrix,
											std::string matrix_name)
	{
		const unsigned int n_rows = full_matrix.rows();
		const unsigned int n_cols = full_matrix.cols();

		FILE *fp;
		std::string filename = "printed_matrices/" + matrix_name;


		fp = fopen(filename.c_str(),"w+");

		AssertThrow(fp != NULL , ExcMessage("Could not open file for writting matrix "));

		for (unsigned int i = 0 ; i < n_rows; i ++)
		{
			for (unsigned int j = 0 ; j < n_cols ; j++)
				fprintf(fp, "%f ",full_matrix.coeffRef(i,j));

			fprintf(fp, "\n");

		}

		fclose(fp);

	}

	void Base_MatrixOpt::print_eigen_mat_to_file(Sparse_matrix &full_matrix,
											std::string matrix_name)
	{
		const unsigned int n_rows = full_matrix.rows();
		const unsigned int n_cols = full_matrix.cols();

		FILE *fp;
		std::string filename = "printed_matrices/" + matrix_name;


		fp = fopen(filename.c_str(),"w+");

		AssertThrow(fp != NULL , ExcMessage("Could not open file for writting matrix "));

		for (unsigned int i = 0 ; i < n_rows; i ++)
		{
			for (unsigned int j = 0 ; j < n_cols ; j++)
				fprintf(fp, "%f ",full_matrix.coeffRef(i,j));

			fprintf(fp, "\n");

		}

		fclose(fp);

	}

	void Base_MatrixOpt::print_dealii_sparse(const TrilinosWrappers::SparseMatrix &matrix,
										  std::string matrix_name)
	  {
	  	FILE *fp;
	  	std::string filename;
	  	filename = "printed_matrices/" + matrix_name;
	  	fp = fopen(filename.c_str(),"w+");

	  	AssertThrow(fp != NULL, ExcMessage("could not open file for writting global matrix"));

	  	typename TrilinosWrappers::SparseMatrix::const_iterator it = matrix.begin();
	  	const typename TrilinosWrappers::SparseMatrix::const_iterator it_end = matrix.end();

	  	for (; it != it_end; it++)
	  		fprintf(fp,"%llu %llu %lf \n",it->row(),it->column(),it->value());
	  	

	  	fclose(fp);


	  }

	  void Base_MatrixOpt::print_dealii_vector(const Vector<double> &vec,
	  										std::string vector_name)
	  {
	  	std::string filename = "printed_matrices/" + vector_name;

	  	FILE *fp;
	  	fp = fopen(filename.c_str(),"w+");

	  	for (unsigned int i = 0 ; i < vec.size() ; i ++)
	  		fprintf(fp, "%f\n",vec(i));

	  	fclose(fp);
	  }

	  // The following function performs the dot product between a dealii sparse matrix and a vector
	  Vector<double> Base_MatrixOpt::Sparse_matrix_dot_Vector(const Sparse_matrix &matrix,
	  														  const Vector<double> &vec)
	  {
	  	const unsigned int num_entries = vec.size();
	  	Assert(matrix.rows() != 0 || matrix.cols() != 0,ExcNotInitialized());
	  	Assert(vec.size() != 0,ExcNotInitialized());
	  	Assert(matrix.IsRowMajor,ExcMessage("Basic assumption failed"));

		Vector<double> result(num_entries);

		for (unsigned int m = 0 ; m < matrix.outerSize(); m++)
		{
			result(m) = 0;
			for (Sparse_matrix::InnerIterator n(matrix,m); n ; ++n)
				result(m) += n.value() * vec(n.col());
		}

		return(result);
	  }

	  // the following function expects result to be initialized before hand
	  Full_matrix Base_MatrixOpt::compute_A_outer_B(Full_matrix &A,Full_matrix &B)
	  {
	  	


	  	const unsigned int rows_A = A.rows();
	  	const unsigned int cols_A = A.cols();
	  	const unsigned int rows_B = B.rows();
	  	const unsigned int cols_B = B.cols();
	  	Full_matrix result(rows_A * rows_B,cols_A * cols_B);


	  	AssertDimension(result.rows(),rows_A * rows_B);
	  	AssertDimension(result.cols(),cols_A * cols_B);

	  	for (unsigned int i = 0 ; i < rows_A ; i ++)
	  		for (unsigned int j = 0 ; j < cols_A ; j++)
	  			result.block(i * rows_B, j * cols_B, rows_B,cols_B) = A(i,j) * B;

	  	return(result);
	  }

	  	  // the following function expects result to be initialized before hand
	  Sparse_matrix Base_MatrixOpt::compute_A_outer_B(Full_matrix &A,Sparse_matrix &B)
	  {
	  	
	  	const unsigned int rows_A = A.rows();
	  	const unsigned int cols_A = A.cols();
	  	const unsigned int rows_B = B.rows();
	  	const unsigned int cols_B = B.cols();
	  	Sparse_matrix result(rows_A * rows_B, cols_A * cols_B);

	  	AssertDimension(result.rows(),rows_A * rows_B);
	  	AssertDimension(result.cols(),cols_A * cols_B);

	  	for (unsigned int i = 0 ; i < rows_A ; i ++)
	  		for (unsigned int j = 0 ; j < cols_A ; j++)
	  			for (unsigned int m = 0 ; m < B.outerSize(); m++)
          			for (Sparse_matrix::InnerIterator n(B,m); n ; ++n)
	  					result.coeffRef(i * rows_B + n.row(),j * cols_B + n.col()) = A(i,j) * n.value();


	  	result.makeCompressed();
	  	return(result);
	  }

	  Sparse_matrix Base_MatrixOpt::compute_A_outer_B(Sparse_matrix &A,Full_matrix &B)
	  {

	  	const unsigned int rows_A = A.rows();
	  	const unsigned int cols_A = A.cols();
	  	const unsigned int rows_B = B.rows();
	  	const unsigned int cols_B = B.cols();
	  	Sparse_matrix result(rows_A * rows_B,cols_A * cols_B);


	  	AssertDimension(result.rows(),rows_A * rows_B);
	  	AssertDimension(result.cols(),cols_A * cols_B);

	  	for (unsigned int i = 0 ; i < A.outerSize(); i++)
          	for (Sparse_matrix::InnerIterator j(A,i); j ; ++j)
	  			for (unsigned int k = 0 ; k < rows_B ; k ++)
          			for (unsigned int l = 0 ; l < cols_B ; l ++)
	  					result.coeffRef(j.row() * rows_B + k,j.col() * cols_B + l) = j.value() * B(k,l);


	  	result.makeCompressed();
	  	return(result);
	  }

	  Sparse_matrix Base_MatrixOpt::compute_A_outer_B(Sparse_matrix &A,Sparse_matrix &B)
	  {
	 

	  	const unsigned int rows_A = A.rows();
	  	const unsigned int cols_A = A.cols();
	  	const unsigned int rows_B = B.rows();
	  	const unsigned int cols_B = B.cols();
	  	Sparse_matrix result(rows_A * rows_B, cols_A * cols_B);


	  	AssertDimension(result.rows(),rows_A * rows_B);
	  	AssertDimension(result.cols(),cols_A * cols_B);

	  	for (unsigned int i = 0 ; i < A.outerSize(); i++)
          	for (Sparse_matrix::InnerIterator j(A,i); j ; ++j)
	  			for (unsigned int k = 0 ; k < B.outerSize(); k++)
          			for (Sparse_matrix::InnerIterator l(B,k); l ; ++l)
	  					result.coeffRef(j.row() * rows_B + l.row(),j.col() * cols_B + l.col())
	  											 = j.value() * l.value();


	  	result.makeCompressed();
	  	return(result);
	  }

	  void Base_MatrixOpt::COO_to_CSR(const TrilinosWrappers::SparseMatrix &matrix, MKL_INT *IA,MKL_INT *JA,double *V)
	  {

	  	std::cout << "Developing CSR " << std::endl;
	  	fflush(stdout);

	  	const unsigned int n_rows = matrix.m();
	  	const unsigned int nnz = matrix.n_nonzero_elements();

	  	// we check the sizes of the various input pointers
	  	Assert(sizeof(IA)/sizeof(MKL_INT) != 0,ExcNotInitialized());
		Assert(sizeof(JA)/sizeof(MKL_INT) != 0,ExcNotInitialized());
		Assert(sizeof(V)/sizeof(double) != 0,ExcNotInitialized());

	  	typename TrilinosWrappers::SparseMatrix::const_iterator it = matrix.begin(),
	  														    it_end = matrix.end();

	  	typename TrilinosWrappers::SparseMatrix::const_iterator it_dummy = matrix.begin();


	  	 unsigned int counter = 1;
	  	 std::vector<unsigned int> row_index(nnz);
	  	 unsigned int row_counter = 1;
	  	 unsigned int rows_captured = 1;
	  	 IA[0] = 0;
	  	 JA[0] = it->column();
	  	 V[0] = it->value();
	  	 row_index[0] = it->row();
	  	 it++;

	  	for (; it != it_end ; it++)
	  	{
	  		Assert(counter <= nnz ,ExcMessage("Inappropriate counter"));
	  		
	  		V[counter] = it->value();
	  		JA[counter] = it->column();
	  		row_index[counter] = it->row();

	  		Assert(row_index[counter-1] <= row_index[counter],ExcMessage("row indices not sorted"));

	  			// if the row of the entries is the same then we increase the counter
	  		if(row_index[counter] == row_index[counter-1])
	  			row_counter++;

	  			// if the rows of the entries is not the same the new add this location to IA
	  		else
	  		{
	  			// set the counter back to default
	  			IA[rows_captured] = row_counter + IA[rows_captured-1];

	  			rows_captured ++;
	  			row_counter = 1;
	  		}

	  		Assert(rows_captured <= n_rows,ExcMessage("Too many rows captured"));
	  		

	  		counter ++;
	  		
	  	}

	  	Assert(rows_captured == n_rows,ExcMessage("Incorrect number of rows captured"));

	  	// the final entry
	  	IA[n_rows] = nnz;
	  	row_index.clear();

	  	std::cout << "Finished Developing CSR "<< std::endl;
	  	
	  }


	  // size the number of rows in the matrix
	  Vector<double> Base_MatrixOpt::Sparse_matrix_dot_Vector(const MKL_INT *ia,const MKL_INT *ja,
                                                                 const double *values,const Vector<double> &vec,
                                                                const double size)
        {
                Assert(sizeof(ia) != 0,ExcNotInitialized());
                Assert(sizeof(ja) != 0,ExcNotInitialized());
                Assert(sizeof(values) != 0,ExcNotInitialized());
                Assert(vec.size() !=0,ExcNotInitialized());

                Vector<double> result(size);
                int k = 0;


                for (unsigned int i = 0 ; i < size ; i++)
                        for (unsigned int j = ia[i] ; j < ia[i + 1] ; j++)
                        {
                        		// this gives us the column number
                                k = ja[j];

                                // values[j] gives us the entry of the matrix and 
                                // vec(k) gives us the entry of the input vecotr
                                result(i) += values[j] * vec(k);


                        }


                return(result);
        }
}
