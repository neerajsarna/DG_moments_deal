using namespace dealii;

// we define certain basic routines which are needed everywhere
// compare tow flow matrices
void Compare_Float_Mat(Full_matrix &A,Full_matrix &B)
{
	// the dimensions should be the same for comparision
	ASSERT_EQ(A.rows(),B.rows());
	ASSERT_EQ(A.cols(),B.cols());

	const unsigned int rows = A.rows();
	const unsigned int cols = A.cols();

	for (unsigned int i = 0 ; i < rows ; i ++)
		for (unsigned int j = 0 ; j < cols ; j++)
			EXPECT_NEAR(A.coeffRef(i,j),B.coeffRef(i,j),1e-5);
}

// same as above but for sparse matrices
void Compare_Float_Mat(Sparse_matrix &A,Sparse_matrix &B)
{
	// the dimensions should be the same for comparision
	ASSERT_EQ(A.rows(),B.rows());
	ASSERT_EQ(A.cols(),B.cols());

	const unsigned int rows = A.rows();
	const unsigned int cols = A.cols();

	for (unsigned int i = 0 ; i < rows ; i ++)
		for (unsigned int j = 0 ; j < cols ; j++)
			EXPECT_NEAR(A.coeffRef(i,j),B.coeffRef(i,j),1e-5);
}

// same as above but for sparse matrices
void Compare_Float_Mat(Sparse_matrix &A,Full_matrix &B)
{
	// the dimensions should be the same for comparision
	ASSERT_EQ(A.rows(),B.rows());
	ASSERT_EQ(A.cols(),B.cols());

	const unsigned int rows = A.rows();
	const unsigned int cols = A.cols();

	for (unsigned int i = 0 ; i < rows ; i ++)
		for (unsigned int j = 0 ; j < cols ; j++)
			EXPECT_NEAR(A.coeffRef(i,j),B.coeffRef(i,j),1e-5);
}

// same as above but for sparse matrices
void Compare_Float_Mat(Full_matrix &A,Sparse_matrix &B)
{
	// the dimensions should be the same for comparision
	ASSERT_EQ(A.rows(),B.rows());
	ASSERT_EQ(A.cols(),B.cols());

	const unsigned int rows = A.rows();
	const unsigned int cols = A.cols();

	for (unsigned int i = 0 ; i < rows ; i ++)
		for (unsigned int j = 0 ; j < cols ; j++)
			EXPECT_NEAR(A.coeffRef(i,j),B.coeffRef(i,j),1e-5);
}

// same as above but for sparse matrices
void Compare_Float_Mat(TrilinosWrappers::SparseMatrix &A,Full_matrix &B)
{
	// the dimensions should be the same for comparision
	ASSERT_EQ(A.m(),B.rows());
	ASSERT_EQ(A.n(),B.cols());

	const unsigned int rows = A.m();
	const unsigned int cols = A.n();

	typename TrilinosWrappers::SparseMatrix::const_iterator it = A.begin();
	const typename TrilinosWrappers::SparseMatrix::const_iterator it_end = A.end();

	for (; it!=it_end;it++)
		EXPECT_NEAR(it->value(),B.coeffRef(it->row(),it->column()),1e-6)<< "value of row: " << it->row() << " value of col: "<< it->column();
}


// compare two float vectors
// Both vectors from dealii
void Compare_Float_Vec(dealii::Vector<double> &vec1,dealii::Vector<double> &vec2)
{
	AssertDimension(vec1.size(),vec2.size());

	for (unsigned int i = 0 ; i < vec1.size() ; i++)
		EXPECT_NEAR(vec1(i),vec2(i),1e-5);
}


