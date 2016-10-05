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