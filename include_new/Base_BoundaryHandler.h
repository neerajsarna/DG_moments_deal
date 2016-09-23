// The following routine handles the boundary conditions

namespace BoundaryHandler
{
	using namespace dealii;

	// develops the boundary matrices needed for characteristic boundary conditions
	template<int dim>
	class
	Base_BoundaryHandler_Char:public MatrixOpt::Base_MatrixOpt
	{
		public:
			Base_BoundaryHandler_Char(const Sparse_matrix Ax,Sparse_matrix &B,const unsigned int nBC);

			const unsigned int nBC;
			// Jacobian in the x-direction
			Sparse_matrix Ax;

			// The boundary matrix from BU = g
			Sparse_matrix *B;

			// The matrix containing the eigenvectors of negative eigenvalues
			Full_matrix X_minus;

			// fix the boundary conditions for vx by introducing the epsilon .
			// If action == true only then B is fixed, else it's not fixed
			void fix_B_vx(	const double epsilon,
					 const bool action);

			// build Inv(B.X_minus)
			Full_matrix build_B_tilde_inv() const;

			// build X_minus.inverse(B_tilde).B
			Full_matrix build_B_hat(const Full_matrix &B_tilde_inv) const;

									
	};

	// the constructor simply stores the matrices which will be needed in the future
	template<int dim>
	Base_BoundaryHandler_Char<dim>
	::
	Base_BoundaryHandler_Char(const Sparse_matrix Ax,Sparse_matrix &B,const unsigned int nBC)
	:
	nBC(nBC),
	Ax(Ax),
	B(&B)
	{
		X_minus = compute_Xminus(Ax,nBC);
	;}

	template<int dim>
	void 
	Base_BoundaryHandler_Char<dim>
	::fix_B_vx(const double epsilon,const bool action)
	{
		Assert(B->cols()!= 0 && B->rows()!= 0,ExcMessage("Matrix not initialized") );
		Assert(epsilon >= 0,ExcMessage("epsilon cant be negative"));

		// we first compute the number of equations in the system
		const unsigned int neqn_local = B->cols();

		// every coefficient apart from the normal velocity itself should be changed 
		if (action)
			for (unsigned int i = 0 ; i < neqn_local ; i++)
				if (i != 1)				// do nothing if the coefficient corresponds to the normal velocity						
					B->coeffRef(0,i) *= epsilon;
	}

	template<int dim>
	Full_matrix
	Base_BoundaryHandler_Char<dim>
	::build_B_tilde_inv() const
	{
		

		Assert(B->rows() !=0 && B->cols() !=0,ExcNotInitialized());
		Assert(X_minus.rows() !=0 && X_minus.cols() !=0,ExcNotInitialized());

		// B.rows = number of boundary conditions = number of negative eigenvalues
		// B.cols = number of eqatuions in the system
		// X.rows = number of equations in the system
		// X.cols = number of negative eigen values in the system
		AssertDimension(B->rows(),X_minus.cols());
		AssertDimension(B->cols(),X_minus.rows());

		return(((*B) * X_minus).inverse());
	}

	template<int dim>
	Full_matrix
	Base_BoundaryHandler_Char<dim>
	::build_B_hat(const Full_matrix &B_tilde_inv) const
	{
		Assert(B->rows() !=0 && B->cols() !=0,ExcNotInitialized());
		Assert(X_minus.rows() !=0 && X_minus.cols() !=0,ExcNotInitialized());

		// B.rows = number of boundary conditions = number of negative eigenvalues
		// B.cols = number of eqatuions in the system
		// X.rows = number of equations in the system
		// X.cols = number of negative eigen values in the system
		AssertDimension(B->rows(),X_minus.cols());
		AssertDimension(B->cols(),X_minus.rows());

		Assert(B_tilde_inv.rows()!=0 && B_tilde_inv.cols() != 0, ExcMessage("Matrix not initialized"));

		return(X_minus * B_tilde_inv * (*B));
	}

}