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

	// we would now like to look into the development of odd Boundary Conditions
	template<int dim>
	class
	Base_BoundaryHandler_Odd
	{
		public:
			Base_BoundaryHandler_Odd(Sparse_matrix &B,const MatrixUI &odd_ID);

			Sparse_matrix B;
			const MatrixUI odd_ID;

			// the boundary conditions will be implemented in the form U = BCU+g
			Sparse_matrix develop_BC();
	};

	template<int dim>
	Base_BoundaryHandler_Odd<dim>::Base_BoundaryHandler_Odd(Sparse_matrix &B,const MatrixUI &odd_ID)
	:
	B(B),
	odd_ID(odd_ID)
	{
		Assert(B.rows() != 0 || B.cols() != 0,ExcNotInitialized());
		Assert(odd_ID.rows() != 0 || odd_ID.cols() != 0,ExcNotInitialized());
	}

	// we now try to develop the BC matrix for the implementation of the odd boundary conditions
	template<int dim>
	Sparse_matrix
	Base_BoundaryHandler_Odd<dim>::develop_BC()
	{
		const unsigned int nEqn = B.cols();
		unsigned int num_odd = 1;
		Sparse_matrix BC;
		BC.resize(nEqn,nEqn);

		for (unsigned int i = 0 ; i < nEqn ; i++)
		{
			//if we encounter an odd variable
			if (i == odd_ID(num_odd - 1,0))
			{
																									
				const double odd_coeff = B.coeffRef(num_odd-1,odd_ID(num_odd - 1,0));

				// We now need to check whether the order in which the odd variables are listed is the
				// same order in which the boundary conditions for these odd variables have been loaded.
				Assert(fabs(odd_coeff) > 1e-5,ExcMessage("Assumption has broken down"));
				
				//If we write the boundary conditions as A.u_odd + C.u_even = g then u_odd = -Inverse(A).C.u_even
				//A = B.coeffRef(num_odd,odd_ID(num_odd)) and C = B.coeffRef(num_odd,j). 
				for (unsigned int j = 0 ; j < B.cols() ; j++)
					if (j != odd_ID(num_odd - 1,0))	
						BC.coeffRef(odd_ID(num_odd - 1,0),j) = -B.coeffRef(num_odd - 1,j)/odd_coeff;
				
				// // // update only if it doest now exceed the size limit
				
				if (num_odd < B.rows())
					num_odd ++;
			}
			else
				BC.coeffRef(i,i) = 1;
		}

		BC.makeCompressed();
		B.makeCompressed();
		AssertDimension(B.rows(),num_odd);

		return(BC);
	}

}