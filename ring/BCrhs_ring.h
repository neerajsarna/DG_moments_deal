// a general setting for BCrhs
namespace BCrhs
{
	using namespace dealii;

	// BCrhs for a ring with characteristic variables
	template<int dim>
	class
	BCrhs_ring_char:public Base_BCrhs
	{
		public:
			BCrhs_ring_char(const constant_data &constants,
					        const Sparse_matrix &B);

			virtual void BCrhs(const Tensor<1,dim,double> p,
							  const Tensor<1,dim,double> normal_vector,
							   Vector<double> &bc_rhs);
	}

	template<int dim>
	BCrhs_ring_char<dim>::BCrhs_ring_char(const constant_data &constants,
										  const Sparse_matrix &B)
	:
	Base_BCrhs(constants,B),
	{;};

	template<int dim>
	void
	BCrhs_ring_char<dim>::BCrhs(const Tensor<1,dim,double> p,
							  const Tensor<1,dim,double> normal_vector,
							   Vector<double> &bc_rhs)
	{
			AssertDimension(bc_rhs.size(),this->constants.nBC);
			Assert(dim > 1,ExcNotImplemented());

			double x_cord = p[0];
			double y_cord = p[1];
			double thetaW;
			const unsigned int ID_theta = dim + 1;

			const double norm = p.norm();

			// temperature of the outer wall
			if (norm > 0.7)
				thetaW = this->constants.theta1;

			// temperature of the inner wall
			else
				thetaW = this->constants.theta0;

			bc_rhs = 0;

			for (unsigned int m = 0 ; m < B.outerSize() ; m++)
				for (Sparse_matrix::InnerIterator n(B,m); n ; ++n)
				{
 			   		// only provide a boundary value for the temperature and no condition for velocity equation
					if (n.col() == ID_theta && n.row() > 0)
						bc_rhs(n.row()) = -sqrt(3.0/2.0) * thetaW * n.value();
				}
	}

	template<int dim>
	class
	BCrhs_ring_odd:public Base_BCrhs<dim>
	{
		public:
			BCrhs_ring_odd(const constant_data &constants,
					       Sparse_matrix &B,
					       const MatrixUI &odd_ID);

			
			const MatrixUI odd_ID;
			virtual void BCrhs(const Tensor<1,dim,double> p,
							  const Tensor<1,dim,double> normal_vector,
							   Vector<double> &bc_rhs);
	}

	template<int dim>
	BCrhs_ring_odd<dim>::BCrhs_ring_odd(const constant_data &constants,
					       				Sparse_matrix &B,
					       				const MatrixUI &odd_ID)
	:
	Base_BCrhs(constants,B),
	odd_ID(odd_ID)
	{;}

	template<int dim>
	void
	BCrhs_ring_odd<dim>>:BCrhs(const Tensor<1,dim,double> p,
							  const Tensor<1,dim,double> normal_vector,
							   Vector<double> &bc_rhs)
	{
		// for this boundary implementation, we need to have a bc_rhs the size 
		// of which equals the number of equations
		AssertDimension(bc_rhs.size(),this->constants.nEqn);
		double x_cord = p[0];
		double y_cord = p[1];
		double thetaW;
		const unsigned int ID_theta = dim + 1;

		const double norm = p.norm();
			
			// temperature of the outer wall
		if (norm > 0.7)
			thetaW = this->constants.theta1;

			// temperature of the inner wall
		else
			thetaW = this->constants.theta0;

			for (unsigned int m = 0 ; m < B.outerSize() ; m++)
				for (Sparse_matrix::InnerIterator n(B,m); n ; ++n)
				{
					const double odd_coeff = B.coeffRef(n.row(),odd_ID(n.row,0));

					// the assumption in this implementation is the order in which the variables
					// appear in the odd_ID is the same order in which the boundary conditions have been specified.

					Assert(fabs(odd_coeff) > 1e-5,ExcMessage("Assumption broken down "));
 			   		// only provide a boundary value for the temperature and no condition for velocity equation
					if (n.col() == ID_theta && n.row() > 0)
						bc_rhs(odd_ID(n.row(),0)) = -sqrt(3.0/2.0) * thetaW * n.value()/odd_coeff;
				}
	}

}