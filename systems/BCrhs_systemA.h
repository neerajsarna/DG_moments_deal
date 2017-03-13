namespace BCrhs_systemA
{
	using namespace dealii;

	// BCrhs for a ring with characteristic variables
	template<int dim>
	class
	BCrhs_ring_char_systemA:public BCrhs::Base_BCrhs<dim>
	{
		public:
			BCrhs_ring_char_systemA(const constant_numerics &constants,const int nBC);

			virtual void BCrhs(const Tensor<1,dim,double> p,
							  const Tensor<1,dim,double> normal_vector,
							   Vector<double> &bc_rhs,
							   const unsigned int b_id);

			const int nBC;
	};

	template<int dim>
	BCrhs_ring_char_systemA<dim>::BCrhs_ring_char_systemA(const constant_numerics &constants,const int nBC)
	:
	BCrhs::Base_BCrhs<dim>(constants),
	nBC(nBC)
	{;}

	template<int dim>
	void
	BCrhs_ring_char_systemA<dim>::BCrhs(const Tensor<1,dim,double> p,
								  		const Tensor<1,dim,double> normal_vector,
							   			Vector<double> &bc_rhs,
							   			const unsigned int b_id)
	{
			AssertDimension((int)bc_rhs.size(),this->constants.nBC);

			// the implementation is only for a 2D setting
			AssertDimension(dim,2);
			Assert(dim > 1,ExcNotImplemented());

			double x_cord = p[0];
			double y_cord = p[1];
			
			const double norm = p.norm();

			// implementation for the outer boundary
			if( norm > 0.7 ) 
			{
		 		bc_rhs(0) = -this->constants.chi * this->constants.theta1;  // is chi \alpha and same with zeta
			}
			

		 	// implementation for the inner boundary
		 	else 
		 	{
		 		bc_rhs(0) = -this->constants.chi * this->constants.theta0;
		 		bc_rhs(1) = -this->constants.uW * normal_vector[1];
		 	}; 
		 		

	}

	// BCrhs with odd variables
	template<int dim>
	class
	BCrhs_ring_odd_systemA:public BCrhs::Base_BCrhs<dim>
	{
		public:
			BCrhs_ring_odd_systemA(const constant_numerics &constants,const int nBC);

			virtual void BCrhs(const Tensor<1,dim,double> p,
							  const Tensor<1,dim,double> normal_vector,
							   Vector<double> &bc_rhs,
							   const unsigned int b_id);

			const int nBC;
	};

	template<int dim>
	BCrhs_ring_odd_systemA<dim>::BCrhs_ring_odd_systemA(const constant_numerics &constants,const int nBC)
	:
	BCrhs::Base_BCrhs<dim>(constants),
	nBC(nBC)
	{;}

	template<int dim>
	void
	BCrhs_ring_odd_systemA<dim>::BCrhs(const Tensor<1,dim,double> p,
								  		const Tensor<1,dim,double> normal_vector,
							   			Vector<double> &bc_rhs,
							   			const unsigned int b_id)
	{
			AssertDimension((int)bc_rhs.size(),this->constants.nEqn);
			Assert(dim > 1,ExcNotImplemented());

			double x_cord = p[0];
			double y_cord = p[1];
			
			const double norm = p.norm();

			if( norm > 0.7 ) 
 		   			bc_rhs(1) = -this->constants.chi * this->constants.theta1;  // is chi \alpha and same with zeta

 		   		else 
 		   		{
 		   			bc_rhs(1) = -this->constants.chi * this->constants.theta0;
 		   			bc_rhs(4) = -this->constants.uW*normal_vector[1];
 		   		}; 


	}

}