namespace BCrhs
{
	using namespace dealii;

	// an abstract class for BCrhs
	template<int dim>
	class
	Base_BCrhs
	{
		public:
			Base_BCrhs(const constant_data &constants);
			const constant_data constants;
			virtual void BCrhs(const Tensor<1,dim,double> p,
							  const Tensor<1,dim,double> normal_vector,
							   Vector<double> &bc_rhs) = 0;
	};

	template<int dim>
	Base_BCrhs<dim>
	::Base_BCrhs(const constant_data &constants)
	:
	constants(constants)
	{;}
}