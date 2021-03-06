namespace ExactSolution
{
	using namespace dealii;

	template<int dim>
	class
	Base_ExactSolution:public Function<dim>
	{
		public:
			Base_ExactSolution(const constant_numerics &constants,const Sparse_matrix &S_half,const int nEqn,
							   const int Ntensors);

			const Sparse_matrix S_half;
			virtual void vector_value(const Point<dim> &p,Vector<double> &value) const = 0;

			const constant_numerics constants;
			const int nEqn;
			const int Ntensors;

			double BI(const int n,const double x) const; 
			double BK(const int n,const double x) const;

		// the following routines are system dependent
			virtual double s_r(const double r,const double phi) const = 0;
			virtual double s_phi(const double r,const double phi) const = 0;
			virtual double thetaP(const double r,const double phi) const  = 0;	

			virtual double R_rr(const double r,const double phi) const = 0;
			virtual double R_thetatheta(const double r,const double phi) const = 0;
			virtual double R_rtheta(const double r,const double phi) const = 0;
			virtual double R_zz(const double r,const double phi) const = 0;

			double C0, C1, C2, C3, C4, K1, K2, K3, K4, K7, K8, 
	     			lambda1, lambda2, lambda3;
	};

	template<int dim>
	Base_ExactSolution<dim>::Base_ExactSolution(const constant_numerics &constants,
											    const Sparse_matrix &S_half,
											    const int nEqn,
											    const int Ntensors)
	:
	Function<dim>(nEqn),
	S_half(S_half),
	constants(constants),
	nEqn(nEqn),
	Ntensors(Ntensors)
	{;}

	template<int dim>
	double
	Base_ExactSolution<dim>::BI(const int n,const double x) const
	{
		return boost::math::cyl_bessel_i(n,x);
	}

	template<int dim>
	double
	Base_ExactSolution<dim>::BK(const int n,const double x) const
	{
		return boost::math::cyl_bessel_k(n,x);
	}

}