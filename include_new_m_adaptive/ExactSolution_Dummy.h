//Exact solution which returns zero value as the output
namespace ExactSolution
{
	using namespace dealii;

	template<int dim>
	class
	ExactSolution_Dummy:public Base_ExactSolution<dim>
	{
		public:
			ExactSolution_Dummy(const constant_numerics &constants, 
								const Sparse_matrix &S_half,const int nEqn,const int Ntensors);

			virtual void vector_value(const Point<dim> &p,Vector<double> &value) const ;

			virtual double s_r(const double r,const double phi) const;
			virtual double s_phi(const double r,const double phi) const;
			virtual double thetaP(const double r,const double phi) const;	
			virtual double R_rr(const double r,const double phi) const;
			virtual double R_thetatheta(const double r,const double phi) const;
			virtual double R_rtheta(const double r,const double phi) const;
			virtual double R_zz(const double r,const double phi) const;
	};

	template<int dim>
	ExactSolution_Dummy<dim>::ExactSolution_Dummy(const constant_numerics &constants,const Sparse_matrix &S_half,
												const int nEqn,const int Ntensors)
	:
	Base_ExactSolution<dim>(constants,S_half,nEqn,Ntensors)
	{
	}

	template<int dim>
	void
	ExactSolution_Dummy<dim>::vector_value(const Point<dim> &p,Vector<double> &value) const
	{
		Assert(value.size() == (unsigned int)this->nEqn,ExcMessage("Incorrect dimension of the input array "));
		Assert(p.norm() >=0 ,ExcMessage("Incorrect point"));

		// we just return zero value 
		value = 0;
	}

	// need to fix the following routines
	template<int dim>
	double
	ExactSolution_Dummy<dim>::s_r(const double r,const double phi)const
	{
		Assert(r >= 0,ExcMessage("Incorrect radial distance"));
		Assert(fabs(phi) >= 0,ExcMessage("Incorrect angle"));

		return 0;

	}

	template<int dim>
	double
	ExactSolution_Dummy<dim>::s_phi(const double r,const double phi)const
	{
		Assert(r >= 0,ExcMessage("Incorrect radial distance"));
		Assert(fabs(phi) >= 0,ExcMessage("Incorrect angle"));

		return 0;

	}

	template<int dim>
	double
	ExactSolution_Dummy<dim>::thetaP(const double r,const double phi)const
	{
		Assert(r >= 0,ExcMessage("Incorrect radial distance"));
		Assert(fabs(phi) >= 0,ExcMessage("Incorrect angle"));

		return 0;

	}

	template<int dim>
	double
	ExactSolution_Dummy<dim>::R_rr(const double r,const double phi)const
	{
		Assert(r >= 0,ExcMessage("Incorrect radial distance"));
		Assert(fabs(phi) >= 0,ExcMessage("Incorrect angle"));

		return 0;

	}
	
	template<int dim>
	double
	ExactSolution_Dummy<dim>::R_thetatheta(const double r,const double phi)const
	{
		Assert(r >= 0,ExcMessage("Incorrect radial distance"));
		Assert(fabs(phi) >= 0,ExcMessage("Incorrect angle"));

		return 0;

	}

	template<int dim>
	double
	ExactSolution_Dummy<dim>::R_rtheta(const double r,const double phi)const
	{
		Assert(r >= 0,ExcMessage("Incorrect radial distance"));
		Assert(fabs(phi) >= 0,ExcMessage("Incorrect angle"));

		return 0;

	}

	template<int dim>
	double
	ExactSolution_Dummy<dim>::R_zz(const double r,const double phi)const
	{
		Assert(r >= 0,ExcMessage("Incorrect radial distance"));
		Assert(fabs(phi) >= 0,ExcMessage("Incorrect angle"));
		
		return 0;

	}

}