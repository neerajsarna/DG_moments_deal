// exact solution for poisson heat conduction
namespace ExactSolution
{
	using namespace dealii;

	template<int dim>
	class
	G20_PoissonHeat:public Base_ExactSolution<dim>
	{
		public:
			G20_PoissonHeat(const constant_data &constants,const Sparse_matrix &S_half);

			virtual void vector_value(const Point<dim> &p,Vector<double> &value) const ;

			virtual double s_r(const double r,const double phi) const;
			virtual double s_phi(const double r,const double phi) const;
			virtual double thetaP(const double r,const double phi) const;	
			virtual double R_rr(const double r,const double phi) const;
			virtual double R_thetatheta(const double r,const double phi) const;
			virtual double R_rtheta(const double r,const double phi) const;
			virtual double R_zz(const double r,const double phi) const;
	};

	template<>
	G20_PoissonHeat<2>::G20_PoissonHeat(const constant_data &constants,const Sparse_matrix &S_half)
	:
	Base_ExactSolution<2>(constants,S_half)
	{

		Assert(constants.nEqn == 13,ExcMessage("Only a solution for G20"));
	}


	template<>
	G20_PoissonHeat<1>::G20_PoissonHeat(const constant_data &constants,const Sparse_matrix &S_half)
	:
	Base_ExactSolution<1>(constants,S_half)
	{

		Assert(constants.nEqn == 6,ExcMessage("Only a solution for G20"));
	}

	template<>
	void
	G20_PoissonHeat<2>::vector_value(const Point<2> &p,Vector<double> &value) const
	{
		// first we check the size of the value vector
		Assert((int)value.size() == this->constants.nEqn,ExcNotInitialized());
		Assert(fabs(this->constants.alpha - 0.816496580927726) < 1e-5,ExcMessage("Exact Solution does not correspond to the given value of alpha"));

		// variables for which we need the 
		const unsigned int ID_theta = this->constants.variable_map.find("theta")->second;
		const unsigned int ID_heat = this->constants.variable_map.find("qy")->second;
		const unsigned int ID_stress = this->constants.variable_map.find("sigmayy")->second;
		const double y = p(1);

		value = 0;

		Assert(ID_theta == 3,ExcMessage("Wrong ID for theta"));

		if(fabs(this->constants.alpha-0.816496580927726) < 1e-5)
		{
			if (fabs(this->constants.tau - 0.1) < 1e-5 )
			{

				value[ID_theta] = -1.2829983743239428 + 0.0003362282005094463*cosh(7.453559924999298*y) - 
   								0.026127890589687234*pow(y,2) + 0.408248290463863*pow(y,4);

				value[ID_stress] = -(pow(5,0.5)*(-0.0013576450198781716 + 
        							0.0004853036051831808*cosh(7.453559924999298*y) - 0.03771236166328255*pow(y,2)
        							))/2.;

			
			}

			if (fabs(this->constants.tau - 0.3) < 1e-5)
			{
				value[ID_theta] = -1.2898753498434943 + 0.02291602260455048*cosh(2.484519974999766*y) - 
   									0.0783836717690617*pow(y,2) + 0.13608276348795437*pow(y,4);

				value[ID_stress] = -(pow(5,0.5)*(-0.03665641553671062 + 0.03307642954873192*cosh(2.484519974999766*y) - 
        								0.1131370849898476*pow(y,2)))/2.;

			}			
		}

		// The above values correspond to a unsymmetric system, therefore we now need to accomodate the symmetric system
		MatrixOpt::Base_MatrixOpt matrix_opt;
		value = matrix_opt.Sparse_matrix_dot_Vector(this->S_half,value);
	}

	template<>
	void
	G20_PoissonHeat<1>::vector_value(const Point<1> &p,Vector<double> &value) const
	{
		// first we check the size of the value vector
		Assert((int)value.size() == this->constants.nEqn,ExcNotInitialized());
		Assert(fabs(this->constants.alpha) < 1e-5,ExcMessage("Exact Solution does not correspond to the given value of alpha"));
		Assert(fabs(this->constants.theta0+1) < 1e-5,ExcMessage("Incorrect temperature value"));
		Assert(fabs(this->constants.theta1-1) < 1e-5,ExcMessage("Incorrect temperature value"));
		Assert(fabs(this->constants.tau-0.1) < 1e-5,ExcMessage("Incorrect tau value"));

		// variables for which we need the exact solution
		const unsigned int ID_theta = this->constants.variable_map_1D.find("theta")->second;
	
		AssertDimension(ID_theta,2);
		const double x = p(0);

		value = 0;

		if(fabs(this->constants.alpha) < 1e-5)
			if (fabs(this->constants.tau - 0.1) < 1e-5 )
				value[ID_theta] = 0.9950211932019228*x + 0.00974532392134874*sinh(4.47213595499958*x);
			
		

		// The above values correspond to a unsymmetric system, therefore we now need to accomodate the symmetric system
		MatrixOpt::Base_MatrixOpt matrix_opt;
		value = matrix_opt.Sparse_matrix_dot_Vector(this->S_half,value);
	}
	// need to fix the following routines
	template<int dim>
	double
	G20_PoissonHeat<dim>::s_r(const double r,const double phi)const
	{
		return 0;

	}

	template<int dim>
	double
	G20_PoissonHeat<dim>::s_phi(const double r,const double phi)const
	{

		return 0;

	}

	template<int dim>
	double
	G20_PoissonHeat<dim>::thetaP(const double r,const double phi)const
	{
		return 0;

	}

	template<int dim>
	double
	G20_PoissonHeat<dim>::R_rr(const double r,const double phi)const
	{
		return 0;

	}
	
	template<int dim>
	double
	G20_PoissonHeat<dim>::R_thetatheta(const double r,const double phi)const
	{
		return 0;

	}

	template<int dim>
	double
	G20_PoissonHeat<dim>::R_rtheta(const double r,const double phi)const
	{
		return 0;

	}

	template<int dim>
	double
	G20_PoissonHeat<dim>::R_zz(const double r,const double phi)const
	{
		return 0;

	}

}