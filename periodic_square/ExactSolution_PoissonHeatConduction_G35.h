// exact solution for poisson heat conduction
namespace ExactSolution
{
	using namespace dealii;

	template<int dim>
	class
	G35_PoissonHeat:public Base_ExactSolution<dim>
	{
		public:
			G35_PoissonHeat(const constant_data &constants,const Sparse_matrix &S_half);

			virtual void vector_value(const Point<dim> &p,Vector<double> &value) const ;

			virtual double s_r(const double r,const double phi) const;
			virtual double s_phi(const double r,const double phi) const;
			virtual double thetaP(const double r,const double phi) const;	

			virtual double R_rr(const double r,const double phi) {return 0;};
			virtual double R_thetatheta(const double r,const double phi) {return 0;};
			virtual double R_rtheta(const double r,const double phi) {return 0;};
			virtual double R_zz(const double r,const double phi) {return 0;};
	};

	template<int dim>
	G35_PoissonHeat<dim>::G35_PoissonHeat(const constant_data &constants,const Sparse_matrix &S_half)
	:
	Base_ExactSolution<dim>(constants,S_half)
	{
		Assert(constants.nEqn == 22,ExcMessage("Only a solution for G35"));
	}

	template<int dim>
	void
	G35_PoissonHeat<dim>::vector_value(const Point<dim> &p,Vector<double> &value) const
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

		Assert(ID_theta == dim+1,ExcMessage("Wrong ID for theta"));

		if(fabs(this->constants.alpha-0.816496580927726) < 1e-5)
		{
			if (fabs(this->constants.tau - 0.1) < 1e-5 )
			{

				value[ID_theta] = -1.2857735434320117 + 0.002601070885430934*cosh(4.662524041201569*y) - 
   									0.18289523412781067*pow(y,2) + 0.408248290463863*pow(y,4) + 
   									2.168404344971009e-19*sinh(4.662524041201569*y);

				value[ID_stress] = -(pow(5,0.5)*(-0.002715290039756343 + 0.00187716121985606*cosh(4.662524041201569*y) - 
        							0.03771236166328255*pow(y,2) - 
        							1.0842021724855044e-19*sinh(4.662524041201569*y)))/2.;

			}

			if (fabs(this->constants.tau - 0.3) < 1e-5)
			{
				value[ID_theta] = -1.3660400339849652 + 0.09916339949943631*cosh(1.554174680400523*y) - 
   									0.548685702383432*pow(y,2) + 0.13608276348795437*pow(y,4) + 
   									6.938893903907228e-18*sinh(1.554174680400523*y);

				value[ID_stress] = -(pow(5,0.5)*(-0.07331283107342124 + 0.07156501924344744*cosh(1.554174680400523*y) - 
        							0.1131370849898476*pow(y,2)))/2.;
				
			}		
		}

		// The above values correspond to a unsymmetric system, therefore we now need to accomodate the symmetric system
		MatrixOpt::Base_MatrixOpt matrix_opt;
		value = matrix_opt.Sparse_matrix_dot_Vector(this->S_half,value);
	}

	// need to fix the following routines
	template<int dim>
	double
	G35_PoissonHeat<dim>::s_r(const double r,const double phi)const
	{
		return 0;

	}

	template<int dim>
	double
	G35_PoissonHeat<dim>::s_phi(const double r,const double phi)const
	{

				return 0;
	}

	template<int dim>
	double
	G35_PoissonHeat<dim>::thetaP(const double r,const double phi)const
	{
		return 0;

	}
}