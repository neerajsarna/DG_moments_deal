// exact solution for poisson heat conduction
namespace ExactSolution
{
	using namespace dealii;

	template<int dim>
	class
	G45_PoissonHeat:public Base_ExactSolution<dim>
	{
		public:
			G45_PoissonHeat(const constant_data &constants,const Sparse_matrix &S_half);

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
	G45_PoissonHeat<dim>::G45_PoissonHeat(const constant_data &constants,const Sparse_matrix &S_half)
	:
	Base_ExactSolution<dim>(constants,S_half)
	{
		//Assert(constants.nEqn == 28,ExcMessage("Only a solution for G45"));
	}

	template<int dim>
	void
	G45_PoissonHeat<dim>::vector_value(const Point<dim> &p,Vector<double> &value) const
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

				value[ID_theta] = -1.2995436998701937 + 0.012092296117067552*cosh(3.849386397950811*y) + 
   0.0015393861867350283*cosh(5.674996469312898*y) + 
   0.000022293298256729303*cosh(10.161219234535377*y) - 
   0.18289523412781067*pow(y,2) + 0.408248290463863*pow(y,4) - 
   8.673617379884035e-19*sinh(3.849386397950811*y) - 
   5.421010862427522e-19*sinh(5.674996469312898*y) - 
   8.470329472543003e-21*sinh(10.161219234535377*y);

				value[ID_stress] = -(pow(5,0.5)*(-0.002715290039756343 + 
        0.0031225329963700333*cosh(3.849386397950811*y) - 
        0.0007338661922650514*cosh(5.674996469312898*y) + 
        0.000033602606211988805*cosh(10.161219234535377*y) - 
        0.03771236166328255*pow(y,2) - 
        1.5178830414797062e-18*sinh(3.849386397950811*y) - 
        1.0842021724855044e-19*sinh(5.674996469312898*y) + 
        1.0164395367051604e-20*sinh(10.161219234535377*y)))/2.;

			}

			if (fabs(this->constants.tau - 0.3) < 1e-5)
			{
				value[ID_theta] = -1.7443580752405021 + 0.39756135791713554*cosh(1.283128799316937*y) + 
   0.07237762337695752*cosh(1.8916654897709662*y) + 
   0.003071248679976172*cosh(3.3870730781784593*y) - 0.548685702383432*pow(y,2) + 
   0.13608276348795437*pow(y,4) - 2.7755575615628914e-17*sinh(1.283128799316937*y) - 
   3.469446951953614e-17*sinh(1.8916654897709662*y) - 
   1.5178830414797062e-18*sinh(3.3870730781784593*y);

				value[ID_stress] = -(pow(5,0.5)*(-0.07331283107342124 + 0.10266027610966066*cosh(1.283128799316937*y) - 
        0.03450433122665434*cosh(1.8916654897709662*y) + 
        0.004629281804058666*cosh(3.3870730781784593*y) - 
        0.1131370849898476*pow(y,2) - 
        4.85722573273506e-17*sinh(1.283128799316937*y) - 
        6.938893903907228e-18*sinh(1.8916654897709662*y) + 
        8.673617379884035e-19*sinh(3.3870730781784593*y)))/2.;
				
			}		
		}

		// The above values correspond to a unsymmetric system, therefore we now need to accomodate the symmetric system
		MatrixOpt::Base_MatrixOpt matrix_opt;
		value = matrix_opt.Sparse_matrix_dot_Vector(this->S_half,value);
	}

	// need to fix the following routines
	template<int dim>
	double
	G45_PoissonHeat<dim>::s_r(const double r,const double phi)const
	{
		return 0;

	}

	template<int dim>
	double
	G45_PoissonHeat<dim>::s_phi(const double r,const double phi)const
	{

				return 0;
	}

	template<int dim>
	double
	G45_PoissonHeat<dim>::thetaP(const double r,const double phi)const
	{
		return 0;

	}
}