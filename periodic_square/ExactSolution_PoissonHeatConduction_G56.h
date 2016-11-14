// exact solution for poisson heat conduction
namespace ExactSolution
{
	using namespace dealii;

	template<int dim>
	class
	G56_PoissonHeat:public Base_ExactSolution<dim>
	{
		public:
			G56_PoissonHeat(const constant_data &constants,const Sparse_matrix &S_half);

			virtual void vector_value(const Point<dim> &p,Vector<double> &value) const ;

			virtual double s_r(const double r,const double phi) const;
			virtual double s_phi(const double r,const double phi) const;
			virtual double thetaP(const double r,const double phi) const;	
	};

	template<int dim>
	G56_PoissonHeat<dim>::G56_PoissonHeat(const constant_data &constants,const Sparse_matrix &S_half)
	:
	Base_ExactSolution<dim>(constants,S_half)
	{
		Assert(constants.nEqn == 34,ExcMessage("Only a solution for G56"));
	}

	template<int dim>
	void
	G56_PoissonHeat<dim>::vector_value(const Point<dim> &p,Vector<double> &value) const
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

				value[ID_theta] = -1.2991542242917795 + 0.008670378733644997*cosh(3.536982581771994*y) + 
   0.004430725803158755*cosh(5.204696760361044*y) + 
   0.000056930062888828204*cosh(10.45417499077846*y) - 0.18289523412781067*pow(y,2) + 
   0.408248290463863*pow(y,4) - 1.3010426069826053e-18*sinh(5.204696760361044*y) + 
   1.0164395367051604e-20*sinh(10.45417499077846*y);

				value[ID_stress] = -0.002715290039756343 + 0.0026403019195414515*cosh(3.536982581771994*y) - 
   0.0006407561824316862*cosh(5.204696760361044*y) + 
   0.00008292874116147389*cosh(10.45417499077846*y) - 0.03771236166328255*pow(y,2) - 
   3.469446951953614e-18*sinh(3.536982581771994*y) - 
   3.2526065174565133e-18*sinh(5.204696760361044*y);

   				value[ID_heat] = -(pow(0.4,0.5)*pow(y,3))/3.;

			}

			if (fabs(this->constants.tau - 0.3) < 1e-5)
			{
				value[ID_theta] = -1.7435149236809255 + 0.27723853580237484*cosh(1.178994193923998*y) + 
   0.18476156620816786*cosh(1.734898920120348*y) + 
   0.008511751878547334*cosh(3.4847249969261536*y) - 0.548685702383432*pow(y,2) + 
   0.13608276348795437*pow(y,4) - 5.551115123125783e-17*sinh(1.734898920120348*y) + 
   8.673617379884035e-19*sinh(3.4847249969261536*y);

				value[ID_stress] = -(pow(5,0.5)*(-0.07331283107342124 + 0.0844246209694861*cosh(1.178994193923998*y) - 
        0.02671957622366166*cosh(1.734898920120348*y) + 
        0.0123988773689772*cosh(3.4847249969261536*y) - 0.1131370849898476*pow(y,2) - 
        1.1102230246251565e-16*sinh(1.178994193923998*y) - 
        1.3357370765021415e-16*sinh(1.734898920120348*y)))/2.;
				
			}		
		}

		// The above values correspond to a unsymmetric system, therefore we now need to accomodate the symmetric system
		MatrixOpt::Base_MatrixOpt matrix_opt;
		value = matrix_opt.Sparse_matrix_dot_Vector(this->S_half,value);
	}

	// need to fix the following routines
	template<int dim>
	double
	G56_PoissonHeat<dim>::s_r(const double r,const double phi)const
	{
		return 0;

	}

	template<int dim>
	double
	G56_PoissonHeat<dim>::s_phi(const double r,const double phi)const
	{

				return 0;
	}

	template<int dim>
	double
	G56_PoissonHeat<dim>::thetaP(const double r,const double phi)const
	{
		return 0;

	}
}