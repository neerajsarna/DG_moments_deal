// exact solution for poisson heat conduction
namespace ExactSolution
{
	using namespace dealii;

	template<int dim>
	class
	G71_PoissonHeat:public Base_ExactSolution<dim>
	{
		public:
			G71_PoissonHeat(const constant_data &constants,const Sparse_matrix &S_half);

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
	G71_PoissonHeat<dim>::G71_PoissonHeat(const constant_data &constants,const Sparse_matrix &S_half)
	:
	Base_ExactSolution<dim>(constants,S_half)
	{
		Assert(constants.nEqn == 43,ExcMessage("Only a solution for G71"));
	}

	template<int dim>
	void
	G71_PoissonHeat<dim>::vector_value(const Point<dim> &p,Vector<double> &value) const
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

				value[ID_theta] = -1.3012568770728508 + 0.014368644008044136*cosh(2.9184816217908978*y) + 
   0.000847741867453309*cosh(3.885676724821849*y) + 
   0.0004461674640337943*cosh(5.896124120768356*y) + 
   0.00014238973643247108*cosh(9.591077444821629*y) - 0.18289523412781067*pow(y,2) + 
   0.408248290463863*pow(y,4) + 8.673617379884035e-19*sinh(3.885676724821849*y) - 
   1.3552527156068805e-20*sinh(9.591077444821629*y);

				value[ID_stress] = -0.002715290039756343 + 0.0018696658823278787*cosh(2.9184816217908978*y) - 
   									0.0006185894211626305*cosh(3.885676724821849*y) + 
   									0.0004424972314741653*cosh(5.896124120768356*y) + 
   0.00010260159078614521*cosh(9.591077444821629*y) - 0.03771236166328255*pow(y,2) - 
   4.445228907190568e-18*sinh(2.9184816217908978*y) + 
   3.7947076036992655e-19*sinh(3.885676724821849*y) + 
   2.439454888092385e-19*sinh(5.896124120768356*y) - 
   4.0657581468206416e-20*sinh(9.591077444821629*y);

   				value[ID_heat] = -(pow(0.4,0.5)*pow(y,3))/3.;

			}

			if (fabs(this->constants.tau - 0.3) < 1e-5)
			{
				value[ID_theta] = -1.747245655093501 + 0.41369103486028996*cosh(0.9728272072636326*y) + 
   0.02470452567666627*cosh(1.2952255749406163*y) + 
   0.021023428196056175*cosh(1.9653747069227854*y) + 
   0.017552620592354434*cosh(3.197025814940543*y) - 0.548685702383432*pow(y,2) + 
   0.13608276348795437*pow(y,4) + 2.6020852139652106e-17*sinh(1.2952255749406163*y) - 
   1.734723475976807e-18*sinh(3.197025814940543*y);

				value[ID_stress] = -0.07331283107342124 + 0.053829993510186594*cosh(0.9728272072636326*y) - 
   0.01802666451326118*cosh(1.2952255749406163*y) + 
   0.020850486695610182*cosh(1.9653747069227854*y) + 
   0.012647869434713875*cosh(3.197025814940543*y) - 0.1131370849898476*pow(y,2) - 
   1.249000902703301e-16*sinh(0.9728272072636326*y) + 
   1.214306433183765e-17*sinh(1.2952255749406163*y) + 
   1.214306433183765e-17*sinh(1.9653747069227854*y) - 
   6.071532165918825e-18*sinh(3.197025814940543*y);

   				value[ID_heat] = -(pow(0.4,0.5)*pow(y,3))/3.;
				
			}		
		}

		// The above values correspond to a unsymmetric system, therefore we now need to accomodate the symmetric system
		MatrixOpt::Base_MatrixOpt matrix_opt;
		value = matrix_opt.Sparse_matrix_dot_Vector(this->S_half,value);
	}

	// need to fix the following routines
	template<int dim>
	double
	G71_PoissonHeat<dim>::s_r(const double r,const double phi)const
	{
		return 0;

	}

	template<int dim>
	double
	G71_PoissonHeat<dim>::s_phi(const double r,const double phi)const
	{

				return 0;
	}

	template<int dim>
	double
	G71_PoissonHeat<dim>::thetaP(const double r,const double phi)const
	{
		return 0;

	}
}