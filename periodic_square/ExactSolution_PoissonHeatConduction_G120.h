// exact solution for poisson heat conduction
namespace ExactSolution
{
	using namespace dealii;

	template<int dim>
	class
	G120_PoissonHeat:public Base_ExactSolution<dim>
	{
		public:
			G120_PoissonHeat(const constant_data &constants,const Sparse_matrix &S_half);

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
	G120_PoissonHeat<dim>::G120_PoissonHeat(const constant_data &constants,const Sparse_matrix &S_half)
	:
	Base_ExactSolution<dim>(constants,S_half)
	{
		Assert(constants.nEqn == 70,ExcMessage("Only a solution for G120"));
	}

	template<int dim>
	void
	G120_PoissonHeat<dim>::vector_value(const Point<dim> &p,Vector<double> &value) const
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

				value[ID_theta] = -1.2995390567873166 + 0.008828936885305286*cosh(2.480359051432615*y) + 
   0.0008633188573549727*cosh(3.028788951922567*y) + 
   0.003055077322774268*cosh(4.434798678173429*y) + 
   0.0011530952298520565*cosh(7.524303732560254*y) + 
   0.000013607923731562518*cosh(12.976581617858354*y) - 0.18289523412781067*pow(y,2) + 
   0.408248290463863*pow(y,4);

				value[ID_stress] = -0.002715290039756343 + 0.0008187280616851396*cosh(2.480359051432615*y) - 
   0.0002604386379633425*cosh(3.028788951922567*y) + 
   0.0016721283891165937*cosh(4.434798678173429*y) - 
   0.00017655344915385124*cosh(7.524303732560254*y) + 
   0.00001945869320751019*cosh(12.976581617858354*y) - 0.03771236166328255*pow(y,2);

   				value[ID_heat] = -(pow(0.4,0.5)*pow(y,3))/3.;

			}

			if (fabs(this->constants.tau - 0.3) < 1e-5)
			{
				value[ID_theta] = -1.7436949415826675 + 0.2523069701617391*cosh(0.8267863504775383*y) + 
   0.02533033765112578*cosh(1.0095963173075224*y) + 
   0.11056926376618245*cosh(1.4782662260578097*y) + 
   0.07739730163732568*cosh(2.508101244186751*y) + 
   0.004179811421577409*cosh(4.325527205952785*y) - 0.548685702383432*pow(y,2) + 
   0.13608276348795437*pow(y,4);

				value[ID_stress] = -0.07331283107342124 + 0.02339701815900208*cosh(0.8267863504775383*y) - 
   0.007641439290718827*cosh(1.0095963173075224*y) + 
   0.060517618827160266*cosh(1.4782662260578097*y) - 
   0.011850504802646761*cosh(2.508101244186751*y) + 
   0.005976934448057989*cosh(4.325527205952785*y) - 0.1131370849898476*pow(y,2);

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
	G120_PoissonHeat<dim>::s_r(const double r,const double phi)const
	{
		return 0;

	}

	template<int dim>
	double
	G120_PoissonHeat<dim>::s_phi(const double r,const double phi)const
	{

				return 0;
	}

	template<int dim>
	double
	G120_PoissonHeat<dim>::thetaP(const double r,const double phi)const
	{
		return 0;

	}
}