// exact solution for poisson heat conduction
namespace ExactSolution
{
	using namespace dealii;

	template<int dim>
	class
	G26_PoissonHeat:public Base_ExactSolution<dim>
	{
		public:
			G26_PoissonHeat(const constant_data &constants,const Sparse_matrix &S_half);

			virtual void vector_value(const Point<dim> &p,Vector<double> &value) const ;

			virtual double s_r(const double r,const double phi) const;
			virtual double s_phi(const double r,const double phi) const;
			virtual double thetaP(const double r,const double phi) const;	

			virtual double R_rr(const double r,const double phi)const {return 0;};
			virtual double R_thetatheta(const double r,const double phi)const {return 0;};
			virtual double R_rtheta(const double r,const double phi)const {return 0;};
			virtual double R_zz(const double r,const double phi)const {return 0;};
	};

	template<>
	G26_PoissonHeat<2>::G26_PoissonHeat(const constant_data &constants,const Sparse_matrix &S_half)
	:
	Base_ExactSolution<2>(constants,S_half)
	{
		Assert(constants.nEqn == 17,ExcMessage("Only a solution for G26"));
	}

	template<>
	G26_PoissonHeat<1>::G26_PoissonHeat(const constant_data &constants,const Sparse_matrix &S_half)
	:
	Base_ExactSolution<1>(constants,S_half)
	{
		Assert(constants.nEqn == 8,ExcMessage("Only a solution for G26"));
	}


	template<>
	void
	G26_PoissonHeat<2>::vector_value(const Point<2> &p,Vector<double> &value) const
	{
		// first we check the size of the value vector
		Assert((int)value.size() == this->constants.nEqn,ExcNotInitialized());
		Assert(fabs(this->constants.alpha - 0.816496580927726) < 1e-5,ExcMessage("Incorrect alpha value"));

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

				value[ID_theta] = -1.2852476333586613 - 0.18289523412781067*pow(y,2) + 0.408248290463863*pow(y,4) + 
				0.0015643311536178803*cosh(6.573421981221795 * y) - 
				2.168404344971009e-19*sinh(6.573421981221795 * y);

				value[ID_stress] = -0.002715290039756343 - 0.03771236166328255*pow(y,2) + 
				0.0011289587658037511*cosh(6.573421981221795*y);

				value[ID_heat] = -0.21081851067789195*pow(y,3);

			}

			if (fabs(this->constants.tau - 0.3) < 1e-5)
			{
				value[ID_theta] = -1.3650564316627332 - 0.548685702383432*pow(y,2) + 0.13608276348795437*pow(y,4) + 
				0.09338201417693054*cosh(2.191140660407265*y);

				value[ID_stress] = -0.07331283107342124 - 0.1131370849898476*pow(y,2) + 
				0.06739266377815037*cosh(2.191140660407265*y);

				value[ID_heat] = -0.21081851067789195*pow(y,3);
			}

			if (fabs(this->constants.tau - 0.5) < 1e-5)
			{
				value[ID_theta] = -1.7298151776227189 - 0.9144761706390533*pow(y,2) + 0.08164965809277261*pow(y,4) + 
				0.45942750968124524*cosh(1.314684396244359*y) - 
				2.7755575615628914e-17*sinh(1.314684396244359*y);

				value[ID_stress] = -0.3394112549695428 - 0.1885618083164127*pow(y,2) + 
				0.33156324548448296*cosh(1.314684396244359*y);

				value[ID_heat] = -0.21081851067789195*pow(y,3);
			}			
		}

		// The above values correspond to a unsymmetric system, therefore we now need to accomodate the symmetric system
		MatrixOpt::Base_MatrixOpt matrix_opt;
		value = matrix_opt.Sparse_matrix_dot_Vector(this->S_half,value);
	}

	template<>
	void
	G26_PoissonHeat<1>::vector_value(const Point<1> &p,Vector<double> &value) const
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
				value[ID_theta] = 0.9918095457622201*x + 0.004077339347408366*sinh(2.517180972479634*x) + 
   								  0.0031858247584759585*sinh(6.71508535912671*x);
			
		

		// The above values correspond to a unsymmetric system, therefore we now need to accomodate the symmetric system
		MatrixOpt::Base_MatrixOpt matrix_opt;
		value = matrix_opt.Sparse_matrix_dot_Vector(this->S_half,value);
	}


	// need to fix the following routines
	template<int dim>
	double
	G26_PoissonHeat<dim>::s_r(const double r,const double phi)const
	{
		Assert(1 == 0, ExcMessage("Should not have reached here"));
		return 0;

	}

	template<int dim>
	double
	G26_PoissonHeat<dim>::s_phi(const double r,const double phi)const
	{
		Assert(1 == 0, ExcMessage("Should not have reached here"));
		return 0;

	}

	template<int dim>
	double
	G26_PoissonHeat<dim>::thetaP(const double r,const double phi)const
	{
		Assert(1 == 0, ExcMessage("Should not have reached here"));
		return 0;

	}
}
