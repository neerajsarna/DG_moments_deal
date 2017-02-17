// exact solution for poisson heat conduction
namespace ExactSolution
{
	using namespace dealii;

	template<int dim>
	class
	G84_PoissonHeat:public Base_ExactSolution<dim>
	{
		public:
			G84_PoissonHeat(const constant_data &constants,const Sparse_matrix &S_half);

			virtual void vector_value(const Point<dim> &p,Vector<double> &value) const ;

			virtual double s_r(const double r,const double phi) const;
			virtual double s_phi(const double r,const double phi) const;
			virtual double thetaP(const double r,const double phi) const;	

			virtual double R_rr(const double r,const double phi)const {return 0;};
			virtual double R_thetatheta(const double r,const double phi)const {return 0;};
			virtual double R_rtheta(const double r,const double phi)const {return 0;};
			virtual double R_zz(const double r,const double phi)const {return 0;};
	};

	template<int dim>
	G84_PoissonHeat<dim>::G84_PoissonHeat(const constant_data &constants,const Sparse_matrix &S_half)
	:
	Base_ExactSolution<dim>(constants,S_half)
	{
		Assert(constants.nEqn == 50,ExcMessage("Only a solution for G84"));
	}

	template<int dim>
	void
	G84_PoissonHeat<dim>::vector_value(const Point<dim> &p,Vector<double> &value) const
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

				value[ID_theta] = -1.3016616131237297 + 0.012498585238508764*cosh(2.867124849692578*y) + 
   0.0026996208249191254*cosh(3.6131625584441966*y) + 
   0.0008429866277644128*cosh(6.002769036811866*y) - 0.18289523412781067*pow(y,2) + 
   0.408248290463863*pow(y,4) - 2.6020852139652106e-18*sinh(2.867124849692578*y) - 
   8.673617379884035e-19*sinh(3.6131625584441966*y) + 
   1.0842021724855044e-19*sinh(6.002769036811866*y);

				value[ID_stress] = -0.002715290039756343 + 0.0018033048153696892*cosh(2.867124849692578*y) - 
   0.0007033685726265145*cosh(3.6131625584441966*y) + 
   0.0008142812934811178*cosh(6.002769036811866*y) - 0.03771236166328255*pow(y,2) + 
   1.951563910473908e-18*sinh(2.867124849692578*y) - 
   1.6263032587282567e-18*sinh(3.6131625584441966*y) + 
   1.0842021724855044e-19*sinh(6.002769036811866*y);

   				value[ID_heat] = -(pow(0.4,0.5)*pow(y,3))/3.;

			}

			if (fabs(this->constants.tau - 0.3) < 1e-5)
			{
				value[ID_theta] = -1.7479896311184409 + 0.36064751628260816*cosh(0.9557082832308594*y) + 
   0.07833003503050456*cosh(1.2043875194813989*y) + 
   0.04090384023127845*cosh(2.000923012270622*y) - 0.548685702383432*pow(y,2) + 
   0.13608276348795437*pow(y,4) - 8.326672684688674e-17*sinh(0.9557082832308594*y) - 
   3.469446951953614e-17*sinh(1.2043875194813989*y) + 
   3.469446951953614e-18*sinh(2.000923012270622*y);

				value[ID_stress] = -0.07331283107342124 + 0.05203448153153864*cosh(0.9557082832308594*y) - 
   0.0204083789933134*cosh(1.2043875194813989*y) + 
   0.0395109849134863*cosh(2.000923012270622*y) - 0.1131370849898476*pow(y,2) + 
   5.551115123125783e-17*sinh(0.9557082832308594*y) - 
   4.683753385137379e-17*sinh(1.2043875194813989*y) + 
   6.938893903907228e-18*sinh(2.000923012270622*y);

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
	G84_PoissonHeat<dim>::s_r(const double r,const double phi)const
	{
		return 0;

	}

	template<int dim>
	double
	G84_PoissonHeat<dim>::s_phi(const double r,const double phi)const
	{

				return 0;
	}

	template<int dim>
	double
	G84_PoissonHeat<dim>::thetaP(const double r,const double phi)const
	{
		return 0;

	}
}
