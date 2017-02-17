// exact solution for poisson heat conduction
namespace ExactSolution
{
	using namespace dealii;

	template<int dim>
	class
	G105_PoissonHeat:public Base_ExactSolution<dim>
	{
		public:
			G105_PoissonHeat(const constant_data &constants,const Sparse_matrix &S_half);

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
	G105_PoissonHeat<dim>::G105_PoissonHeat(const constant_data &constants,const Sparse_matrix &S_half)
	:
	Base_ExactSolution<dim>(constants,S_half)
	{
		Assert(constants.nEqn == 62,ExcMessage("Only a solution for G105"));
	}

	template<int dim>
	void
	G105_PoissonHeat<dim>::vector_value(const Point<dim> &p,Vector<double> &value) const
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

				value[ID_theta] = -1.2997563950037463 + 0.009276947895287067*cosh(2.492986999351348*y) + 
   0.00034793324875276896*cosh(3.170853933340145*y) + 
   0.0019572547988298645*cosh(4.388700162109401*y) + 
   0.0017525338749780097*cosh(5.365272100745118*y) + 
   0.0006825478821284116*cosh(7.820516007402679*y) + 
   6.90472273800691e-7*cosh(10.44464164243648*y) + 
   7.880009718028203e-6*cosh(13.163405846932237*y) - 0.18289523412781067*pow(y,2) + 
   0.408248290463863*pow(y,4);

				value[ID_stress] = -0.002715290039756343 + 0.0008236716794642021*cosh(2.492986999351348*y) - 
   0.00023601868281006479*cosh(3.170853933340145*y) + 
   0.0011001213990571828*cosh(4.388700162109401*y) + 
   0.0006230403314441149*cosh(5.365272100745118*y) - 
   0.00020406994398414134*cosh(7.820516007402679*y) + 
   2.59443373639468e-6*cosh(10.44464164243648*y) + 
   0.000010993753807602042*cosh(13.163405846932237*y) - 0.03771236166328255*pow(y,2);

   				value[ID_heat] = -(pow(0.4,0.5)*pow(y,3))/3.;

			}

			if (fabs(this->constants.tau - 0.3) < 1e-5)
			{
				value[ID_theta] = -1.7441170080467425 + 0.2660463749885166*cosh(0.8309956664504494*y) + 
   0.00981560623079519*cosh(1.0569513111133817*y) + 
   0.07058859198806446*cosh(1.462900054036467*y) + 
   0.07231192665512148*cosh(1.7884240335817059*y) + 
   0.049617730356627546*cosh(2.6068386691342265*y) + 
   0.00009537413540004746*cosh(3.4815472141454937*y) + 
   0.0025738046900230706*cosh(4.387801948977413*y) - 0.548685702383432*pow(y,2) + 
   0.13608276348795437*pow(y,4);

				value[ID_stress] = -0.07331283107342124 + 0.023621439613073676*cosh(0.8309956664504494*y) - 
   0.006658364677359996*cosh(1.0569513111133817*y) + 
   0.03967598936112579*cosh(1.462900054036467*y) + 
   0.025707489820209463*cosh(1.7884240335817059*y) - 
   0.014834838286982257*cosh(2.6068386691342265*y) + 
   0.0003583661268530281*cosh(3.4815472141454937*y) + 
   0.0035908299765454493*cosh(4.387801948977413*y) - 0.1131370849898476*pow(y,2);

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
	G105_PoissonHeat<dim>::s_r(const double r,const double phi)const
	{
		return 0;

	}

	template<int dim>
	double
	G105_PoissonHeat<dim>::s_phi(const double r,const double phi)const
	{

				return 0;
	}

	template<int dim>
	double
	G105_PoissonHeat<dim>::thetaP(const double r,const double phi)const
	{
		return 0;

	}
}
