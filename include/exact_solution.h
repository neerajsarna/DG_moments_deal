namespace ExactSolution
{
	using namespace std;
	using namespace dealii;
	using namespace EquationGenerator;

	template<int dim> class Base_ExactSolution:public Function<dim>,
												protected Base_Basics
	{
		public:
			Base_ExactSolution(const unsigned int system_id,const unsigned int nEqn,
							  const system_matrix S_half,
							  const System_Type system_type,
							  physical_data &physical_constants,
							  string &output_dir);

			virtual void vector_value(const Point<dim> &p,Vector<double> &value) const;

		protected:
			const System_Type system_type;
			const system_matrix S_half;

			double BI(const int n,const double x) const; 
			double BK(const int n,const double x) const;
			double s_r(const double ,const double ) const;
			double s_phi(const double r,const double ) const;
			double thetaP(const double r,const double) const ;			
			double C0, C1, C2, C3, C4, K1, K2, K3, K4, K7, K8, 
        			lambda1, lambda2, lambda3;
			const unsigned int nEqn;
			const unsigned int system_id;
	};

	template<int dim> 
	Base_ExactSolution<dim>
	::Base_ExactSolution(const unsigned int system_id,
						 const unsigned int nEqn,
					     const system_matrix S_half,
					     const System_Type system_type,
					     physical_data &physical_constants,
					     string &output_dir)
	:
	Function<dim>(nEqn),
	Base_Basics(physical_constants,output_dir),
	system_type(system_type),
	S_half(S_half),
	nEqn(nEqn),
	system_id(system_id)
	{
		switch(nEqn)
		{
			case 6:
			{
				assert(nEqn == 6);
				if( fabs(this->tau - 0.01) < 1e-5 )
				{
					C1 = -258.6317250754188;
					C2 = 0.1403951449098079;
					C3 = -23.22459252429678;
					C4 = 6.490985249358112;
					K1 = 9.713664158637732e28;
					K2 = 1.0821414602962837e-127;	
				}


				if( fabs(this->tau - 0.1) < 1e-5 )
				{
					C1 = -1.9505992348859245;
					C2 = 0.143798955700584;
					C3 = 1.8448329709116922;
					C4 = 0.11342119104343287;
					K1 = 53.97166647969362;
					K2 = 1.26255082813026e-15;

				}

				if( fabs(this->tau - 1.0) < 1e-5 )
				{
					C1 = -0.14869287988905205;
					C2 = 0.17257205989623575;
					C3 = 1.2976971448738395;
					C4 = -0.10358631917431671;
					K1 = 0.18803597175199138;
					K2 = 0.003604064889912796;

				}

				if( fabs(this->tau - 10.0) < 1e-5 )
				{
					C1 = -0.1885054386029859;
					C2 = 1.4351918021177739;
					C3 = 0.11716270300600619;
					C4 = -0.004296145849336291;
					K1 = 0.18900239036185373;
					K2 = -0.6486989954482233;

				}

				lambda1 = sqrt(2 * this->zeta);
				break;

			}

			case 10:
			{
				assert(nEqn == 10);
				if( fabs(tau - 0.01) < 1e-5 )
				{
					C0 = -23.216704086856847;
					C1= 0.; 
					C2 = 0.;
					C3 = 6.489101153848549;
					K1 = 0.;
					K2 = 0.;
					K3 = -1.2300374620579039e18;
					K4 = -4.1431183805377506e-81;
					K7 = 0.;
					K8 = 0.;
				}


				if( fabs(tau - 0.1) < 1e-5 )
				{
					C0 = 1.8572099974915603;
					C1 = 0.;
					C2 = 0.;
					C3 = 0.10322985310279324;
					K1 = 0.;
					K2 = 0.;
					K3 = -0.8821623999826842;
					K4 = -3.124572666491822e-9;
					K7 = 0.;
					K8 = 0.;
				}

				if( fabs(tau - 1.0) < 1e-5 )
				{

					C0 = 1.8378924710991469;
					C1 = 0.;
					C2 = 0.;
					C3 = -0.1680773374037424;
					K1 = 0.;
					K2 = 0.;
					K3 = 0.13618194331758396;
					K4 = -0.4397130699425769;
					K7 = 0.;
					K8 = 0.;
				}

				if( fabs(tau - 10.0) < 1e-5 )
				{
					C0 = 642.4671720371373;
					C1 = 0.;
					C2 = 0.;
					C3 = -0.0279282371834904;
					K1 = 0.;
					K2 = 0.;
					K3 = 0.021731572883326996;
					K4 = -641.1814793927525;
					K7 = 0.;
					K8 = 0.;
				}

				lambda1 = sqrt(15*zeta/(15+16*zeta));
				lambda2 = sqrt(5.0/6.0);
				lambda3 = sqrt(3.0/2.0);
				break;
			}
		}
	}


template<int dim> 
double 
Base_ExactSolution<dim>::
s_r(const double r,const double phi) const
	{
		switch(nEqn)
		{
			case 6:
			{
				Assert(nEqn == 6,ExcNotImplemented());

				return (this->A0 * r * this->tau)/2. - C4*pow(r,-1) + cos(phi) * (-C2 + (2 * this->A1 * r * this->tau)/3. + C1*pow(r,-2) - 
						K2*BI(1,r*lambda1)*pow(2,0.5)*pow(r,-1) + K1*BK(1,r*lambda1)*pow(2,0.5)*pow(r,-1))
						+ (this->A2 * pow(r,3) * pow(this->tau,3))/4.;

				break;
			}

			case 10:
			{
				Assert(nEqn == 10,ExcNotImplemented());

		      return (A0*r*tau)/2. - C3*pow(r,-1) + (A2*pow(r,3)*pow(tau,3))/4. + 
      			cos(phi) * (K7*BI(1,lambda1*r)*pow(r,-1) + K8 * BK(1,lambda1*r)*pow(r,-1)
      			 - (2*A1*r*tau*(-9 + pow(r,2)))/27. - C2*pow(zeta,-1) + 
      			C1*pow(r,-2)*pow(zeta,-1));


				break;
			}
		}

		assert(1 == 0);
		return 0;

	}

template<int dim> 
double 
Base_ExactSolution<dim>::
s_phi(const double r,const double phi) const
	{
		switch(nEqn)
		{
			case 6:
			{
				Assert(nEqn == 6,ExcNotImplemented());

				return (C2 - (this->A1 * r * this->tau)/3. + C1*pow(r,-2) + 
					K2*(BI(0,lambda1 * r) * pow(this->zeta,0.5) + BI(2,lambda1 * r) * pow(this->zeta,0.5)) + 
					K1*(BK(0,lambda1 * r) * pow(this->zeta,0.5) + BK(2,lambda1 * r) * pow(this->zeta,0.5)))*sin(phi);

				break;
			}

			case 10:
			{
				Assert(nEqn == 10,ExcNotImplemented());

		      return  -((K7*lambda1*(this->BI(0,lambda1*r) + this->BI(2,lambda1*r)))/2. - 
      				(K8*lambda1*(this->BK(0,lambda1*r) + this->BK(2,lambda1*r)))/2. 
      				- (A1*r*tau*(-18 + pow(r,2)))/54. - C2*pow(zeta,-1) - 
			        C1*pow(r,-2)*pow(zeta,-1))*sin(phi);

				break;
			}
		}

		assert(1 == 0);
		return 0;
	}

template<int dim> 
double 
Base_ExactSolution<dim>::
thetaP(const double r,const double phi) const
	{
		switch(nEqn)
		{
			case 6:
			{
				assert(nEqn == 6);
				return C3 + C4 * this->zeta * log(r) - (this->A0 * this->tau * this->zeta * pow(r,2))/4. + 
				cos(phi)*(C2 * r * this->zeta + C1 * this->zeta * pow(r,-1) - (this->A1 * this->tau * (-2 + this->zeta * pow(r,2)))/3.) + 
				(2 * this->A2 * pow(r,2) * pow(this->tau,3))/3. - (this->A2 * this->zeta * pow(r,4) * pow(this->tau,3))/16.;

				break;
			}

			case 10:
			{
				assert(nEqn == 10);

      			return C0 + K4*this->BI(0,lambda2*r) + K3*this->BK(0,lambda2*r) + C3*zeta*log(r) - (A0*tau*zeta*pow(r,2))/4. + 
			      cos(phi)*(C2*r + K1*this->BI(1,lambda2*r) + K2*this->BK(1,lambda2*r) + C1*pow(r,-1) + 
      			(A1*tau*(-20 + zeta*(-18 + pow(r,2)))*pow(r,2))/54.) + (A2*pow(r,2)*(64 - 3*zeta*pow(r,2))*pow(tau,3))/48.;
				break;
			}
		}

		assert(1 == 0);
		return 0;

	}

template<int dim> 
void 
Base_ExactSolution<dim>::
vector_value(const Point<dim> &p,Vector<double> &value) const
	{

		double r = sqrt(p.square());
		double phi = atan2(p[1],p[0]);
		r /= this->tau;

		value[0] =  thetaP(r,phi);		
		value[1] = cos(phi) * s_r(r,phi) - sin(phi) * s_phi(r,phi);
		value[2] = sin(phi) * s_r(r,phi) + cos(phi) * s_phi(r,phi);								

		for (unsigned int i = 3 ; i < nEqn ; i++)
			value[i] = 0; 

		switch(system_type)
		{
			case un_symmetric:
				break;

			case symmetric:
			{
				Assert(nEqn != 10, ExcNotImplemented());
				Sparse_matrix_dot_Vector(S_half,value);
				break;
			}
		}

	}

 template<int dim> 
 double
 Base_ExactSolution<dim>::
 BI(const int n,const double x) const
	{
				return boost::math::cyl_bessel_i(n,x);
	}

 template<int dim> 
 double
 Base_ExactSolution<dim>
 ::BK(const int n,const double x) const
	{
		return boost::math::cyl_bessel_k(n,x);
	}

template<int dim> class systemA_period_sqr:public Base_ExactSolution<dim>
	{
		public:
			systemA_period_sqr(const unsigned int system_id,const unsigned int nEqn,
							  const system_matrix S_half,
							  const System_Type system_type,
							  physical_data &physical_constants,
							  string &output_dir);

			virtual void vector_value(const Point<dim> &p,Vector<double> &value) const;

	};

template<int dim> systemA_period_sqr<dim>
::systemA_period_sqr(const unsigned int system_id,
							  const unsigned int nEqn,
							  const system_matrix S_half,
							  const System_Type system_type,
							  physical_data &physical_constants,
							  string &output_dir)
:
Base_ExactSolution<dim>(system_id,
				   nEqn,
				   S_half,
				   system_type,
				   physical_constants,
				   output_dir)
{
	Assert(nEqn == 6,ExcMessage("number of equations larger or smaller than expected"));
}

template<int dim> 
void
systemA_period_sqr<dim>
::vector_value(const Point<dim> &p,Vector<double> &value) const
{
	// temprature different between the walls
	const double delta_theta_wall = this->theta0 - this->theta1; 			// thetaW(1)-thetaW(-1)

	// mean value of the temperature between the walls
	const double mean_theta_wall = (this->theta0 + this->theta1)/2;				// (thetaW(1) + thetaW(-1))/2
	const double ycord = p(1);
	
	Assert(value.size() == this->nEqn, ExcMessage("incorrect size"));
	// value of temperature
	value(0) = delta_theta_wall * ycord / (2 * (1 + this->tau)) + mean_theta_wall;

	// value of sx
	value(1) = 0;

	// value of sy
	value(2) = -delta_theta_wall * this->tau / (2 * (1 + this->tau));

	// value of Rxx
	value(3) = 0;

	// value of Rxy
	value(4) = 0;

	/// value of Ryy
	value(5) = 0;

	switch(this->system_type)
	{
		case un_symmetric:
		{
			Assert(1 == 0, ExcMessage("dont use unsymmetric system"));
			break;
		}

		case symmetric:
		{
			this->Sparse_matrix_dot_Vector(this->S_half,value);
			break;
		}
	}
	
	//cout << "exact solution is " << value << endl;

}

/*the exact solution of this particular system has not been implemented yet, 
so rather we will just have a zero function as the exact solution. So the error estimate will 
simply provide us with a norm of the solution*/
template<int dim> class G26_period_sqr: public Base_ExactSolution<dim>
{
	public:
		G26_period_sqr(const unsigned int system_id,const unsigned int nEqn,
							  const system_matrix S_half,
							  const System_Type system_type,
							  physical_data &physical_constants,
							  string &output_dir);

		virtual void vector_value(const Point<dim> &p,Vector<double> &value) const;
};

template<int dim> 
G26_period_sqr<dim>
::G26_period_sqr(const unsigned int system_id,
							  const unsigned int nEqn,
							  const system_matrix S_half,
							  const System_Type system_type,
							  physical_data &physical_constants,
							  string &output_dir)
:
Base_ExactSolution<dim>(system_id,
				   nEqn,
				   S_half,
				   system_type,
				   physical_constants,
				   output_dir)
{
	Assert(nEqn == 17, ExcMessage("not the desired number of equations"));
}

// presently the solution is not known
template<int dim>
void 
G26_period_sqr<dim>
::vector_value(const Point<dim> &p,Vector<double> &value) const
{
	Assert(value.size() == this->nEqn, ExcMessage("incorrect size"));
	value = 0;
}


}
