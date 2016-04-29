namespace ExactSolution
{
	using namespace std;
	using namespace dealii;
	using namespace Basics;

	template<int dim> class Base_ExactSolution:public Function<dim>
	{
		public:
			Base_ExactSolution(const unsigned int nEqn);
			virtual void vector_value(const Point<dim> &p,Vector<double> &value) const = 0;

		protected:
			double BI(const int n,const double x) const; 
			double BK(const int n,const double x) const;
			virtual double s_r(const double ,const double ) const = 0;
			virtual double s_phi(const double r,const double ) const = 0;
			virtual double thetaP(const double r,const double) const = 0;			
	};

	template<int dim> Base_ExactSolution<dim>::Base_ExactSolution(const unsigned int nEqn)
	:
	Function<dim>(nEqn)
	{}


	template<int dim> double Base_ExactSolution<dim>::BI(const int n,const double x) const
	{
		return boost::math::cyl_bessel_i(n,x);
	}

	template<int dim> double Base_ExactSolution<dim>::BK(const int n,const double x) const
	{
		return boost::math::cyl_bessel_k(n,x);
	}

	template<int dim> class exact_solutionA:public Base_ExactSolution<dim>,protected constants_SystemA
	{
		public:
			exact_solutionA(const unsigned int nEqn);
			virtual void vector_value(const Point<dim> &p,Vector<double> &value) const;

		private:
			double C1, C2, C3, C4, K1, K2, lambda1;
			virtual double s_r(const double ,const double ) const;
			virtual double s_phi(const double r,const double ) const;
			virtual double thetaP(const double r,const double) const;		
	};

template<int dim> exact_solutionA<dim>::exact_solutionA(const unsigned int nEqn)
	:
	Base_ExactSolution<dim>(nEqn)
	{

		if( fabs(tau - 0.01) < 1e-5 )
		{
			C1 = -258.6317250754188;
			C2 = 0.1403951449098079;
			C3 = -23.22459252429678;
			C4 = 6.490985249358112;
			K1 = 9.713664158637732e28;
			K2 = 1.0821414602962837e-127;	
		}


		if( fabs(tau - 0.1) < 1e-5 )
		{
			C1 = -1.9505992348859245;
			C2 = 0.143798955700584;
			C3 = 1.8448329709116922;
			C4 = 0.11342119104343287;
			K1 = 53.97166647969362;
			K2 = 1.26255082813026e-15;

		}

		if( fabs(tau - 1.0) < 1e-5 )
		{
			C1 = -0.14869287988905205;
			C2 = 0.17257205989623575;
			C3 = 1.2976971448738395;
			C4 = -0.10358631917431671;
			K1 = 0.18803597175199138;
			K2 = 0.003604064889912796;

		}

		if( fabs(tau - 10.0) < 1e-5 )
		{
			C1 = -0.1885054386029859;
			C2 = 1.4351918021177739;
			C3 = 0.11716270300600619;
			C4 = -0.004296145849336291;
			K1 = 0.18900239036185373;
			K2 = -0.6486989954482233;

		}

		lambda1 = sqrt(2 * zeta);

	}



template<int dim> double exact_solutionA<dim>::thetaP(const double r,const double phi) const
	{
		return C3 + C4 * zeta * log(r) - (A0 * tau * zeta * pow(r,2))/4. + 
		cos(phi)*(C2 * r * zeta + C1 * zeta * pow(r,-1) - (A1 * tau * (-2 + zeta * pow(r,2)))/3.) + 
		(2 * A2 * pow(r,2) * pow(tau,3))/3. - (A2 * zeta * pow(r,4) * pow(tau,3))/16.;
	}

template<int dim> double exact_solutionA<dim>::s_r(const double r,const double phi) const
	{
		return (A0 * r * tau)/2. - C4*pow(r,-1) + cos(phi) * (-C2 + (2 * A1 * r * tau)/3. + C1*pow(r,-2) - 
			K2*this->BI(1,r*lambda1)*pow(2,0.5)*pow(r,-1) + K1*this->BK(1,r*lambda1)*pow(2,0.5)*pow(r,-1))
		+ (A2 * pow(r,3) * pow(tau,3))/4.;

	}

template<int dim> double exact_solutionA<dim>::s_phi(const double r,const double phi) const
	{
		return (C2 - (A1 * r * tau)/3. + C1*pow(r,-2) + 
			K2*(this->BI(0,lambda1 * r) * pow(zeta,0.5) + this->BI(2,lambda1 * r) * pow(zeta,0.5)) + 
			K1*(this->BK(0,lambda1 * r) * pow(zeta,0.5) + this->BK(2,lambda1 * r) * pow(zeta,0.5)))*sin(phi);

	}


 template<int dim> void exact_solutionA<dim>::vector_value(const Point<dim> &p,Vector<double> &value) const
	{
		assert(value.size() != 0);

		double r = sqrt(p.square());
		double phi = atan2(p[1],p[0]);
		r /=tau;



		value[0] =  thetaP(r,phi);
											// theta
		value[1] = cos(phi) * s_r(r,phi) - sin(phi) * s_phi(r,phi);
								// qx
		value[2] = sin(phi) * s_r(r,phi) + cos(phi) * s_phi(r,phi);								// qy

	// allocating zero value for higher order moments(values not needed during error evaluation)
		for (unsigned int i = 3 ; i < this->nEqn ; i++)
			value[i] = 0; // a value for higher order moments is not needed

   }

	template<int dim> class exact_solutionB:public Base_ExactSolution<dim>,protected constants_SystemB
	{
		public:
			exact_solutionB(const unsigned int nEqn);
			virtual void vector_value(const Point<dim> &p,Vector<double> &value) const;

		private:       
			double C0, C1, C2, C3, K1, K2, K3, K4, K7, K8, 
        			lambda1, lambda2, lambda3;
			virtual double s_r(const double ,const double ) const;
			virtual double s_phi(const double r,const double ) const;
			virtual double thetaP(const double r,const double) const;
		
	};


	template<int dim> exact_solutionB<dim>::exact_solutionB(const unsigned int nEqn)
	:
	Base_ExactSolution<dim>(nEqn)
	{

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
	}

	template<int dim> double exact_solutionB<dim>::thetaP(const double r,const double phi) const
	{
      return C0 + K4*this->BI(0,lambda2*r) + K3*this->BK(0,lambda2*r) + C3*zeta*log(r) - (A0*tau*zeta*pow(r,2))/4. + 
      cos(phi)*(C2*r + K1*this->BI(1,lambda2*r) + K2*this->BK(1,lambda2*r) + C1*pow(r,-1) + 
      (A1*tau*(-20 + zeta*(-18 + pow(r,2)))*pow(r,2))/54.) + (A2*pow(r,2)*(64 - 3*zeta*pow(r,2))*pow(tau,3))/48.;
	}

	template<int dim> double exact_solutionB<dim>::s_r(const double r,const double phi) const
	{
      return (A0*r*tau)/2. - C3*pow(r,-1) + (A2*pow(r,3)*pow(tau,3))/4. + 
      cos(phi)*(K7*this->BI(1,lambda1*r)*pow(r,-1) + K8*this->BK(1,lambda1*r)*pow(r,-1) - (2*A1*r*tau*(-9 + pow(r,2)))/27. - C2*pow(zeta,-1) + 
      C1*pow(r,-2)*pow(zeta,-1));

	}

	template<int dim> double exact_solutionB<dim>::s_phi(const double r,const double phi) const
	{
      return  -((K7*lambda1*(this->BI(0,lambda1*r) + this->BI(2,lambda1*r)))/2. - 
      (K8*lambda1*(this->BK(0,lambda1*r) + this->BK(2,lambda1*r)))/2. - (A1*r*tau*(-18 + pow(r,2)))/54. - C2*pow(zeta,-1) - 
      C1*pow(r,-2)*pow(zeta,-1))*sin(phi);

	}

	template<int dim> void exact_solutionB<dim>::vector_value(const Point<dim> &p,Vector<double> &value) const
	{
		assert(value.size() == this->nEqn);

		double r = sqrt(p.square());
		double phi = atan2(p[1],p[0]);
		r /=tau;



		value[0] =  thetaP(r,phi);
											// theta
		value[1] = cos(phi) * s_r(r,phi) - sin(phi) * s_phi(r,phi);
								// qx
		value[2] = sin(phi) * s_r(r,phi) + cos(phi) * s_phi(r,phi);								// qy

		// allocating zero value for higher order moments(values not needed during error evaluation)
		for (unsigned int i = 3 ; i < this->nEqn ; i++)
			value[i] = 0; // a value for higher order moments is not needed

   }

}
