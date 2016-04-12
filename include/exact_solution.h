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
			virtual double s_r(const double ,const double ) const = 0;
			virtual double s_phi(const double r,const double ) const = 0;
			virtual double thetaP(const double r,const double) const = 0;
			virtual double BI(const int n,const double x) const = 0;
			virtual double BK(const int n,const double x) const = 0;			
	};

	template<int dim> Base_ExactSolution<dim>::Base_ExactSolution(const unsigned int nEqn)
	:
	Function<dim>(nEqn)
	{;};

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
			virtual double BI(const int n,const double x) const;
			virtual double BK(const int n,const double x) const;			
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

template<int dim> double exact_solutionA<dim>::BI(const int n,const double x) const
	{
		return boost::math::cyl_bessel_i(n,x);
	}

template<int dim> double exact_solutionA<dim>::BK(const int n,const double x) const
	{
		return boost::math::cyl_bessel_k(n,x);
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
			K2*BI(1,r*lambda1)*pow(2,0.5)*pow(r,-1) + K1*BK(1,r*lambda1)*pow(2,0.5)*pow(r,-1))
		+ (A2 * pow(r,3) * pow(tau,3))/4.;

	}

template<int dim> double exact_solutionA<dim>::s_phi(const double r,const double phi) const
	{
		return (C2 - (A1 * r * tau)/3. + C1*pow(r,-2) + 
			K2*(BI(0,lambda1 * r) * pow(zeta,0.5) + BI(2,lambda1 * r) * pow(zeta,0.5)) + 
			K1*(BK(0,lambda1 * r) * pow(zeta,0.5) + BK(2,lambda1 * r) * pow(zeta,0.5)))*sin(phi);

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
	value[3] = 0;								//
	value[4] = 0;
	value[5] = 0;

}

}
