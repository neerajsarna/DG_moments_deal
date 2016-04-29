namespace ExactSolution
{
	using namespace std;
	using namespace dealii;
	using namespace EquationGenerator;

	template<int system_type,int num_flux,int dim> class Base_ExactSolution:public Function<dim>,
																			protected Base_Basics
	{
		public:
			Base_ExactSolution(const unsigned int system_id,const unsigned int nEqn);
			virtual void vector_value(const Point<dim> &p,Vector<double> &value) const;

		protected:
			double BI(const int n,const double x) const; 
			double BK(const int n,const double x) const;
			double s_r(const double ,const double ) const;
			double s_phi(const double r,const double ) const;
			double thetaP(const double r,const double) const ;			
			double C1, C2, C3, C4, K1, K2, lambda1;
			const unsigned int nEqn;
			const unsigned int system_id;
	};

	template<int system_type,int num_flux,int dim> 
	Base_ExactSolution<system_type,num_flux,dim>
	::Base_ExactSolution(const unsigned int system_id,const unsigned int nEqn)
	:
	nEqn(nEqn),
	system_id(system_id),
	Function<dim>(nEqn)
	{
		switch(system_id)
		{
			case 0:
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
		}
	}


template<int system_type,int num_flux,int dim> 
double 
Base_ExactSolution<system_type,num_flux,dim>::
s_r(const double r,const double phi) const
	{
		switch(system_id)
		{
			case 0:
			{
				assert(nEqn == 6);

				return (this->A0 * r * this->tau)/2. - C4*pow(r,-1) + cos(phi) * (-C2 + (2 * this->A1 * r * this->tau)/3. + C1*pow(r,-2) - 
						K2*BI(1,r*lambda1)*pow(2,0.5)*pow(r,-1) + K1*BK(1,r*lambda1)*pow(2,0.5)*pow(r,-1))
						+ (this->A2 * pow(r,3) * pow(this->tau,3))/4.;

				break;
			}
		}

		assert(1 == 0);
		return 0;

	}

template<int system_type,int num_flux,int dim> 
double 
Base_ExactSolution<system_type,num_flux,dim>::
s_phi(const double r,const double phi) const
	{
		switch(system_id)
		{
			case 0:
			{
				assert(nEqn == 6);

				return (C2 - (this->A1 * r * this->tau)/3. + C1*pow(r,-2) + 
					K2*(BI(0,lambda1 * r) * pow(this->zeta,0.5) + BI(2,lambda1 * r) * pow(this->zeta,0.5)) + 
					K1*(BK(0,lambda1 * r) * pow(this->zeta,0.5) + BK(2,lambda1 * r) * pow(this->zeta,0.5)))*sin(phi);

				break;
			}
		}

		assert(1 == 0);
		return 0;
	}

template<int system_type,int num_flux,int dim> 
double 
Base_ExactSolution<system_type,num_flux,dim>::
thetaP(const double r,const double phi) const
	{
		switch(system_id)
		{
			case 0:
			{
				assert(nEqn == 6);
				return C3 + C4 * this->zeta * log(r) - (this->A0 * this->tau * this->zeta * pow(r,2))/4. + 
				cos(phi)*(C2 * r * this->zeta + C1 * this->zeta * pow(r,-1) - (this->A1 * this->tau * (-2 + this->zeta * pow(r,2)))/3.) + 
				(2 * this->A2 * pow(r,2) * pow(this->tau,3))/3. - (this->A2 * this->zeta * pow(r,4) * pow(this->tau,3))/16.;

				break;
			}
		}

		assert(1 == 0);
		return 0;

	}

template<int system_type,int num_flux,int dim> 
void 
Base_ExactSolution<system_type,num_flux,dim>::
vector_value(const Point<dim> &p,Vector<double> &value) const
	{
		switch(system_id)
		{
			case 0:
			{
				double r = sqrt(p.square());
				double phi = atan2(p[1],p[0]);
				r /= this->tau;

				value[0] =  thetaP(r,phi);		
				value[1] = cos(phi) * s_r(r,phi) - sin(phi) * s_phi(r,phi);
				value[2] = sin(phi) * s_r(r,phi) + cos(phi) * s_phi(r,phi);								

				for (unsigned int i = 3 ; i < nEqn ; i++)
					value[i] = 0; 

				break;
			}

		}
	}

 template<int system_type,int num_flux,int dim> 
 double
 Base_ExactSolution<system_type,num_flux,dim>::
 BI(const int n,const double x) const
	{
				return boost::math::cyl_bessel_i(n,x);
	}

 template<int system_type,int num_flux,int dim> 
 double
 Base_ExactSolution<system_type,num_flux,dim>::BK(const int n,const double x) const
	{
		return boost::math::cyl_bessel_k(n,x);
	}



}
