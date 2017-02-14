namespace ExactSolution
{
	using namespace dealii;

	template<int dim>
	class
	ExactSolution_SystemA_ring: public Base_ExactSolution<dim>
	{
		public:
			ExactSolution_SystemA_ring(const constant_data &constants,const Sparse_matrix &S_half);

			virtual void vector_value(const Point<dim> &p,Vector<double> &value) const ;

			virtual double s_r(const double r,const double phi) const;
			virtual double s_phi(const double r,const double phi) const;
			virtual double thetaP(const double r,const double phi) const;			
	};

	template<int dim>
	ExactSolution_SystemA_ring<dim>::ExactSolution_SystemA_ring(const constant_data &constants,const Sparse_matrix &S_half)
	:
	Base_ExactSolution<dim>(constants,S_half)
	{
		Assert(this->constants.nEqn == 6,ExcMessage("incorrect equation number"));

		bool entered = false;

		Assert(fabs(this->constants.A0 - 0.0) < 1e-5 &&
			fabs(this->constants.A1 - 0.2) < 1e-5 &&
			fabs(this->constants.A2 - 0.1) < 1e-5 &&
			fabs(this->constants.uW - this->constants.tau) < 1e-5, ExcMessage("Please load the exact solution for this test problem"));

		// different test cases
		if (fabs(this->constants.A0 - 0.0) < 1e-5 &&
			fabs(this->constants.A1 - 0.2) < 1e-5 &&
			fabs(this->constants.A2 - 0.1) < 1e-5 &&
			fabs(this->constants.uW - this->constants.tau) < 1e-5)
		{

			if (fabs(this->constants.tau - 0.01) < 1e-5)
			{
					this->C1 = -258.6317250754188;
					this->C2 = 0.1403951449098079;
					this->C3 = -23.22459252429678;
					this->C4 = 6.490985249358112;
					this->K1 = 9.713664158637732e28;
					this->K2 = 1.0821414602962837e-127;					
					entered = true;
			}

				if( fabs(this->constants.tau - 0.1) < 1e-5 )
				{
					this->C1 = -1.9505992348859245;
					this->C2 = 0.143798955700584;
					this->C3 = 1.8448329709116922;
					this->C4 = 0.11342119104343287;
					this->K1 = 53.97166647969362;
					this->K2 = 1.26255082813026e-15;
					entered = true;
				}

				if( fabs(this->constants.tau - 1.0) < 1e-5 )
				{
					this->C1 = -0.14869287988905205;
					this->C2 = 0.17257205989623575;
					this->C3 = 1.2976971448738395;
					this->C4 = -0.10358631917431671;
					this->K1 = 0.18803597175199138;
					this->K2 = 0.003604064889912796;
					entered = true;

				}

				if( fabs(this->constants.tau - 10.0) < 1e-5 )
				{
					this->C1 = -0.1885054386029859;
					this->C2 = 1.4351918021177739;
					this->C3 = 0.11716270300600619;
					this->C4 = -0.004296145849336291;
					this->K1 = 0.18900239036185373;
					this->K2 = -0.6486989954482233;
					entered = true;
				}

				this->lambda1 = sqrt(2 * this->constants.zeta);
				
		}

		Assert(fabs(this->lambda1-sqrt(2 * this->constants.zeta)) < 1e-5 , ExcMessage("Parameters not initialized for the exact solution"));
		Assert(entered,ExcMessage("Parameters for exact solution not initialized"));

	}

	template<int dim>
	double
	ExactSolution_SystemA_ring<dim>::s_r(const double r,const double phi)const
	{
		return (this->constants.A0 * r * this->constants.tau)/2. - this->C4*pow(r,-1) + cos(phi) * 
				(-this->C2 + (2 * this->constants.A1 * r * this->constants.tau)/3. + this->C1*pow(r,-2) - 
			this->K2*this->BI(1,r*this->lambda1)*pow(2,0.5)*pow(r,-1) + this->K1*this->BK(1,r*this->lambda1)*pow(2,0.5)*pow(r,-1))
		+ (this->constants.A2 * pow(r,3) * pow(this->constants.tau,3))/4.;

	}

	template<int dim>
	double
	ExactSolution_SystemA_ring<dim>::s_phi(const double r,const double phi)const
	{

				return (this->C2 - (this->constants.A1 * r * this->constants.tau)/3. + this->C1*pow(r,-2) + 
					this->K2*(this->BI(0,this->lambda1 * r) * pow(this->constants.zeta,0.5) + 
						this->BI(2,this->lambda1 * r) * pow(this->constants.zeta,0.5)) + 
					this->K1*(this->BK(0,this->lambda1 * r) * pow(this->constants.zeta,0.5) 
						+ this->BK(2,this->lambda1 * r) * pow(this->constants.zeta,0.5)))*sin(phi);

	}

	template<int dim>
	double
	ExactSolution_SystemA_ring<dim>::thetaP(const double r,const double phi)const
	{
		return (this->C3 + this->C4 * this->constants.zeta * log(r) 
			- (this->constants.A0 * this->constants.tau * this->constants.zeta * pow(r,2))/4. + 
		cos(phi)*(this->C2 * r * this->constants.zeta + 
			this->C1 * this->constants.zeta * pow(r,-1) - (this->constants.A1 * this->constants.tau * (-2 + this->constants.zeta * pow(r,2)))/3.) + 
		(2 * this->constants.A2 * pow(r,2) * pow(this->constants.tau,3))/3.
		 - (this->constants.A2 * this->constants.zeta * pow(r,4) * pow(this->constants.tau,3))/16.);

	}

	template<int dim> 
	void 
	ExactSolution_SystemA_ring<dim>::
	vector_value(const Point<dim> &p,Vector<double> &value) const
	{
		// check whether value has been initialized or not
		Assert((int)value.size() == this->constants.nEqn,ExcNotInitialized());

		double r = sqrt(p.square());
		double phi = atan2(p[1],p[0]);
		r /= this->constants.tau;

		// we only provide the exact solution for theta, qx and qy
		value[0] =  thetaP(r,phi);		
		value[1] = cos(phi) * s_r(r,phi) - sin(phi) * s_phi(r,phi);
		value[2] = sin(phi) * s_r(r,phi) + cos(phi) * s_phi(r,phi);								

		// put all the other values to zero
		for (int i = 3 ; i < this->constants.nEqn ; i++)
			value[i] = 0; 


		// The above values correspond to a unsymmetric system, therefore we now need to accomodate the symmetric system
		MatrixOpt::Base_MatrixOpt matrix_opt;
		value = matrix_opt.Sparse_matrix_dot_Vector(this->S_half,value);

	}

}