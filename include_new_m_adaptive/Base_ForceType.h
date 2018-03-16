// in the following namespace we declare the different forces we can have in the system
namespace ForceType
{
	using namespace dealii;

	template<int dim>
	class
	Base_ForceType
	{
		public:
			Base_ForceType(const constant_numerics &constants,const int nEqn);

			const constant_numerics constants;
			const int nEqn;

		public:
		// The source term
		// The term factor appears due to symmetrization
			virtual void source_term(const std::vector<Point<dim>> &p,
									 std::vector<Vector<double>> &value,
									 const double factor) = 0;
	};

	template<int dim>
	Base_ForceType<dim>::
	Base_ForceType(const constant_numerics &constants,const int nEqn)
	:
	constants(constants),
	nEqn(nEqn)
	{;}

	template<int dim>
	class 
	ForceType1:public Base_ForceType<dim>
	{
		public:
			ForceType1(const constant_numerics &constants,const int nEqn);

			// the force will be acting on the first equation
			const unsigned int var_rho = 0;

			// The source term
			virtual void source_term(const std::vector<Point<dim>> &p,
									 std::vector<Vector<double>> &value,
									 const double factor);
	};

	template<int dim>
	ForceType1<dim>
	::ForceType1(const constant_numerics &constants,const int nEqn)
	:
	Base_ForceType<dim>(constants,nEqn)
	{;}

	// since the forcing only depends upon the x-coordinate so we do not need to specialize the following
	// force type.
	template<int dim>
	void 
	ForceType1<dim>
	::source_term(const std::vector<Point<dim>> &p,
				  std::vector<Vector<double>> &value,
			      const double factor)
	{

		Assert(p.size()!=0,ExcNotInitialized());
		Assert(value.size()!=0,ExcNotInitialized());
		AssertDimension(p.size(),value.size());


		for (unsigned int i = 0 ; i < value.size() ; i++)
			{
				const double norm = sqrt(p[i].square());
				const double x_coord = p[i][0];

				value[i] = 0;
				AssertDimension((int)value[i].size(),this->nEqn);

				// the factor comes from the symmetrizer
				value[i](var_rho) = factor * (this->constants.A0 + this->constants.A2*norm*norm 
									+ this->constants.A1*x_coord/norm);	
			}

	}


	// Force corresponding to exact solution of systemB
	template<int dim>
	class 
	ForceType2:public Base_ForceType<dim>
	{
		public:
			ForceType2(const constant_numerics &constants,const int nEqn);

			const unsigned int var_rho = 0;
		// The source term
		virtual void source_term(const std::vector<Point<dim>> &p,
								 std::vector<Vector<double>> &value,
								 const double factor);
	};

	template<int dim>
	ForceType2<dim>
	::ForceType2(const constant_numerics &constants,const int nEqn)
	:
	Base_ForceType<dim>(constants,nEqn)
	{;}

	// since we only have the x-coordinate in the force definition therefore 
	// we do not need to specialize the template parameter. The same force 
	// can be used for 1D, 2D or 3D case. 
	template<int dim>
	void 
	ForceType2<dim>
	::source_term(const std::vector<Point<dim>> &p,
				 std::vector<Vector<double>> &value,
				 const double factor)
	{
		Assert(p.size()!=0,ExcNotInitialized());
		Assert(value.size()!=0,ExcNotInitialized());
		AssertDimension(p.size(),value.size());
		

		for (unsigned int i = 0 ; i < value.size() ; i++)
			{
				const double norm = sqrt(p[i].square());
				const double x_coord = p[i][0];

				value[i] = 0;
				AssertDimension((int)value[i].size(),this->nEqn);
				value[i](var_rho) = factor * (this->constants.A0 + this->constants.A2 * norm * norm + 
							this->constants.A1*(1.0-5.0/18*norm*norm/(this->constants.tau*this->constants.tau))*x_coord / norm);
			}

	}

	//forcing for poisson heat conduction
	template<int dim>
	class 
	ForceType3:public Base_ForceType<dim>
	{
		public:
			ForceType3(const constant_numerics &constants,const int nEqn);

			unsigned int var_theta;

		// The source term
		virtual void source_term(const std::vector<Point<dim>> &p,
								 std::vector<Vector<double>> &value,
								 const double factor);
	};

	template<int dim>
	ForceType3<dim>
	::ForceType3(const constant_numerics &constants,const int nEqn)
	:
	Base_ForceType<dim>(constants,nEqn)
	{;}

	// specialization for 1D
	template<>
	void 
	ForceType3<1>
	::source_term(const std::vector<Point<1>> &p,
				 std::vector<Vector<double>> &value,
				 const double factor)
	{
		Assert(p.size()!=0,ExcNotInitialized());
		Assert(value.size()!=0,ExcNotInitialized());
		AssertDimension(p.size(),value.size());
		var_theta = 2;
		

		// now we need to check whether we captured the correct variable or not for the poisson heat conduction problem
		Assert(var_theta == 2,ExcMessage("Forcing not being applied to the energy equation"));

		for (unsigned int i = 0 ; i < value.size() ; i++)
			{
				const double x_cord = p[i][0];
				// initialize the variable
				value[i] = 0;
				AssertDimension((int)value[i].size(),this->nEqn);

				// the source term for the energy equation. The factor appears from the symmetrizer
				value[i](var_theta) = -factor * this->constants.alpha * pow(x_cord,2);
			}

	}

	
	// specialization for 2D
	template<>
	void 
	ForceType3<2>
	::source_term(const std::vector<Point<2>> &p,
				 std::vector<Vector<double>> &value,
				 const double factor)
	{
		Assert(p.size()!=0,ExcNotInitialized());
		Assert(value.size()!=0,ExcNotInitialized());
		AssertDimension(p.size(),value.size());
		var_theta = 3;
		

		// now we need to check whether we captured the correct variable or not for the poisson heat conduction problem
		Assert(var_theta == 3,ExcMessage("Forcing not being applied to the energy equation; check var_theta."));

		for (unsigned int i = 0 ; i < value.size() ; i++)
			{
				const double y_cord = p[i][1];
				// initialize the variable
				value[i] = 0;
				AssertDimension((int)value[i].size(),this->nEqn);

				// the source term for the energy equation.
				// The factor appears from the symmetrizer
				value[i](var_theta) = -factor * this->constants.alpha * pow(y_cord,2);
			}

	}



	// same as above but for the 3D case
	template<>
	void 
	ForceType3<3>
	::source_term(const std::vector<Point<3>> &p,
				 std::vector<Vector<double>> &value,
				 const double factor)
	{
		Assert(p.size()!=0,ExcNotInitialized());
		Assert(value.size()!=0,ExcNotInitialized());
		AssertDimension(p.size(),value.size());
		var_theta = 4;
		

		// now we need to check whether we captured the correct variable or not for the poisson heat conduction problem
		Assert(var_theta == 4,ExcMessage("Forcing not being applied to the energy equation"));

		for (unsigned int i = 0 ; i < value.size() ; i++)
			{
				const double x_cord = p[i][0];
				// initialize the variable
				value[i] = 0;
				AssertDimension((int)value[i].size(),this->nEqn);

				// the source term for the energy equation. The factor appears from the symmetrizer
				value[i](var_theta) = -factor * this->constants.alpha * pow(x_cord,2);
			}

	}
}