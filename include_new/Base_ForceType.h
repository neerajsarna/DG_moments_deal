// in the following namespace we declare the different forces we can have in the system
namespace ForceType
{
	using namespace dealii;

	template<int dim>
	class
	Base_ForceType
	{
		public:
			Base_ForceType(const constant_data &constants);

			const constant_data constants;
		public:
		// The source term
		// The term factor appears due to symmetrization
			virtual void source_term(const std::vector<Point<dim>> &p,
									 std::vector<Vector<double>> &value,
									 const double factor) = 0;
	};

	template<int dim>
	Base_ForceType<dim>::
	Base_ForceType(const constant_data &constants)
	:
	constants(constants)
	{;}

	template<int dim>
	class 
	ForceType1:public Base_ForceType<dim>
	{
		public:
			ForceType1(const constant_data &constants);

			// the force will be acting on the first equation
			const unsigned int var_rho = 0;

			// The source term
			virtual void source_term(const std::vector<Point<dim>> &p,
									 std::vector<Vector<double>> &value,
									 const double factor);
	};

	template<int dim>
	ForceType1<dim>
	::ForceType1(const constant_data &constants)
	:
	Base_ForceType<dim>(constants)
	{;}

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
				double norm = sqrt(p[i].square());
				value[i] = 0;
				AssertDimension((int)value[i].size(),this->constants.nEqn);
				value[i](var_rho) = factor * (this->constants.A0 + this->constants.A2*norm*norm + this->constants.A1*p[i][0]/norm);	
			}

	}


	// Force corresponding to exact solution of systemB
	template<int dim>
	class 
	ForceType2:public Base_ForceType<dim>
	{
		public:
			ForceType2(const constant_data &constants);

			const unsigned int var_rho = 0;
		// The source term
		virtual void source_term(const std::vector<Point<dim>> &p,
								 std::vector<Vector<double>> &value,
								 const double factor);
	};

	template<int dim>
	ForceType2<dim>
	::ForceType2(const constant_data &constants)
	:
	Base_ForceType<dim>(constants)
	{;}

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
				double norm = sqrt(p[i].square());
				value[i] = 0;
				AssertDimension((int)value[i].size(),this->constants.nEqn);
				value[i](var_rho) = factor * (this->constants.A0 + this->constants.A2 * norm * norm + 
							this->constants.A1*(1.0-5.0/18*norm*norm/(this->constants.tau*this->constants.tau))*p[i][0] / norm);
			}

	}

	// force corresponding to couette flow
	template<int dim>
	class 
	ForceType3:public Base_ForceType<dim>
	{
		public:
			ForceType3(const constant_data &constants);

			unsigned int var_theta;

		// The source term
		virtual void source_term(const std::vector<Point<dim>> &p,
								 std::vector<Vector<double>> &value,
								 const double factor);
	};

	template<int dim>
	ForceType3<dim>
	::ForceType3(const constant_data &constants)
	:
	Base_ForceType<dim>(constants)
	{;}

	template<int dim>
	void 
	ForceType3<dim>
	::source_term(const std::vector<Point<dim>> &p,
				 std::vector<Vector<double>> &value,
				 const double factor)
	{
		Assert(p.size()!=0,ExcNotInitialized());
		Assert(value.size()!=0,ExcNotInitialized());
		AssertDimension(p.size(),value.size());
		var_theta = this->constants.variable_map.find(this->constants.force_variable)->second;

		// now we need to check whether we captured the correct variable or not for the poisson heat conduction problem
		Assert(var_theta == dim + 1,ExcMessage("Forcing not being applied to the energy equation"));

		for (unsigned int i = 0 ; i < value.size() ; i++)
			{
				const double y_cord = p[i][1];
				// initialize the variable
				value[i] = 0;
				AssertDimension((int)value[i].size(),this->constants.nEqn);

				// the source term for the energy equation. The factor appears from the symmetrizer
				value[i](var_theta) = -factor * this->constants.alpha * pow(y_cord,2);
			}

	}
}