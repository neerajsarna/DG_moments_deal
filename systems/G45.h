namespace G45
{
	using namespace dealii;

		// In the following class we will construct all the matrices associated with the systemA
	template<int dim>
	class
	Base_G45:public EquationGenerator::Base_EquationGenerator<dim>
	{
		public:
			Base_G45(const constant_data &constants,
						std::string &folder_name);
			

			// prepare varIdx for systemA, since we do not have a one to one correspondence with 
			const unsigned int Ntensors = 11;
	};

	// folder name is the name of the folder in which the system_matrices have been kept
	template<int dim>
	Base_G45<dim>
	::Base_G45(const constant_data &constants,
				   std::string &folder_name)
	:
	// we need to first initialize the constant structure
	EquationGenerator::Base_EquationGenerator<dim>(constants)
	{

		// initialize the system
		this->reinit_system(folder_name);
	}



	// and the other variation is the boundary conditions
	template<int dim>
	class
	G45:public Base_G45<dim>
	{
		public:
			// we need the input parameters to initialize the base class
			G45(const constant_data &constants,
					std::string &folder_name);


			// The following boundary routines are system dependent, that is 
			// they depend upon the boundary matrix B. Therefore they cannot be
			// a member of the base class.
			BCrhs::BCrhs_wall<dim> bcrhs_wall;
			BCrhs::BCrhs_inflow<dim> bcrhs_inflow;


			virtual void reinit_BCrhs();

	};

	template<int dim>
	G45<dim>
	::G45(const constant_data &constants,
			  std::string &folder_name)
	:
	Base_G45<dim>(constants,folder_name),
	// CAUTION sending B without relaxational normal velocity
	bcrhs_wall(constants,this->system_data.B.matrix),
	bcrhs_inflow(constants,this->system_data.Binflow.matrix)
	{

		// we reinitialize all the data for base_tensorinfo for this particular system
		this->base_tensorinfo.reinit(this->Ntensors);

		//initialize the forcing term for this system
		this->reinit_force();

		this->reinit_BCrhs();

		// initialize the boundary matrices for this system
		this->reinit_BoundaryMatrices();

		// we will always send Ax independent of symmetric or unsymmetric system
		this->reinit_Aminus1D();

   		this->reinit_Xminus();

   		this->symmetrize_system();

   		this->force_factor = this->forcing_factor();
	}

	template<int dim>
	void
	G45<dim>
	::reinit_BCrhs()
	{


		this->base_bcrhs = &bcrhs_wall;
		this->base_bcrhs_inflow = &bcrhs_inflow;
	}
}
