namespace Develop_System
{	
	using namespace dealii;

		// In the following class we will construct all the matrices associated with the systemA
	template<int dim>
	class
	Base_DevelopSystem:public EquationGenerator::Base_EquationGenerator<dim>
	{
		public:
			Base_DevelopSystem(const constant_numerics &constants,
								const int nEqn,
								const int nBC,
								const int Ntensors,
								std::string &folder_name);
		
	};

	// folder name is the name of the folder in which the system_matrices have been kept
	template<int dim>
	Base_DevelopSystem<dim>
	::Base_DevelopSystem(const constant_numerics &constants,
								const int nEqn,
								const int nBC,
								const int Ntensors,
								std::string &folder_name)
	:
	// we first initialize the base class 
	EquationGenerator::Base_EquationGenerator<dim>(constants,nEqn,nBC,Ntensors)
	{
		// initialize the system
		this->reinit_system(folder_name);
	}



	// and the other variation is the boundary conditions
	template<int dim>
	class
	System:public Base_DevelopSystem<dim>
	{
		public:
			// we need the input parameters to initialize the base class
			System(const constant_numerics &constants,
								const int nEqn,
								const int nBC,
								const int Ntensors,
								std::string &folder_name);


			// The following boundary routines are system dependent, that is 
			// they depend upon the boundary matrix B. Therefore they cannot be
			// a member of the base class.
			BCrhs::BCrhs_wall<dim> bcrhs_wall;
			BCrhs::BCrhs_inflow<dim> bcrhs_inflow;
			BCrhs::BCrhs_wall_kinetic<dim> bcrhs_wall_kinetic;
			BCrhs::BCrhs_inflow_kinetic<dim> bcrhs_inflow_kinetic;
			BCrhs_systemA::BCrhs_ring_char_systemA<dim> bcrhs_ring_char_systemA;


			void reinit_BCrhs();
			void initialize_system();

	};

	template<int dim>
	System<dim>
	::System(const constant_numerics &constants,
								const int nEqn,
								const int nBC,
								const int Ntensors,
								std::string &folder_name)
	:
	Base_DevelopSystem<dim>(constants,nEqn,nBC,Ntensors,folder_name),
	// CAUTION sending B without relaxational normal velocity
	bcrhs_wall(constants,nBC,this->system_data.B.matrix),
	bcrhs_inflow(constants,nBC,this->system_data.Binflow.matrix),
	bcrhs_wall_kinetic(constants,nEqn,this->system_data.rhoW.matrix),
	bcrhs_inflow_kinetic(constants,nEqn,this->system_data.rhoInflow.matrix),
	bcrhs_ring_char_systemA(constants,nBC)
	{
		// if system A then no odd boundary implementation
		if(Ntensors == 3)
			Assert(constants.bc_type == characteristic,ExcNotImplemented());
	}


	template<int dim>
	void 
	System<dim>::initialize_system()
	{
		this->initialized_system = true;

		// else
		this->base_tensorinfo.reinit();

		//initialize the forcing term for this system
		this->reinit_force();

		// initialize the boundary matrices for this system
		this->reinit_BoundaryMatrices();

		// we will always send Ax independent of symmetric or unsymmetric system
		this->reinit_Aminus1D();

   		this->reinit_Xminus();

   		this->symmetrize_system();

   		this->force_factor = this->forcing_factor();

   		this->reinit_BCrhs();		
	}

	template<int dim>
	void 
	System<dim>::reinit_BCrhs()
	{
		if (this->Ntensors == 3 && dim == 2)
			this->base_bcrhs = &bcrhs_ring_char_systemA;

		else
			this->base_bcrhs = &bcrhs_wall;
	}



}