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

			void reinit_BCrhs();


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
	bcrhs_inflow(constants,nBC,this->system_data.Binflow.matrix)
	{

		// we reinitialize all the data for base_tensorinfo for this particular system
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

	}

	template<int dim>
	void 
	System<dim>::reinit_BCrhs()
	{
		this->base_bcrhs = &bcrhs_wall;
	}


}