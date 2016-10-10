// In the following file we declare all the attributes of the A system
namespace SystemA
{
	using namespace dealii;


	// In the following class we will construct all the matrices associated with the systemA
	template<int dim>
	class
	Base_SystemA:public EquationGenerator::Base_EquationGenerator<dim>
	{
		public:
			Base_SystemA(const constant_data &constants,
						std::string &folder_name);
		

			MatrixUI varIdx;

			// special function to create varIdx for this particular system
			MatrixUI build_varIdx();

			// prepare varIdx for systemA, since we do not have a one to one correspondence with 
			const unsigned int Ntensors = 3;
	};

	// folder name is the name of the folder in which the system_matrices have been kept
	template<int dim>
	Base_SystemA<dim>
	::Base_SystemA(const constant_data &constants,
				   std::string &folder_name)
	:
	// we need to first initialize the constant structure
	EquationGenerator::Base_EquationGenerator<dim>(constants)
	{

		// initialize the system
		this->reinit_system(folder_name);

		varIdx = build_varIdx();
	}


	template<int dim>
	MatrixUI
	Base_SystemA<dim>
	::build_varIdx()
	{
		MatrixUI result;
		result.resize(Ntensors,2);

		Assert(result.rows() !=0 ,ExcNotInitialized());
		Assert(result.cols() !=0 ,ExcNotInitialized());
		
		// theta
		result(0,0) = 1;
		result(0,1) = 0;

		// q
		result(1,0) = 0;
		result(1,1) = 1;

		// R
		result(2,0) = 1;
		result(2,1) = 2;

		return(result);
	}

	// accomodation of different test case.
	// the test cases vary in the sense of the force type
	// and the other variation is the boundary conditions
	template<int dim>
	class
	SystemA:public Base_SystemA<dim>
	{
		public:
			// we need the input parameters to initialize the base class
			SystemA(const constant_data &constants,
					std::string &folder_name);

			virtual void reinit_BCrhs();



			// boundary routines for system A
			BCrhs_systemA::BCrhs_ring_char_systemA<dim> bcrhs_ring_char_systemA;
			BCrhs_systemA::BCrhs_ring_odd_systemA<dim> bcrhs_ring_odd_systemA;
			
			BCrhs_systemA::BCrhs_periodic_char_systemA<dim> bcrhs_periodic_char_systemA;
			BCrhs_systemA::BCrhs_periodic_odd_systemA<dim> bcrhs_periodic_odd_systemA;

	};

	template<int dim>
	SystemA<dim>
	::SystemA(const constant_data &constants,
			  std::string &folder_name)
	:
	Base_SystemA<dim>(constants,folder_name)
	bcrhs_ring_char_systemA(constants),
	bcrhs_ring_odd_systemA(constants),
	bcrhs_periodic_char_systemA(constants),
	bcrhs_periodic_odd_systemA(constants)
	{

		// we reinitialize all the data for base_tensorinfo for this particular system
		this->base_tensorinfo.reinit(this->varIdx,this->Ntensors);

		//initialize the forcing term for this system
		this->reinit_force();

		// initialize the boundary matrices for this system
		this->reinit_BoundaryMatrices();

		reinit_BCrhs();

		// we will always send Ax independent of symmetric or unsymmetric system
		this->reinit_Aminus1D();

   		this->reinit_Xminus();

   		this->symmetrize_system();

   		this->force_factor = this->forcing_factor();
	}
	

	template<int dim>
	void
	SystemA<dim>
	::reinit_BCrhs()
	{
					//the following implementation is mesh dependent
		switch(this->constants.mesh_type)
		{
			case ring:
			{
				switch (this->constants.bc_type)
				{
					case characteristic:
					{
						this->base_bcrhs = &bcrhs_ring_char_systemA;
						break;
					}
					case odd:
					{
						this->base_bcrhs = &this->bcrhs_ring_odd_systemA;
						break;
					}
				}
				break;
			}

			case periodic_square:
			{
				switch(constants.bc_type)
				{
					case characteristic:
					{
						this->base_bcrhs = &bcrhs_periodic_char_systemA;
						break;
					}

					case odd:
					{
						this->base_bcrhs = &bcrhs_periodic_odd_systemA;
						break;
					}

					default:
					{
						Assert(1 == 0, ExcMessage("Should not have reached here"));
						break;
					}
				}

				break;
			}
			default:
			{
				Assert(1 ==0,ExcNotImplemented());
				break;
			}
		}
		
	}
}


