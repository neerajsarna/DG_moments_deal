using namespace dealii;

// the following routine checks the number of active cell and the corresponding mesh

void check_active_cells(const enum Mesh_type &mesh_type,
						const enum Problem_type &problem_type,
						const unsigned long int n_cells)
{
		if (mesh_type == NACA5012)
			AssertDimension(n_cells,53);

		if (mesh_type == square_domain && problem_type == inflow_outflow)				
			AssertDimension(n_cells,61);


		if (mesh_type == square_circular_cavity && problem_type == inflow_outflow)
			AssertDimension(n_cells,476);




}


// to perform our tests completely we would like to save the values of the error and then in the end check whether 
// or not we can reproduce them
// the 2D case
void set_error_values(const enum Mesh_type &mesh_type,
					 const enum Problem_type &problem_type,
					 double &error_value)
{

		if (mesh_type == square_domain && problem_type == heat_conduction)
			error_value =  1.8584849033516051;

		if (mesh_type == square_domain && problem_type == lid_driven_cavity)
			error_value =  1.2248747382201992;

		if (mesh_type == square_domain && problem_type == inflow_outflow)
			error_value =  1.8776535235692213;

		if (mesh_type == square_circular_cavity && problem_type == inflow_outflow)
			error_value =  1.7575244416510396;

		if (mesh_type == NACA5012)
			error_value =  1.1675197098722727;	
}

// error values for the 1d heat conduction case
void set_error_values(const enum Mesh_type &mesh_type,
					 const enum Problem_type &problem_type,
					 std::vector<double> &error_value)
{
	error_value[0] = 4.5432e-03;
	error_value[1]	= 1.1759e-03; 
	error_value[2] = 2.4787e-04;
}


TEST(DISABLED_SolverSingleSystem,HandlesSolverSingleSystem)
{
		const unsigned int dim = 1;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);


		AssertDimension(constants.constants_sys.total_systems,1);


	// 	// we construct a vector of all the system data being considered 
		std::vector<Develop_System::System<dim>> System;

		// initialize the vector containing all the systems
		// initialize the vector containing all the systems
		for(int i = 0 ; i < constants.constants_sys.total_systems ; i++)
		{
			System.push_back(Develop_System::System<dim>(constants.constants_num,constants.constants_sys.nEqn[i],
				constants.constants_sys.nBC[i],constants.constants_sys.Ntensors[i],folder_name));

		}
	
		for (int i = 0 ; i < constants.constants_sys.total_systems ; i++)
			System[i].initialize_system();

		AssertDimension(constants.constants_sys.Ntensors[0],6);

		// we check the test case for which the triangulation has been saved 
		if (dim == 2)
		{
			Assert(constants.constants_num.problem_type == heat_conduction || constants.constants_num.problem_type == inflow_outflow ||
			constants.constants_num.problem_type == lid_driven_cavity,ExcNotImplemented());
			Assert(constants.constants_num.mesh_type == square_domain || constants.constants_num.mesh_type == NACA5012 || 
				constants.constants_num.mesh_type == square_circular_cavity ,ExcNotImplemented());
			AssertDimension(constants.constants_num.initial_refinement,1);

			if (constants.constants_num.mesh_type == square_domain && constants.constants_num.problem_type != inflow_outflow)
			{
				AssertDimension(constants.constants_num.part_x,10);
				AssertDimension(constants.constants_num.part_y,10);
			}

			int system_to_solve = 0;

			// in case of hp fe method we initialize with the highest possible moment system
			if (constants.constants_num.assembly_type == manuel)
				system_to_solve = constants.constants_sys.total_systems-1;

					// the exact solution can only be created for one of the systems
			ExactSolution::ExactSolution_Dummy<dim>  exactsolution_dummy(constants.constants_num,
																		 System[system_to_solve].base_tensorinfo.S_half,
																		 constants.constants_sys.nEqn[system_to_solve],
																		 constants.constants_sys.Ntensors[system_to_solve]);

			// finite element solver for a single system
			FEM_Solver::Run_Problem_FE<dim> fe_solver("grid",
											 	 constants.constants_num,
											 	 System,
											 	 &exactsolution_dummy,
											 	 constants.constants_sys.nEqn,
											 	 constants.constants_sys.nBC,
											 	 system_to_solve);

		// finite element solver for hp data structures, used for m adaptivity
		FEM_Solver::Run_Problem_hp_FE<dim> hp_solver("grid",
											 	 constants.constants_num,
											 	 System,
											 	 &exactsolution_dummy,
											 	 constants.constants_sys.nEqn,
											 	 constants.constants_sys.nBC);

		check_active_cells(fe_solver.constants.mesh_type,fe_solver.constants.problem_type,fe_solver.fe_data_structure.triangulation.n_active_cells());


			double error_manuel;
			set_error_values(fe_solver.constants.mesh_type,fe_solver.constants.problem_type,error_manuel);


		if (constants.constants_num.assembly_type == meshworker)
		{
			fe_solver.run();
			EXPECT_NEAR(fe_solver.error_per_itr[0],error_manuel,1e-10);	
		}

		if (constants.constants_num.assembly_type == manuel)
		{
			hp_solver.run_distribution_deviation();
			EXPECT_NEAR(hp_solver.error_per_itr[0],error_manuel,1e-10);	

		}


		}

		if (dim == 1) 
		{
			Assert(constants.constants_num.mesh_type == line ,
				  ExcMessage("Violates the only possible geometry in the 1D case."));

			AssertDimension(constants.constants_num.refine_cycles,3);
			AssertDimension(constants.constants_num.initial_refinement,1);
			// only add boundary
			AssertDimension(constants.constants_num.bc_type,1);

			const int system_to_solve = 0;

			// we take the solution of the Boltzmann equation as the reference
			ExactSolution::PoissonHeat<dim>  PoissonHeat(constants.constants_num,System[system_to_solve].base_tensorinfo.S_half,
																			 constants.constants_sys.nEqn[system_to_solve],
																			 constants.constants_sys.Ntensors[system_to_solve]);

			// finite element solver for a single system
			FEM_Solver::Run_Problem_FE<dim> fe_solver("grid",
											 	 constants.constants_num,
											 	 System,
											 	 &PoissonHeat,
											 	 constants.constants_sys.nEqn,
											 	 constants.constants_sys.nBC,
											 	 system_to_solve);

		// finite element solver for hp data structures, used for m adaptivity
		FEM_Solver::Run_Problem_hp_FE<dim> hp_solver("grid",
											 	 constants.constants_num,
											 	 System,
											 	 &PoissonHeat,
											 	 constants.constants_sys.nEqn,
											 	 constants.constants_sys.nBC);


			std::vector<double> error_manuel(constants.constants_num.refine_cycles);
			set_error_values(fe_solver.constants.mesh_type,fe_solver.constants.problem_type,error_manuel);


		if (constants.constants_num.assembly_type == meshworker)
		{
			fe_solver.run();
			

			for (int i = 0 ; i < constants.constants_num.refine_cycles ; i++)
				{
					std::cout << "Error difference " << fabs(fe_solver.error_per_itr[i] - error_manuel[i]) << std::endl;
					EXPECT_NEAR(fe_solver.error_per_itr[i],error_manuel[i],1e-7);		
				}	
		}

		if (constants.constants_num.assembly_type == manuel)
		{
			hp_solver.run_distribution_deviation();
			
			for (int i = 0 ; i < constants.constants_num.refine_cycles ; i++)
				{
					std::cout << "Error difference " << fabs(hp_solver.error_per_itr[i] - error_manuel[i]) << std::endl;
					EXPECT_NEAR(hp_solver.error_per_itr[i],error_manuel[i],1e-7);		
				}	

		}

		}
}

TEST(DISABLED_RunSystemA,HandlesSystemA)
{
		const int dim = 2;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);


	// 	// we construct a vector of all the system data being considered 
		std::vector<Develop_System::System<dim>> System;

		// initialize the vector containing all the systems
		// initialize the vector containing all the systems
		for(int i = 0 ; i < constants.constants_sys.total_systems ; i++)
		{
			System.push_back(Develop_System::System<dim>(constants.constants_num,constants.constants_sys.nEqn[i],
				constants.constants_sys.nBC[i],constants.constants_sys.Ntensors[i],folder_name));

		}
	
		for (int i = 0 ; i < constants.constants_sys.total_systems ; i++)
			System[i].initialize_system();


		// the exact solution can only be created for one of the systems
		// ExactSolution::ExactSolution_Dummy<dim>  dummy(constants.constants_num,System[constants.constants_sys.total_systems-1].base_tensorinfo.S_half,
		// 											constants.constants_sys.nEqn[constants.constants_sys.total_systems-1],constants.constants_sys.Ntensors[constants.constants_sys.total_systems-1]);


		const int system_to_solve = 0;
		ExactSolution::ExactSolution_SystemA_ring<dim>  exact_solution_systemA(constants.constants_num,
																			System[system_to_solve].base_tensorinfo.S_half,
																			constants.constants_sys.nEqn[system_to_solve],
																			constants.constants_sys.Ntensors[system_to_solve]);

		
		// finite element solver for a single system
		FEM_Solver::Run_Problem_FE<dim> fe_solver("grid",
											constants.constants_num,
											System,
											&exact_solution_systemA,
											constants.constants_sys.nEqn,
											constants.constants_sys.nBC,
											system_to_solve);


		fe_solver.run();

		AssertThrow(constants.constants_num.refine_cycles == 2,ExcMessage("The refine cycles requested for have not been implemented"));
		AssertThrow(constants.constants_num.p == 1,ExcNotImplemented());
		AssertThrow(constants.constants_num.mapping_order == 2 , ExcNotImplemented());
		AssertThrow(fabs(constants.constants_num.initial_refinement - 1) < 1e-5,ExcNotImplemented());

		Vector<double> exact_error(constants.constants_num.refine_cycles);

		if (fabs(constants.constants_num.tau - 0.1) < 1e-5)
		{
			if(constants.constants_num.bc_type == characteristic)
			{
				
				exact_error(0) = 1.6451e-01;
				exact_error(1) = 2.8996e-02;
			}
		}


		for (int i = 0 ; i < constants.constants_num.refine_cycles ; i++)
				EXPECT_NEAR(fe_solver.error_per_itr[i],exact_error(i),1e-5);	
}


// run without any restrictions
TEST(RunSystem,HandlesRunSystem)
{
		const unsigned int dim = 2;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);


		std::cout << "Creating Systems....." << std::endl;
	// 	// we construct a vector of all the system data being considered 
		std::vector<Develop_System::System<dim>> System;

		// initialize the vector containing all the systems
		// initialize the vector containing all the systems
		for(int i = 0 ; i < constants.constants_sys.total_systems ; i++)
		{
			System.push_back(Develop_System::System<dim>(constants.constants_num,constants.constants_sys.nEqn[i],
				constants.constants_sys.nBC[i],constants.constants_sys.Ntensors[i],folder_name));

		}
	
		for (int i = 0 ; i < constants.constants_sys.total_systems ; i++)
			System[i].initialize_system();


		std::cout << "Done creating systems....." << std::endl;

		
		// We initialize the exact solution with the system we wish to solve.
		// Incase of hp this can be changed to the system with 	
		// we change the system which we wish to solve
		const int system_to_solve = constants.constants_sys.total_systems-1;

		Assert(system_to_solve < System.size(),
			ExcMessage("You have asked for a system which has not been loaded"));

		ExactSolution::ExactSolution_Dummy<dim>  dummy(constants.constants_num,
												System[system_to_solve].base_tensorinfo.S_half,
												constants.constants_sys.nEqn[system_to_solve],
												constants.constants_sys.Ntensors[system_to_solve]);


			// finite element solver for a single system
		FEM_Solver::Run_Problem_FE<dim> fe_solver("grid",
			constants.constants_num,
			System,
			&dummy,
			constants.constants_sys.nEqn,
			constants.constants_sys.nBC,
			system_to_solve);


		fe_solver.run();

		


}


