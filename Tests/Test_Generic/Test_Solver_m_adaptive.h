using namespace dealii;

TEST(MAdaptiveSolver,HandlesMAdaptiveSolver)
{
		const unsigned int dim = 1;

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
		

		// the tests have not been implemented for more than one system
		AssertDimension(constants.constants_sys.total_systems,1);
		Assert(constants.constants_num.problem_type != periodic,ExcNotImplemented());

		if(constants.constants_num.problem_type != periodic)
		{
			if (dim == 2)
			{

				// the exact solution can only be created for one of the systems
				ExactSolution::ExactSolution_Dummy<dim>  exactsolution_dummy(constants.constants_num,System[constants.constants_sys.total_systems-1].base_tensorinfo.S_half,
																			 constants.constants_sys.nEqn[constants.constants_sys.total_systems-1],constants.constants_sys.Ntensors[constants.constants_sys.total_systems-1]);

				FEM_Solver::Base_Solver<dim> base_solver("grid",
											 	 constants.constants_num,
											 	 System,
											 	 &exactsolution_dummy,
											 	 constants.constants_sys.nEqn,
											 	 constants.constants_sys.nBC);

			
				base_solver.run();


				Assert(constants.constants_num.problem_type == heat_conduction || constants.constants_num.problem_type == inflow_outflow ||
					  constants.constants_num.problem_type == lid_driven_cavity,ExcNotImplemented());
				Assert(constants.constants_num.mesh_type == square_domain || constants.constants_num.mesh_type == NACA5012 || constants.constants_num.mesh_type == square_circular_cavity,ExcNotImplemented());
				AssertDimension(constants.constants_num.refine_cycles,1);
				AssertDimension(constants.constants_num.initial_refinement,1);

				// value only stored for these number of cells
				if (constants.constants_num.mesh_type == NACA5012)
					AssertDimension(base_solver.triangulation.n_active_cells(),53);

				if (constants.constants_num.mesh_type == square_domain && constants.constants_num.problem_type == inflow_outflow)				
					AssertDimension(base_solver.triangulation.n_active_cells(),61);


				if (constants.constants_num.mesh_type == square_circular_cavity && constants.constants_num.problem_type == inflow_outflow)
					AssertDimension(base_solver.triangulation.n_active_cells(),476);


				if (constants.constants_num.mesh_type == square_domain && constants.constants_num.problem_type != inflow_outflow)
				{
					AssertDimension(constants.constants_num.part_x,10);
					AssertDimension(constants.constants_num.part_y,10);
				}

				if (dim == 2)
					AssertDimension(constants.constants_sys.Ntensors[0],6);

				// for the heat conduction problem, this is the l2 norm of the temperature
				double error_manuel;

				if (constants.constants_num.mesh_type == square_domain && constants.constants_num.problem_type == heat_conduction)
					error_manuel =  1.8584849033516051;

				if (constants.constants_num.mesh_type == square_domain && constants.constants_num.problem_type == lid_driven_cavity)
					error_manuel =  1.2248747382201992;

				if (constants.constants_num.mesh_type == square_domain && constants.constants_num.problem_type == inflow_outflow)
					error_manuel =  1.8776535235692213;

				if (constants.constants_num.mesh_type == square_circular_cavity && constants.constants_num.problem_type == inflow_outflow)
					error_manuel =  1.7575244416510396;


				if (constants.constants_num.mesh_type == NACA5012)
					error_manuel =  1.1675197098722727;


				EXPECT_NEAR(base_solver.error_per_itr[0],error_manuel,1e-10);		

			}

			// incase of the 1D problem we have the exact solution since we solve the poisson heat conduction problem
			if (dim == 1)
			{

				ExactSolution::PoissonHeat<dim>  PoissonHeat(constants.constants_num,System[constants.constants_sys.total_systems-1].base_tensorinfo.S_half,
													constants.constants_sys.nEqn[constants.constants_sys.total_systems-1],constants.constants_sys.Ntensors[constants.constants_sys.total_systems-1]);

				FEM_Solver::Base_Solver<dim> base_solver("grid",
											 	 constants.constants_num,
											 	 System,
											 	 &PoissonHeat,
											 	 constants.constants_sys.nEqn,
											 	 constants.constants_sys.nBC);

			
				// the number of cells which should be initially present
				base_solver.run();

				fflush(stdout);

				AssertDimension(constants.constants_num.refine_cycles,3);
				AssertDimension(constants.constants_num.initial_refinement,1);
				std::vector<double> error_manuel(constants.constants_num.refine_cycles);

				if(constants.constants_sys.Ntensors[0] == 6)
				{
					error_manuel[0] = 4.5432e-03;
					error_manuel[1]	= 1.1759e-03; 
					error_manuel[2] = 2.4787e-04;
				}

				if (constants.constants_sys.Ntensors[0] == 8)
				{
					error_manuel[0] = 7.2349e-03;
					error_manuel[1]	= 2.1590e-03; 
					error_manuel[2] = 5.0865e-04;		
				}

				if (constants.constants_sys.Ntensors[0] == 9)
				{
					error_manuel[0] = 2.8423e-03;
					error_manuel[1]	= 6.8201e-04 ; 
					error_manuel[2] = 1.4148e-04 ;		
				}


				if (constants.constants_sys.Ntensors[0] == 11)
				{
					error_manuel[0] = 4.0009e-03;
					error_manuel[1]	= 1.0194e-03; 
					error_manuel[2] = 2.1659e-04;		
				}

				for (int i = 0 ; i < constants.constants_num.refine_cycles ; i++)
				{
					std::cout << "Error difference " << fabs(base_solver.error_per_itr[i] - error_manuel[i]) << std::endl;
					EXPECT_NEAR(base_solver.error_per_itr[i],error_manuel[i],1e-5);		
				}	
			}



		}


}