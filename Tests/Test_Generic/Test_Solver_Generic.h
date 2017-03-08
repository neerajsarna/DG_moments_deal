using namespace dealii;
	TEST(Solver,HandlesSolver)
	{
		const unsigned int dim = 2;
		

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		Develop_System::System<dim> System(constants.constants,folder_name);

		if(constants.constants.problem_type != periodic)
		{
			if (dim == 2)
			{

				ExactSolution::ExactSolution_Dummy<dim>  exactsolution_dummy(constants.constants,System.base_tensorinfo.S_half);

				FEM_Solver::Base_Solver<dim> base_solver("grid",
											 	 constants.constants,
											 	 &System,
											 	 &exactsolution_dummy);

			
				base_solver.run();
			}

			// incase of the 1D problem we have the exact solution since we solve the poisson heat conduction problem
			if (dim == 1)
			{

				ExactSolution::PoissonHeat<dim>  PoissonHeat(constants.constants,System.base_tensorinfo.S_half);

				FEM_Solver::Base_Solver<dim> base_solver("grid",
											 	 constants.constants,
											 	 &System,
											 	 &PoissonHeat);

			
				base_solver.run();
			}



		}

		else
		{
			Assert(dim == 2,ExcNotImplemented());

			ExactSolution::PoissonHeat<dim>  PoissonHeat(constants.constants,System.base_tensorinfo.S_half);



			FEM_Solver::Base_Solver<dim> base_solver("grid",
													 constants.constants,
													 &System,
													 &PoissonHeat);

			base_solver.print_mesh_info();
			fflush(stdout);

			base_solver.run_periodic();


			Assert(constants.constants.refine_cycles = 3,ExcNotImplemented());
			Assert(fabs(constants.constants.tau - 0.1) < 1e-5,ExcNotImplemented());
			Assert(constants.constants.part_y == 100,ExcNotImplemented());

			Vector<double> exact_error(3);

			if (constants.constants.Ntensors == 6)
			{
				if (constants.constants.bc_type == odd)
				{
					exact_error(0) = 2.108024e-06;
					exact_error(1) = 5.096096e-07;
					exact_error(2) =  2.236106e-07;				
				}
				else
				{
					exact_error(0) = 1.923850e-06;
					exact_error(1) = 4.848456e-07;
					exact_error(2) = 2.160810e-07;
				}



				for (int i = 0 ; i < constants.constants.refine_cycles ; i++)
				{
					std::cout << "Error difference " << fabs(base_solver.error_per_itr[i] - exact_error(i)) << std::endl;
					EXPECT_NEAR(base_solver.error_per_itr[i],exact_error(i),1e-10);		
				}				
			}

		}
		

	}