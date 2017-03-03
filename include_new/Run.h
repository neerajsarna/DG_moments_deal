using namespace dealii;
void Run()
{
	const unsigned int dim = 2;
	Constants::Base_Constants constants(input_file);
	std::string folder_name = "../system_matrices/";

	switch(constants.constants.nEqn)
	{
		case 6:
		{
			// first develop the class which loads the systems
			SystemA::SystemA<dim> systemA(constants.constants,folder_name);			

			// routine to develop the exact solution
			ExactSolution::ExactSolution_SystemA_ring<dim>  exact_solution_systemA(constants.constants,systemA.base_tensorinfo.S_half);

			// the solver which combines the above two classes
			FEM_Solver::Base_Solver<dim> base_solver(
													 "grid",
													 constants.constants,
													 &systemA,
													 &exact_solution_systemA);

			AssertThrow(constants.constants.mesh_type == ring,ExcMessage("Incorrect mesh type"));
			// now run the simulations 
			base_solver.run();
			break;
		}

		case 13:
		{

			G20::G20<dim> G20(constants.constants,folder_name);

			if (G20.constants.problem_type == periodic)
			{

					ExactSolution::G20_PoissonHeat<dim>  G20_PoissonHeat(constants.constants,G20.base_tensorinfo.S_half);

					FEM_Solver::Base_Solver<dim> base_solver(
															 "grid",
															 constants.constants,
															 &G20,
															 &G20_PoissonHeat);

					AssertThrow(constants.constants.mesh_type == square_domain,ExcMessage("Incorrect mesh type"));
					base_solver.run_periodic();

			}
			else 
			{
					// then we create a dummy exact solution
					ExactSolution::ExactSolution_Dummy<dim>  exactsolution_dummy(constants.constants,G20.base_tensorinfo.S_half);

					// we intialize the solver class
					FEM_Solver::Base_Solver<dim> base_solver(
													 "grid",
													 constants.constants,
													 &G20,
													 &exactsolution_dummy);
					base_solver.run();
					
			}
			break;

		}


		case 17:
		{

			G26::G26<dim> G26(constants.constants,folder_name);

			if (G26.constants.problem_type == periodic)
			{
					ExactSolution::G26_PoissonHeat<dim>  G26_PoissonHeat(constants.constants,G26.base_tensorinfo.S_half);

					FEM_Solver::Base_Solver<dim> base_solver(
						"grid",
						constants.constants,
						&G26,
						&G26_PoissonHeat);

					AssertThrow(constants.constants.mesh_type == square_domain,ExcMessage("Incorrect mesh type"));

					base_solver.run_periodic();
			}

			else
			{
					ExactSolution::ExactSolution_Dummy<dim>  exactsolution_dummy(constants.constants,G26.base_tensorinfo.S_half);

					FEM_Solver::Base_Solver<dim> base_solver(
													 "grid",
													 constants.constants,
													 &G26,
													 &exactsolution_dummy);


					base_solver.run();

			}
			break;
		}

		case 22:
		{

			G35::G35<dim> G35(constants.constants,folder_name);

			if(G35.constants.problem_type == periodic)
			{
					ExactSolution::G35_PoissonHeat<dim>  G35_PoissonHeat(constants.constants,G35.base_tensorinfo.S_half);

					FEM_Solver::Base_Solver<dim> base_solver(
													 "grid",
													 constants.constants,
													 &G35,
													 &G35_PoissonHeat);

					AssertThrow(constants.constants.mesh_type == square_domain,ExcMessage("Incorrect mesh type"));

					base_solver.run_periodic();
			}

			else
			{
					ExactSolution::ExactSolution_Dummy<dim>  exactsolution_dummy(constants.constants,G35.base_tensorinfo.S_half);

					FEM_Solver::Base_Solver<dim> base_solver(
													 		 "grid",
													 		  constants.constants,
													 		 &G35,
													 		 &exactsolution_dummy);
					base_solver.run();
			}

			break;

		}

		case 28:
		{
			G45::G45<dim> G45(constants.constants,folder_name);

			if (G45.constants.problem_type == periodic)
			{
					ExactSolution::G45_PoissonHeat<dim>  G45_PoissonHeat(constants.constants,
														G45.base_tensorinfo.S_half);

					FEM_Solver::Base_Solver<dim> base_solver(
													 "grid",
													 constants.constants,
													 &G45,
													 &G45_PoissonHeat);

					AssertThrow(constants.constants.mesh_type == square_domain,ExcMessage("Incorrect mesh type"));

					base_solver.run_periodic();
			}
			else
			{
					ExactSolution::ExactSolution_Dummy<dim>  exactsolution_dummy(constants.constants,
																G45.base_tensorinfo.S_half);

					FEM_Solver::Base_Solver<dim> base_solver(
													 "grid",
													 constants.constants,
													 &G45,
													 &exactsolution_dummy);
					base_solver.run();
			}

			break;
		}

		case 34:
		{

			G56::G56<dim> G56(constants.constants,folder_name);

			if (G56.constants.problem_type == periodic)
			{
					ExactSolution::G56_PoissonHeat<dim>  G56_PoissonHeat(constants.constants,
															G56.base_tensorinfo.S_half);

					FEM_Solver::Base_Solver<dim> base_solver(
													 "grid",
													 constants.constants,
													 &G56,
													 &G56_PoissonHeat);

					AssertThrow(constants.constants.mesh_type == square_domain,ExcMessage("Incorrect mesh type"));

					base_solver.run_periodic();
			}

			else
			{
					ExactSolution::ExactSolution_Dummy<dim>  exactsolution_dummy(constants.constants,G56.base_tensorinfo.S_half);

					FEM_Solver::Base_Solver<dim> base_solver(
													 "grid",
													 constants.constants,
													 &G56,
													 &exactsolution_dummy);
					base_solver.run();
			}

			break;
		}

		case 43:
		{

			G71::G71<dim> G71(constants.constants,folder_name);

			if (G71.constants.problem_type == periodic)
			{
					ExactSolution::G71_PoissonHeat<dim>  G71_PoissonHeat(constants.constants,
															G71.base_tensorinfo.S_half);

					FEM_Solver::Base_Solver<dim> base_solver(
													 "grid",
													 constants.constants,
													 &G71,
													 &G71_PoissonHeat);

					AssertThrow(constants.constants.mesh_type == square_domain,ExcMessage("Incorrect mesh type"));

					base_solver.run_periodic();				
			}
			else
			{
					ExactSolution::ExactSolution_Dummy<dim>  exactsolution_dummy(constants.constants,
																	G71.base_tensorinfo.S_half);

					FEM_Solver::Base_Solver<dim> base_solver(
													 "grid",
													 constants.constants,
													 &G71,
													 &exactsolution_dummy);
					base_solver.run();
			}


			break;
		}

		case 50:
		{

			G84::G84<dim> G84(constants.constants,folder_name);

			if (G84.constants.problem_type == periodic)
			{
					ExactSolution::G84_PoissonHeat<dim>  G84_PoissonHeat(constants.constants,
															G84.base_tensorinfo.S_half);

					FEM_Solver::Base_Solver<dim> base_solver(
													 "grid",
													 constants.constants,
													 &G84,
													 &G84_PoissonHeat);

					AssertThrow(constants.constants.mesh_type == square_domain,ExcMessage("Incorrect mesh type"));


					base_solver.run_periodic();				
			}
			else
			{
					ExactSolution::ExactSolution_Dummy<dim>  exactsolution_dummy(constants.constants,G84.base_tensorinfo.S_half);

					FEM_Solver::Base_Solver<dim> base_solver(
													 "grid",
													 constants.constants,
													 &G84,
													 &exactsolution_dummy);
					base_solver.run();
			}

		
			break;
		}

		case 62:
		{


			G105::G105<dim> G105(constants.constants,folder_name);

			if (G105.constants.problem_type == periodic)
			{
					ExactSolution::G105_PoissonHeat<dim>  G105_PoissonHeat(constants.constants,G105.base_tensorinfo.S_half);

					FEM_Solver::Base_Solver<dim> base_solver(
													 "grid",
													 constants.constants,
													 &G105,
													 &G105_PoissonHeat);

					AssertThrow(constants.constants.mesh_type == square_domain,ExcMessage("Incorrect mesh type"));

					base_solver.run_periodic();				
			}
			else
			{
					ExactSolution::ExactSolution_Dummy<dim>  exactsolution_dummy(constants.constants,G105.base_tensorinfo.S_half);

					FEM_Solver::Base_Solver<dim> base_solver(
													 "grid",
													 constants.constants,
													 &G105,
													 &exactsolution_dummy);

					base_solver.run();
			}

			break;
		}

		case 70:
		{

			G120::G120<dim> G120(constants.constants,folder_name);

			if(G120.constants.problem_type == periodic)
			{
					ExactSolution::G120_PoissonHeat<dim>  G120_PoissonHeat(constants.constants,G120.base_tensorinfo.S_half);

					FEM_Solver::Base_Solver<dim> base_solver(
													 "grid",
													 constants.constants,
													 &G120,
													 &G120_PoissonHeat);

					AssertThrow(constants.constants.mesh_type == square_domain,ExcMessage("Incorrect mesh type"));

					base_solver.run_periodic();				
			}

			else
			{
					ExactSolution::ExactSolution_Dummy<dim>  exactsolution_dummy(constants.constants,G120.base_tensorinfo.S_half);

					FEM_Solver::Base_Solver<dim> base_solver(
															 "grid",
															 constants.constants,
															 &G120,
															 &exactsolution_dummy);
					base_solver.run();				
			}

			break;
		}


		default:
		{
			AssertThrow(1 == 0,ExcMessage("Should not have reached here"));
			break;
		}
	}


}
