using namespace dealii;
void Run()
{
	const unsigned int dim = 2;
	Constants::Base_Constants constants(input_file);
	std::string folder_name = "../system_matrices/";

	// first we handle the 2D case
	if (dim == 2)
	{
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

			default:
			{
				Develop_System::System<dim> System(constants.constants,folder_name);



				if (constants.constants.problem_type == periodic)
				{
				// prepares the solution for the poisson heat conduction problem
					ExactSolution::PoissonHeat<dim>  PoissonHeat(constants.constants,System.base_tensorinfo.S_half);

					FEM_Solver::Base_Solver<dim> base_solver(
						"grid",
						constants.constants,
						&System,
						&PoissonHeat);

					AssertThrow(constants.constants.mesh_type == square_domain,ExcMessage("Incorrect mesh type"));
					base_solver.run_periodic();
				}
				if (constants.constants.problem_type != periodic)
				{
				// prepare the dummy solutoin which returns zero value
					ExactSolution::ExactSolution_Dummy<dim>  exactsolution_dummy(constants.constants,System.base_tensorinfo.S_half);

				// we intialize the solver class
					FEM_Solver::Base_Solver<dim> base_solver(
						"grid",
						constants.constants,
						&System,
						&exactsolution_dummy);
					
					base_solver.run();
				}

				break
			}
		}	
	}
	

}
