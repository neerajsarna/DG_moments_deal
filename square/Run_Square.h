template<int dim>
void
Base_Solver<dim>::run_square()
{
	const unsigned int refine_cycles = constants.refine_cycles;
	error_per_itr.resize(refine_cycles);

	TimerOutput timer (std::cout, TimerOutput::summary,
                   TimerOutput::wall_times);

	Assert(constants.mesh_type == Mesh_type::square, ExcMessage("Inappropriate mesh type"));

	for (unsigned int i = 0 ; i < refine_cycles ; i ++)
	{
	
		AssertDimension((int)error_per_itr.size(),constants.refine_cycles);

		Assert(constants.matrix_type == Trilinos_Mat, ExcMessage("Algorithm only for an eigen system matrix"));
		//now we distribute the dofs and allocate the memory for the global matrix. We have 
		// already created the mesh so we can directly distribute the degrees of freedom now.

		std::cout << "Distributing dof " << std::endl;
		timer.enter_subsection("Dof Distribution");
		distribute_dof_allocate_matrix();
		timer.leave_subsection();

		// the following routine assembles
		std::cout << "Assembling" << std::endl;
		timer.enter_subsection("Assembly");
		switch (constants.assembly_type)
		{
			case meshworker:
			{
				std::cout << "Using meshwoker " << std::endl;
				assemble_system_meshworker();
				break;		
			}
			case manuel:
			{
				switch(constants.bc_type)
				{
					case characteristic:
					{
						std::cout << "Using characteristic boundary " << std::endl;
						assemble_system_char();
						break;
					}

					case odd:
					{
						std::cout << "Using odd boundary "<< std::endl;
						assemble_system_odd();
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
				Assert(1==0,ExcMessage("Should not have reached here"));
				break;
			}
		}
		timer.leave_subsection();
		

		// we initialize the object which will solve our system
		// We do int the following way so as to keep the solver independent of all the other implementations.
		// This makes the code highly reusable. So one can directly copy the following class and use it somewhere
		// else is one wants to.

		timer.enter_subsection("Linear Solver");
		LinearSolver::LinearSolver linear_solver;
		linear_solver.develop_pardiso_data(global_matrix,sparsity_pattern);
		double residual = linear_solver.solve_with_pardiso(system_rhs,solution);
		timer.leave_subsection();

		// we do not have the exact solution so we dont do post processing.
		// timer.enter_subsection("Post Processing");
		PostProc::Base_PostProc<dim> postproc(constants,base_exactsolution,&dof_handler, &mapping);

		// now we compute the error due to computation
		const double l2_norm = postproc.compute_L2_norm(solution,this->triangulation.n_active_cells());
		std::cout << "l2_norm of error variable " << l2_norm << "Residual " << residual <<  std::endl;

		 postproc.print_options(this->triangulation,solution,i,refine_cycles,convergence_table);		
		// timer.leave_subsection();

		// Grid refinement should be done in the end.
		this->refinement_handling(i,refine_cycles);

	}

}