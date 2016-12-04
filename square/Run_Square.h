template<int dim>
void
Base_Solver<dim>::run_square()
{
	const unsigned int refine_cycles = constants.refine_cycles;
	error_per_itr.resize(refine_cycles);

	TimerOutput timer (std::cout, TimerOutput::summary,
                   TimerOutput::wall_times);

	for (unsigned int i = 0 ; i < refine_cycles ; i ++)
	{
	
		AssertDimension((int)error_per_itr.size(),constants.refine_cycles);

		//now we distribute the dofs and allocate the memory for the global matrix. We have 
		// already created the mesh so we can directly distribute the degrees of freedom now.
		timer.enter_subsection("distribute dof");
		distribute_dof_allocate_matrix();
		timer.leave_subsection();

		// the following routine assembles
		timer.enter_subsection("Assembly");
		assemble_system_meshworker();
		timer.leave_subsection();

		// we initialize the object which will solve our system
		// We do int the following way so as to keep the solver independent of all the other implementations.
		// This makes the code highly reusable. So one can directly copy the following class and use it somewhere
		// else is one wants to.
		timer.enter_subsection("Linear Solver");
		LinearSolver::LinearSolver linear_solver(global_matrix,system_rhs,solution);
		timer.leave_subsection();

		timer.enter_subsection("Post Processing");
		PostProc::Base_PostProc<dim> postproc(constants,base_exactsolution,&dof_handler, &mapping);

		// we do not have a reference yet so no error evaluation
		postproc.print_options(this->triangulation,solution,i,refine_cycles,convergence_table);		
		timer.leave_subsection();

		// Grid refinement should be done in the end.
		this->refinement_handling(i,refine_cycles);

	}

}