template<int dim>
void
Base_Solver<dim>::run_ring(const unsigned int refine_cycles)
{
	for (unsigned int i = 0 ; i < refine_cycles ; i ++)
	{
	

		//now we distribute the dofs and allocate the memory for the global matrix. We have 
		// already created the mesh so we can directly distribute the degrees of freedom now.
		distribute_dof_allocate_matrix();

		// the following routine assembles
		switch(constants.assembly_type)
		{
			case meshworker:
			{
				assemble_system_meshworker();
				break;
			}

			// case manuel:
			// {
			// 	assemble_system();
			// 	break;
			// }
			default :
			{
				Assert(1 == 0,ExcMessage("Should not have reached here"));
				break;
			}
		}

		std::cout << "Memory Before Compression " << global_matrix.memory_consumption() << std::endl;

		if(!global_matrix.is_compressed())
			global_matrix.compress();

		std::cout << "Memory After Compression " << global_matrix.memory_consumption() << std::endl;

		// we initialize the object which will solve our system
		// We do int the following way so as to keep the solver independent of all the other implementations.
		// This makes the code highly reusable. So one can directly copy the following class and use it somewhere
		// else is one wants to.
		LinearSolver::LinearSolver linear_solver(global_matrix,system_rhs,solution);

		 PostProc::Base_PostProc<dim> postproc(constants,base_exactsolution,&dof_handler, &mapping);

		// // now we compute the error due to computation
		postproc.error_evaluation_QGauss(solution,this->triangulation.n_active_cells(),
										GridTools::maximal_cell_diameter(this->triangulation),convergence_table);

		postproc.print_options(this->triangulation,solution,i,refine_cycles,convergence_table);		


		// Grid refinement should be done in the end.
		this->refinement_handling(i,refine_cycles);

	}

}