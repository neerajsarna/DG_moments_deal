template<int dim>
void
Base_Solver<dim>::distribute_dof_allocate_matrix_ring()
{
    dof_handler.distribute_dofs(finite_element);

    DynamicSparsityPattern dsp(dof_handler.n_dofs(),dof_handler.n_dofs());


    DoFTools::make_flux_sparsity_pattern (dof_handler, dsp);

    sparsity_pattern.copy_from(dsp);
 
    global_matrix.reinit(sparsity_pattern);   
 
    solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());

}

template<int dim>
void
Base_Solver<dim>::run_ring(const unsigned int refine_cycles)
{
	for (unsigned int i = 0 ; i < refine_cycles ; i ++)
	{
	

		//now we distribute the dofs and allocate the memory for the global matrix. We have 
		// already created the mesh so we can directly distribute the degrees of freedom now.
		distribute_dof_allocate_matrix_ring();

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

		LinearSolver::LinearSolver(global_matrix,system_rhs,solution,LinearSolver::LinearSolver::Pardiso);

		PostProc::Base_PostProc<dim> postproc(constants,base_exactsolution,&dof_handler, &mapping);

		// now we compute the error due to computation
		postproc.error_evaluation_QGauss(solution,this->triangulation.n_active_cells(),this->triangulation.maximum_cell,this->triangulation.diameter());


	}

}