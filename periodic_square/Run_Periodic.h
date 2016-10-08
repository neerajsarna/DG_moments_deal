template<int dim>
void 
Base_Solver<dim>::distribute_dof_allocate_matrix_periodic_box()
{
    dof_handler.distribute_dofs(finite_element);

    DynamicSparsityPattern dsp(dof_handler.n_dofs(),dof_handler.n_dofs());


    DoFTools::make_flux_sparsity_pattern (dof_handler, dsp);
  
    this->add_periodic_sparsity(dsp);

    sparsity_pattern.copy_from(dsp);
 
    global_matrix.reinit(sparsity_pattern);   
 
    solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());

}

template<int dim>
void
Base_Solver<dim>::run_periodic()
{
	
		

		// first we develop the periodic faces using internal functions of dealii
		this->develop_periodic_faces(dof_handler);

		// now we construct the required data structure
		this->divide_periodicity();


		//now we distribute the dofs and allocate the memory for the global matrix. We have 
		// already created the mesh so we can directly distribute the degrees of freedom now.
		distribute_dof_allocate_matrix_periodic_box();



		// cannot use meshworker for the present problem
		switch(constants.bc_type)
		{
			case characteristic:
			{
				assemble_system_periodic_char();
				break;
			}

			case odd:
			{
				assemble_system_periodic_odd();
				break;
			}
		}


		// we initialize the object which will solve our system
		// We do int the following way so as to keep the solver independent of all the other implementations.
		// This makes the code highly reusable. So one can directly copy the following class and use it somewhere
		// else is one wants to.
		LinearSolver::LinearSolver linear_solver(global_matrix,system_rhs,solution);

		PostProc::Base_PostProc<dim> postproc(constants,base_exactsolution,&dof_handler, &mapping);

		const double L2_error = postproc.L2_error_QGauss(solution,this->triangulation.n_active_cells());

		const double Linf_error = postproc.Linfty_error_QGauss(solution,this->triangulation.n_active_cells());

		std::cout << "\n>>>>>>>>>>>>>>>>>>Error Details<<<<<<<<<<<<<<<<<<<\n" << std::endl;
		std::cout << "L2_error " << L2_error << " Linf_error " << Linf_error << "\n" << std::endl;

		// now we print the convergence table
		postproc.print_options(this->triangulation,solution,0,constants.refine_cycles,L2_error,Linf_error,this->triangulation.n_active_cells(),
								GridTools::maximal_cell_diameter(this->triangulation));		


}