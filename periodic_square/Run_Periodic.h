// distribution of dof but for an eigen matrix
template<int dim>
void 
Base_Solver<dim>::distribute_dof_allocate_matrix_periodic_box_eigen()
{
    dof_handler.distribute_dofs(finite_element);

    DynamicSparsityPattern dsp(dof_handler.n_dofs(),dof_handler.n_dofs());


    DoFTools::make_flux_sparsity_pattern (dof_handler, dsp);
  
    this->add_periodic_sparsity(dsp);
 
 	global_matrix_eigen.resize(dof_handler.n_dofs(),dof_handler.n_dofs());
 	global_matrix_eigen.setZero();
    global_matrix_eigen(dsp.n_nonzero_elements());   
 
    solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());

}

template<int dim>
void 
Base_Solver<dim>::distribute_dof_allocate_matrix_periodic_box()
{
    dof_handler.distribute_dofs(finite_element);

    DynamicSparsityPattern dsp(dof_handler.n_dofs(),dof_handler.n_dofs());


    DoFTools::make_flux_sparsity_pattern (dof_handler, dsp);
  
    this->add_periodic_sparsity(dsp);

//    sparsity_pattern.copy_from(dsp);
 
    global_matrix.reinit(dsp);   
 
  //  solution.reinit (dof_handler.n_dofs());
  //  system_rhs.reinit (dof_handler.n_dofs());

}

template<int dim>
void
Base_Solver<dim>::run_periodic()
{

	error_per_itr.resize(constants.refine_cycles);

	TimerOutput timer (std::cout, TimerOutput::summary,
                   TimerOutput::wall_times);
	
	Assert(constants.mesh_type == periodic_square,ExcMessage("Incorrect mesh"));
	for (int i = 0 ; i < constants.refine_cycles ; i++)
	{

		AssertDimension((int)error_per_itr.size(),constants.refine_cycles);

		timer.enter_subsection("Develop Periodicity");
		// first we develop the periodic faces using internal functions of dealii
		this->develop_periodic_faces(dof_handler);

		// now we construct the required data structure
		this->divide_periodicity();

		timer.leave_subsection();

		//now we distribute the dofs and allocate the memory for the global matrix. We have 
		// already created the mesh so we can directly distribute the degrees of freedom now.
		timer.enter_subsection("Distribute Dof");
		distribute_dof_allocate_matrix_periodic_box();
		allocate_vectors();
		timer.leave_subsection();



		// cannot use meshworker for the present problem
		timer.enter_subsection("Assemble");
		Assert(constants.assembly_type == meshworker,ExcNotImplemented());
		Assert(constants.matrix_type == Trilinos_Mat,ExcMessage("Routine only for Trilinos matrix"));
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

		timer.enter_subsection("Post Proc");
		PostProc::Base_PostProc<dim> postproc(constants,base_exactsolution,&dof_handler, &mapping);

		residual_strong_form = postproc.compute_residual(solution,finite_element,system_info,this->triangulation.n_active_cells()); 
		std::cout << "********Residual*********: " << residual_strong_form << std::endl;

		// now we compute the error due to computation
		postproc.error_evaluation_QGauss(solution,this->triangulation.n_active_cells(),
										 error_per_itr[i],
										GridTools::maximal_cell_diameter(this->triangulation),
										convergence_table,residual_strong_form);

		postproc.print_options(this->triangulation,solution,i,constants.refine_cycles,convergence_table);	

		timer.leave_subsection();

		// now we refine the grid by changing the number of parts in the y direction
		//But we do not want to wast time in refining the grid in the last refinement cycle
		if (i < constants.refine_cycles - 1)
			this->mesh_internal_periodic_square(constants.part_x,constants.part_y + 100 * (i + 1));

	}
 }

template<int dim>
void
Base_Solver<dim>::run_periodic_unsymmetric()
{

	error_per_itr.resize(constants.refine_cycles);
	
	for (int i = 0 ; i < constants.refine_cycles ; i++)
	{

		AssertDimension((int)error_per_itr.size(),constants.refine_cycles);

		// first we develop the periodic faces using internal functions of dealii
		this->develop_periodic_faces(dof_handler);

		// now we construct the required data structure
		this->divide_periodicity();


		//now we distribute the dofs and allocate the memory for the global matrix. We have 
		// already created the mesh so we can directly distribute the degrees of freedom now.
		distribute_dof_allocate_matrix_periodic_box();



		// cannot use meshworker for the present problem
		Assert(constants.assembly_type == meshworker,ExcNotImplemented());
		Assert(constants.matrix_type == Trilinos_Mat,ExcMessage("Routine only for Trilinos matrix, Change the matrix type in the input file"));
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

		postproc.compute_error_unsymmetric(this->triangulation,
										   this->solution,
										   system_info->base_tensorinfo.S_half_inv,
										   error_per_itr[i]);

		// now we refine the grid by changing the number of parts in the y direction
		//But we do not want to wast time in refining the grid in the last refinement cycle
		if (i < constants.refine_cycles - 1)
			this->mesh_internal_periodic_square(constants.part_x,constants.part_y + 50 * (i + 1));

	}


}



template<int dim>
void
Base_Solver<dim>::run_periodic_eigen()
{

	error_per_itr.resize(constants.refine_cycles);

	TimerOutput timer (std::cout, TimerOutput::summary,
                   TimerOutput::wall_times);
	
	for (int i = 0 ; i < constants.refine_cycles ; i++)
	{

		AssertDimension((int)error_per_itr.size(),constants.refine_cycles);

		timer.enter_subsection("Develop Periodicity");
		// first we develop the periodic faces using internal functions of dealii
		this->develop_periodic_faces(dof_handler);

		// now we construct the required data structure
		this->divide_periodicity();

		timer.leave_subsection();

		//now we distribute the dofs and allocate the memory for the global matrix. We have 
		// already created the mesh so we can directly distribute the degrees of freedom now.
		timer.enter_subsection("Distribute Dof");
		distribute_dof_allocate_matrix_periodic_box_eigen();
		timer.leave_subsection();



		// cannot use meshworker for the present problem
		timer.enter_subsection("Assemble");
		Assert(constants.assembly_type == meshworker,ExcNotImplemented());
		Assert(constants.matrix_type == Eigen_Mat,ExcMessage("Routine only for Trilinos matrix"));
		switch(constants.bc_type)
		{
			case characteristic:
			{
				assemble_system_periodic_char_eigen();
				break;
			}

			case odd:
			{
				assemble_system_periodic_odd_eigen();
				break;
			}
		}
		timer.leave_subsection();

        
		// we initialize the object which will solve our system
		// We do int the following way so as to keep the solver independent of all the other implementations.
		// This makes the code highly reusable. So one can directly copy the following class and use it somewhere
		// else is one wants to.
		timer.enter_subsection("System Solver");
		LinearSolver::LinearSolver linear_solver;
		linear_solver.solve_eigen(global_matrix_eigen,system_rhs,solution);
		timer.leave_subsection();

		timer.enter_subsection("Post Proc");
		PostProc::Base_PostProc<dim> postproc(constants,base_exactsolution,&dof_handler, &mapping);

		// now we compute the error due to computation
		postproc.error_evaluation_QGauss(solution,this->triangulation.n_active_cells(),
										 error_per_itr[i],
										GridTools::maximal_cell_diameter(this->triangulation),convergence_table);

		postproc.print_options(this->triangulation,solution,i,constants.refine_cycles,convergence_table);	

		timer.leave_subsection();

		// now we refine the grid by changing the number of parts in the y direction
		//But we do not want to wast time in refining the grid in the last refinement cycle
		if (i < constants.refine_cycles - 1)
			this->mesh_internal_periodic_square(constants.part_x,constants.part_y + 100 * (i + 1));

	}
 }

template<int dim>
void
Base_Solver<dim>::run_periodic_unsymmetric_eigen()
{

	error_per_itr.resize(constants.refine_cycles);
	
	for (int i = 0 ; i < constants.refine_cycles ; i++)
	{

		AssertDimension((int)error_per_itr.size(),constants.refine_cycles);

		// first we develop the periodic faces using internal functions of dealii
		this->develop_periodic_faces(dof_handler);

		// now we construct the required data structure
		this->divide_periodicity();


		//now we distribute the dofs and allocate the memory for the global matrix. We have 
		// already created the mesh so we can directly distribute the degrees of freedom now.
		distribute_dof_allocate_matrix_periodic_box();



		// cannot use meshworker for the present problem
		Assert(constants.assembly_type == meshworker,ExcNotImplemented());
		Assert(constants.matrix_type == Eigen_Mat,ExcMessage("Routine only for Eigen matrix, Change the matrix type in the input file"));
		switch(constants.bc_type)
		{
			case characteristic:
			{
				assemble_system_periodic_char_eigen();
				break;
			}

			case odd:
			{
				assemble_system_periodic_odd_eigen();
				break;
			}
		}


		// we initialize the object which will solve our system
		// We do int the following way so as to keep the solver independent of all the other implementations.
		// This makes the code highly reusable. So one can directly copy the following class and use it somewhere
		// else is one wants to.
		LinearSolver::LinearSolver linear_solver(global_matrix_eigen,system_rhs,solution);

		PostProc::Base_PostProc<dim> postproc(constants,base_exactsolution,&dof_handler, &mapping);

		postproc.compute_error_unsymmetric(this->triangulation,
										   this->solution,
										   system_info->base_tensorinfo.S_half_inv,
										   error_per_itr[i]);

		// now we refine the grid by changing the number of parts in the y direction
		//But we do not want to wast time in refining the grid in the last refinement cycle
		if (i < constants.refine_cycles - 1)
			this->mesh_internal_periodic_square(constants.part_x,constants.part_y + 50 * (i + 1));

	}


}

