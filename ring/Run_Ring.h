template<int dim>
void
Base_Solver<dim>::run_ring()
{
	const unsigned int refine_cycles = constants.refine_cycles;
	error_per_itr.resize(refine_cycles);

	TimerOutput timer (std::cout, TimerOutput::summary,
                   TimerOutput::wall_times);

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
		linear_solver.solve_trilinos(global_matrix,system_rhs,solution);
		timer.leave_subsection();

		timer.enter_subsection("Post Processing");
		PostProc::Base_PostProc<dim> postproc(constants,base_exactsolution,&dof_handler, &mapping);

		// // now we compute the error due to computation
		postproc.error_evaluation_QGauss(solution,this->triangulation.n_active_cells(),
										 error_per_itr[i],
										GridTools::maximal_cell_diameter(this->triangulation),convergence_table);

		postproc.print_options(this->triangulation,solution,i,refine_cycles,convergence_table);		
		timer.leave_subsection();

		// Grid refinement should be done in the end.
		this->refinement_handling(i,refine_cycles);

	}

}


template<int dim>
void
Base_Solver<dim>::run_ring_eigen()
{
	const unsigned int refine_cycles = constants.refine_cycles;
	error_per_itr.resize(refine_cycles);

	for (unsigned int i = 0 ; i < refine_cycles ; i ++)
	{
	
		AssertDimension((int)error_per_itr.size(),constants.refine_cycles);

		//now we distribute the dofs and allocate the memory for the global matrix. We have 
		// already created the mesh so we can directly distribute the degrees of freedom now.
		std::cout << "dof distribution " << std::endl;
		distribute_dof_allocate_matrix_eigen();

		Assert(constants.assembly_type == manuel,ExcMessage("Meshworker cant be used with eigen matrices"));
		Assert(constants.matrix_type == Eigen_Mat, ExcMessage("Algorithm only for an eigen system matrix"));
		// the following routine assembles

		std::cout << "assembly " << std::endl;
		switch (constants.assembly_type)
		{
			case manuel:
			{
				switch(constants.bc_type)
				{
					case characteristic:
					{
						assemble_system_char_eigen();
						break;
					}

					case odd:
					{
						assemble_system_odd_eigen();
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
		

		// we initialize the object which will solve our system
		// We do int the following way so as to keep the solver independent of all the other implementations.
		// This makes the code highly reusable. So one can directly copy the following class and use it somewhere
		// else is one wants to.
		std::cout << "linear solve " << std::endl;
		LinearSolver::LinearSolver linear_solver;
		linear_solver.solve_eigen(global_matrix_eigen,system_rhs,solution);

		std::cout << "Post Proc" << std::endl;
		PostProc::Base_PostProc<dim> postproc(constants,base_exactsolution,&dof_handler, &mapping);

		// // now we compute the error due to computation
		postproc.error_evaluation_QGauss(solution,this->triangulation.n_active_cells(),
										 error_per_itr[i],
										GridTools::maximal_cell_diameter(this->triangulation),convergence_table);

		postproc.print_options(this->triangulation,solution,i,refine_cycles,convergence_table);		


		// Grid refinement should be done in the end.
		this->refinement_handling(i,refine_cycles);

	}

}