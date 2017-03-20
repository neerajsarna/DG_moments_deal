template<int dim>
void
Base_Solver<dim>::run()
{
	const unsigned int refine_cycles = constants.refine_cycles;
	error_per_itr.resize(refine_cycles);

	for (unsigned int i = 0 ; i < refine_cycles ; i ++)
	{
		this->print_mesh_info();

        	TimerOutput timer (std::cout, TimerOutput::summary,
                	   TimerOutput::wall_times);
	
		AssertDimension((int)error_per_itr.size(),constants.refine_cycles);

		//now we distribute the dofs and allocate the memory for the global matrix. We have 
		// already created the mesh so we can directly distribute the degrees of freedom now.

		std::cout << "Distributing dof " << std::endl;
		timer.enter_subsection("Dof Distribution");
		distribute_dof_allocate_matrix();
		allocate_vectors();
		timer.leave_subsection();

		std::cout << "#CELLS " << this->triangulation.n_active_cells() << std::endl;
		std::cout << "#DOFS " << dof_handler.n_dofs() << std::endl;
		std::cout << "Memory by dof handler(Gb) " << dof_handler.memory_consumption()/pow(10,9)<< std::endl;	

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
				//AssertThrow(1 == 0,ExcMessage("Should not have reached here"));
				switch(constants.bc_type)
				{
					case odd:
					{
						std::cout << "Using odd boundary " << std::endl;
						assemble_system_odd();
						break;
					}

					case characteristic:
					default:
					{
						AssertThrow(1 == 0, ExcMessage("Should not have reached here"));
						break;
					}
		
				}
				break;
			}

			default:
			{
				AssertThrow(1==0,ExcMessage("Should not have reached here"));
				break;
			}
		}
		timer.leave_subsection();

		
		// we initialize the object which will solve our system
		// We do int the following way so as to keep the solver independent of all the other implementations.
		// This makes the code highly reusable. So one can directly copy the following class and use it somewhere
		// else if one wants to.

		timer.enter_subsection("Linear Solver");
		LinearSolver::LinearSolver linear_solver;
		std::cout << "Preparing data for pardiso " << std::endl;
		linear_solver.develop_pardiso_data(global_matrix,sparsity_pattern);
		residual = linear_solver.solve_with_pardiso(system_rhs,solution);
		timer.leave_subsection();

		timer.enter_subsection("Post Processing");



		PostProc::Base_PostProc<dim> postproc(constants,base_exactsolution,&dof_handler, &mapping,nEqn,nBC);

		const double residual_weak_form = postproc.compute_residual(residual,this->triangulation.n_active_cells());

		// // now we compute the error due to computation
		postproc.error_evaluation_QGauss(solution,this->triangulation.n_active_cells(),
										 error_per_itr[i],
										GridTools::maximal_cell_diameter(this->triangulation),convergence_table,
										residual_weak_form);



		
		postproc.print_options(this->triangulation,solution,i,refine_cycles,convergence_table,
								system_info[0].base_tensorinfo.S_half_inv,this->finite_element);		
		
		timer.leave_subsection();

		// Grid refinement should be done in the end.
		this->refinement_handling(i,refine_cycles);

	}

}
