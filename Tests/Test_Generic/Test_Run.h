using namespace dealii;

TEST(Solver,HandlesSolver)
{
		const unsigned int dim = 1;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);

	// 	// we construct a vector of all the system data being considered 
		std::vector<Develop_System::System<dim>> System;

		// initialize the vector containing all the systems
		// initialize the vector containing all the systems
		for(int i = 0 ; i < constants.constants_sys.total_systems ; i++)
		{
			System.push_back(Develop_System::System<dim>(constants.constants_num,constants.constants_sys.nEqn[i],
				constants.constants_sys.nBC[i],constants.constants_sys.Ntensors[i],folder_name));

		}
	
		for (int i = 0 ; i < constants.constants_sys.total_systems ; i++)
			System[i].initialize_system();
		

		// the tests have not been implemented for more than one system
		AssertDimension(constants.constants_sys.total_systems,3);
		

				// the exact solution can only be created for one of the systems
				ExactSolution::PoissonHeat<dim>  PoissonHeat(constants.constants_num,System[constants.constants_sys.total_systems-1].base_tensorinfo.S_half,
													constants.constants_sys.nEqn[constants.constants_sys.total_systems-1],constants.constants_sys.Ntensors[constants.constants_sys.total_systems-1]);

				FEM_Solver::Base_Solver<dim> base_solver("grid",
											 	 constants.constants_num,
											 	 System,
											 	 &PoissonHeat,
											 	 constants.constants_sys.nEqn,
											 	 constants.constants_sys.nBC);

		base_solver.print_mesh_info();
		// tolerance per cell
		const double tolerance = 1.955894e-03/base_solver.triangulation.n_active_cells();

		base_solver.error_per_itr.resize(constants.constants_num.refine_cycles);
		for (unsigned int cycle = 0 ; cycle < constants.constants_num.refine_cycles ; cycle++)
		{
			base_solver.distribute_dof_allocate_matrix(base_solver.error_per_cell_VelocitySpace,
													  tolerance,
													  cycle,
													  constants.constants_num.refine_cycles);

			base_solver.allocate_vectors();
			base_solver.assemble_system_odd();

		
			LinearSolver::LinearSolver linear_solver;
			linear_solver.develop_pardiso_data(base_solver.global_matrix);
			base_solver.residual = linear_solver.solve_with_pardiso(base_solver.system_rhs,base_solver.solution);


			PostProc::Base_PostProc<dim> postproc(base_solver.constants,&PoissonHeat,base_solver.nEqn,base_solver.nBC);
			postproc.reinit(base_solver.dof_handler);

			// // now we compute the error due to computation
			base_solver.error_per_cell_VelocitySpace =  postproc.return_error_per_cell(base_solver.solution,
										 						base_solver.triangulation.n_active_cells(),
										 				base_solver.mapping,base_solver.dof_handler);

		// // now we compute the error due to computation
			postproc.error_evaluation_QGauss(base_solver.solution,
										 	base_solver.triangulation.n_active_cells(),
										   base_solver.error_per_itr[cycle],
										GridTools::maximal_cell_diameter(base_solver.triangulation),
										base_solver.convergence_table,
										0,base_solver.mapping,base_solver.dof_handler);

			postproc.print_options(base_solver.triangulation,base_solver.solution,cycle,constants.constants_num.refine_cycles,
								base_solver.convergence_table,
								base_solver.system_info[constants.constants_sys.total_systems-1].base_tensorinfo.S_half_inv,
								base_solver.dof_handler);

			base_solver.refinement_handling(cycle,constants.constants_num.refine_cycles);		

			std::cout << "Error: " << base_solver.error_per_cell_VelocitySpace.l2_norm() << std::endl;
		}

		typename hp::DoFHandler<1>::active_cell_iterator cell = base_solver.dof_handler.begin_active(),
													     endc = base_solver.dof_handler.end();
		

}