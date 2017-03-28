using namespace dealii;

// TEST(SetTolerance,HandlesToleranceSetting)
// {
// 		const unsigned int dim = 1;

// 		std::string folder_name = "../system_matrices/";
// 		Constants::Base_Constants constants(input_file);

// 	// 	// we construct a vector of all the system data being considered 
// 		std::vector<Develop_System::System<dim>> System;

// 		// initialize the vector containing all the systems
// 		// initialize the vector containing all the systems
// 		for(int i = 0 ; i < constants.constants_sys.total_systems ; i++)
// 		{
// 			System.push_back(Develop_System::System<dim>(constants.constants_num,constants.constants_sys.nEqn[i],
// 				constants.constants_sys.nBC[i],constants.constants_sys.Ntensors[i],folder_name));

// 		}
	
// 		for (int i = 0 ; i < constants.constants_sys.total_systems ; i++)
// 			System[i].initialize_system();
		
		

// 		// the exact solution can only be created for one of the systems
// 		ExactSolution::PoissonHeat<dim>  PoissonHeat(constants.constants_num,System[constants.constants_sys.total_systems-1].base_tensorinfo.S_half,
// 													constants.constants_sys.nEqn[constants.constants_sys.total_systems-1],constants.constants_sys.Ntensors[constants.constants_sys.total_systems-1]);

// 		FEM_Solver::Base_Solver<dim> base_solver("grid",
// 											 	 constants.constants_num,
// 											 	 System,
// 											 	 &PoissonHeat,
// 											 	 constants.constants_sys.nEqn,
// 											 	 constants.constants_sys.nBC);

// 		base_solver.print_mesh_info();

// 		base_solver.VelocitySpace_error_per_cell.reinit(base_solver.triangulation.n_active_cells());

// 		for (unsigned int i = 0 ; i < base_solver.triangulation.n_active_cells() ; i++)
// 			base_solver.VelocitySpace_error_per_cell(i) = i;

// 		base_solver.set_tolerance_bands();

// 		for (unsigned long int i = 0 ; i < base_solver.VelocitySpace_error_tolerance.size() ; i++)
// 			std::cout << "Tolerance: " <<  base_solver.VelocitySpace_error_tolerance[i] << std::endl;
// }

// // run the system without any restrictions
TEST(Solver,HandlesSolver)
{
		const unsigned int dim = 2;

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
		
		
		if (dim == 2)
		{
		// the exact solution can only be created for one of the systems
		ExactSolution::ExactSolution_Dummy<dim>  exactsolution_dummy(constants.constants_num,System[constants.constants_sys.total_systems-1].base_tensorinfo.S_half,
																			 constants.constants_sys.nEqn[constants.constants_sys.total_systems-1],constants.constants_sys.Ntensors[constants.constants_sys.total_systems-1]);

		FEM_Solver::Base_Solver<dim> base_solver("grid",
											 	 constants.constants_num,
											 	 System,
											 	 &exactsolution_dummy,
											 	 constants.constants_sys.nEqn,
											 	 constants.constants_sys.nBC);

		base_solver.print_mesh_info();
		



		base_solver.error_per_itr.resize(constants.constants_num.refine_cycles);
		for (int cycle = 0 ; cycle < constants.constants_num.refine_cycles ; cycle++)
		{
			base_solver.distribute_dof_allocate_matrix(cycle,constants.constants_num.refine_cycles);

			base_solver.allocate_vectors();
			base_solver.assemble_system_odd();
		
			LinearSolver::LinearSolver linear_solver;
			linear_solver.develop_pardiso_data(base_solver.global_matrix);
			base_solver.residual = linear_solver.solve_with_pardiso(base_solver.system_rhs,base_solver.solution);


			PostProc::Base_PostProc<dim> postproc(base_solver.constants,&exactsolution_dummy,base_solver.nEqn,base_solver.nBC);
			postproc.reinit(base_solver.dof_handler);


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


			base_solver.compute_equilibrium_deviation();

			// we now set the tolerance
			// we do not set the tolerance for the last iteration
			if (cycle != constants.constants_num.refine_cycles -1)
				base_solver.set_tolerance_bands();

		

			FILE *fp;
			std::string filename;

			filename = "residual_per_cell_" + std::to_string(cycle);
			fp = fopen(filename.c_str(),"w+");

			AssertThrow(fp != NULL,ExcMessage("could not open the file"));
			
			typename hp::DoFHandler<dim>::active_cell_iterator cell = base_solver.dof_handler.begin_active(),
													  endc = base_solver.dof_handler.end();


			unsigned int counter  = 0;

			for (; cell != endc; cell++)
			{
				if (dim == 1)
					fprintf(fp, "%f %0.30f %u\n",cell->center()(0),base_solver.VelocitySpace_error_per_cell(counter),cell->active_fe_index());

				if (dim == 2)
					fprintf(fp, "%f %f %f %u\n",cell->center()(0),cell->center()(1),base_solver.VelocitySpace_error_per_cell(counter),cell->active_fe_index());

				counter++;
			}

			fclose(fp);		




	

		}			
		}

		if (dim == 1)
		{
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
		



		base_solver.error_per_itr.resize(constants.constants_num.refine_cycles);
		for (int cycle = 0 ; cycle < constants.constants_num.refine_cycles ; cycle++)
		{
			base_solver.distribute_dof_allocate_matrix(cycle,constants.constants_num.refine_cycles);

			base_solver.allocate_vectors();
			base_solver.assemble_system_odd();
		
			LinearSolver::LinearSolver linear_solver;
			linear_solver.develop_pardiso_data(base_solver.global_matrix);
			base_solver.residual = linear_solver.solve_with_pardiso(base_solver.system_rhs,base_solver.solution);


			PostProc::Base_PostProc<dim> postproc(base_solver.constants,&PoissonHeat,base_solver.nEqn,base_solver.nBC);
			postproc.reinit(base_solver.dof_handler);


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


			base_solver.compute_equilibrium_deviation();

			// we now set the tolerance
			// we do not set the tolerance for the last iteration
			if (cycle != constants.constants_num.refine_cycles -1)
				base_solver.set_tolerance_bands();

		

			FILE *fp;
			std::string filename;

			filename = "residual_per_cell_" + std::to_string(cycle);
			fp = fopen(filename.c_str(),"w+");

			AssertThrow(fp != NULL,ExcMessage("could not open the file"));
			
			typename hp::DoFHandler<dim>::active_cell_iterator cell = base_solver.dof_handler.begin_active(),
													  endc = base_solver.dof_handler.end();


			unsigned int counter  = 0;

			for (; cell != endc; cell++)
			{
				if (dim == 1)
					fprintf(fp, "%f %0.30f %u\n",cell->center()(0),base_solver.VelocitySpace_error_per_cell(counter),cell->active_fe_index());

				if (dim == 2)
					fprintf(fp, "%f %f %f %u\n",cell->center()(0),cell->center()(1),base_solver.VelocitySpace_error_per_cell(counter),cell->active_fe_index());

				counter++;
			}

			fclose(fp);		




	

		}
	}

	
}