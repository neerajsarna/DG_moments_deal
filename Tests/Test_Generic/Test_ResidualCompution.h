using namespace dealii;

// run the system without any restrictions
TEST(ResidualComputation,HandlesResidualComputation)
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
		for (int cycle = 0 ; cycle < constants.constants_num.refine_cycles ; cycle++)
		{
			base_solver.distribute_dof_allocate_matrix(cycle,constants.constants_num.refine_cycles);

			base_solver.allocate_vectors();
			base_solver.assemble_system_odd();

			// we put a dummy solution
			for (unsigned int i = 0 ; i < base_solver.solution.size() ; i++)
				base_solver.solution(i)= i ;

			Vector<double> residual_dummy;
			residual_dummy = base_solver.matrix_opt.Sparse_matrix_dot_Vector(base_solver.global_matrix,base_solver.solution);
			residual_dummy -= base_solver.system_rhs;


			base_solver.compute_residual();

			for (unsigned int i = 0 ; i < residual_dummy.size() ; i++)
				EXPECT_NEAR(residual_dummy(i),base_solver.VelocitySpace_residual(i),1e-10);

			std::string mat_name = "residual_dummy";
			base_solver.matrix_opt.print_dealii_vector(residual_dummy,mat_name);

			std::string vec_name = "residual";
			base_solver.matrix_opt.print_dealii_vector(base_solver.VelocitySpace_residual,vec_name);
		}

		

}
