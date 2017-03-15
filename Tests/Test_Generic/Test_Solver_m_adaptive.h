using namespace dealii;

TEST(MAdaptiveSolver,HandlesMAdaptiveSolver)
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
		ExactSolution::ExactSolution_Dummy<dim>  exactsolution_dummy(constants.constants_num,System[0].base_tensorinfo.S_half,
			constants.constants_sys.nEqn[0],constants.constants_sys.Ntensors[0]);

		FEM_Solver::Base_Solver<dim> base_solver("grid",
			constants.constants_num,
			System,
			&exactsolution_dummy,
			constants.constants_sys.nEqn,
			constants.constants_sys.nBC);

		base_solver.print_mesh_info();
		fflush(stdout);
		base_solver.distribute_dof_allocate_matrix();


		for (unsigned long int i = 0 ; i < base_solver.nEqn.size(); i++)
			std::cout << "System: " << i << ", Dofs: " << base_solver.dofs_per_cell[i] << std::endl;


		typename hp::DoFHandler<dim>::active_cell_iterator cell = base_solver.dof_handler.begin_active(),
															endc = base_solver.dof_handler.end();


		for (; cell != endc ; cell++)
			std::cout << "center " << cell->center() << " active fe index" << cell->active_fe_index() << std::endl;

		int count = 0;

		// loop over a particular component
		for(int i = 0 ; i < base_solver.nEqn[0] ; i++)

		// loop over all the dofs of that particular component
			for (unsigned int j = 0 ; j < base_solver.dofs_per_component ; j++)
			{
				// compute the dof to which the it corresponds to 
				
				EXPECT_EQ(count,base_solver.finite_element[0].component_to_system_index(i,j));
				count ++;
			}

		count = 0;
		
		// loop over a particular component
		for(int i = 0 ; i < base_solver.nEqn[1] ; i++)

		// loop over all the dofs of that particular component
			for (unsigned int j = 0 ; j < base_solver.dofs_per_component ; j++)
			{
				// compute the dof to which the it corresponds to 
				
				EXPECT_EQ(count,base_solver.finite_element[1].component_to_system_index(i,j));
				count ++;
			}

		base_solver.assemble_system_odd();


}