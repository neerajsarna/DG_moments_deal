namespace Test_Solver
{
	using namespace dealii;



	// TEST(ChecksPeriodicGrid,HandlesPeridicGrid)
	// {
	// 	const unsigned int dim = 2;
	// 	ASSERT_EQ(dim,2) << "3D not implemented" << std::endl;

	// 	std::string folder_name = "../system_matrices/";
	// 	Constants::Base_Constants constants(input_file);
	// 	G20::G20<dim> G20(constants.constants,folder_name);

	// 	ExactSolution::G20_PoissonHeat<dim>  G20_PoissonHeat(constants.constants,G20.base_tensorinfo.S_half);
	// 	FEM_Solver::Base_Solver<dim> base_solver("mesh_file_name","grid",
	// 											constants.constants,&G20,
	// 										 	&G20_PoissonHeat);


	// 	if (constants.constants.mesh_type == periodic_square)
	// 	{
	// 		base_solver.print_grid(0);

	// 		// triangulation should have some number of cell
	// 		EXPECT_NE(base_solver.triangulation.n_active_cells(),0);

	// 	}


	// }

	// TEST(PeriodicityHandler,HandlesPeriodicityManagement)
	// {
	// 	const unsigned int dim = 2;
	// 	ASSERT_EQ(dim,2) << "3D not implemented" << std::endl;

	// 	std::string folder_name = "../system_matrices/";
	// 	Constants::Base_Constants constants(input_file);
	// 	G26::G26<dim> G26(constants.constants,folder_name);

	// 	ExactSolution::G26_PoissonHeat<dim>  G26_PoissonHeat(constants.constants,G26.base_tensorinfo.S_half);

	// 	FEM_Solver::Base_Solver<dim> base_solver("mesh_file_name",
	// 										 "grid",
	// 										 constants.constants,
	// 										 &G26,
	// 										 &G26_PoissonHeat);


	// 	// first we check the boundary IDs of the periodic mesh
	// 	typename Triangulation<dim>::active_cell_iterator cell = base_solver.triangulation.begin_active(),
	// 													   endc = base_solver.triangulation.end();

	// 	if (constants.constants.mesh_type == periodic_square)
	// 	{

	// 	// now we check the routines of periodicity handler
	// 	FESystem<dim> finite_element(FE_DGQ<dim>(2),constants.constants.nEqn);
	//     DoFHandler<dim> dof_handler(base_solver.triangulation);
	//     dof_handler.distribute_dofs(finite_element);

	// 	base_solver.develop_periodic_faces(dof_handler);

	// 	DynamicSparsityPattern dsp(dof_handler.n_dofs(),dof_handler.n_dofs());
	// 	base_solver.add_periodic_sparsity(dsp);

	// 	std::vector<types::global_dof_index> local_dof_cell1;
 //    	std::vector<types::global_dof_index> local_dof_cell2;

	// 	// now we check whether the periodicity vector is okay or not
	// 	for (int i = 0 ; i < (int)base_solver.periodicity_vector.size() ; i++)
	// 	{
	// 		    const typename DoFHandler<dim>::cell_iterator cell_1 = base_solver.periodicity_vector[i].cell[0];
 //    			const typename DoFHandler<dim>::cell_iterator cell_2 = base_solver.periodicity_vector[i].cell[1];

 //    			const unsigned int dof_cell1 = cell_1->get_fe().dofs_per_cell;
 //    			const unsigned int dof_cell2 = cell_2->get_fe().dofs_per_cell;

 //    			local_dof_cell1.resize(dof_cell1);
 //    			local_dof_cell2.resize(dof_cell2);

 //    			cell_1->get_dof_indices(local_dof_cell1);
 //    			cell_2->get_dof_indices(local_dof_cell2);


 //    			for (unsigned int i = 0 ; i < dof_cell1 ; i++)
 //    				for (unsigned int j = 0 ; j < dof_cell2 ; j++)
 //    				{
 //    					const bool check1 = dsp.exists(local_dof_cell1[i],local_dof_cell2[j]);
 //    					const bool check2 = dsp.exists(local_dof_cell2[j],local_dof_cell1[i]);

 //    					if (!check1)
 //    						ASSERT_EQ(1,0)<<"Incorrect sparsity pattern" << std::endl;

 //    					if (!check2)
 //    						ASSERT_EQ(1,0) << "Incorrect sparsity pattern" << std::endl;
 //    				}
 //    			// the y coordinate should be equal
 //    			ASSERT_NEAR(cell_1->center()(1),cell_2->center()(1),1e-5) << "Incorrect periodicity vector ";
	// 	}

	// 	base_solver.divide_periodicity();
	// 	//check on the right boundary
	// 	if (constants.constants.part_x == 2 && constants.constants.part_y == 2)
	// 	{
	// 		const double x1 = 1.0;
	// 		const double y1 = 0.5;

	// 		typename DoFHandler<dim>::active_cell_iterator it = base_solver.get_periodic_neighbor(x1,y1);

	// 		EXPECT_NEAR(it->center()(1),y1,1e-5) << "problem with finding neighbor " ;

	// 	//check on the left boundary
	// 		const double x2 = -1.0;
	// 		const double y2 = -0.5;

	// 		it = base_solver.get_periodic_neighbor(x2,y2);

	// 		EXPECT_NEAR(it->center()(1),y2,1e-5) << "problem with finding neighbor " ;			
	// 	}

	// 	}
	// }

	// understanding the grid refinement for a periodic square
	// TEST(GridRefinementPeriodicSquare,HandlesRefinement)
	// {
	// 	const unsigned int dim = 2;
	// 	ASSERT_EQ(dim,2) << "3D not implemented" << std::endl;

	// 	std::string folder_name = "../system_matrices/";
	// 	Constants::Base_Constants constants(input_file);
	// 	G20::G20<dim> G20(constants.constants,folder_name);

	// 	ExactSolution::G20_PoissonHeat<dim>  G20_PoissonHeat(constants.constants,G20.base_tensorinfo.S_half);

	// 	FEM_Solver::Base_Solver<dim> base_solver("mesh_file_name",
	// 											"grid",
	// 											constants.constants,
	// 											&G20,
	// 											&G20_PoissonHeat);		

	// 	for (int i = 0 ; i < constants.constants.refine_cycles ; i++)
	// 	{

	// 	// first we develop the periodic faces using internal functions of dealii
	// 		base_solver.develop_periodic_faces(base_solver.dof_handler);

	// 	// now we construct the required data structure
	// 		base_solver.divide_periodicity();

			
	// 		//The first entry and the ycoord of the second entry should match
	// 			for (std::map<double, typename DoFHandler<dim>::cell_iterator>::iterator it = base_solver.xplus1_set.begin();
	// 					 it != base_solver.xplus1_set.end() ;it++)
	// 				EXPECT_NEAR(it->first,it->second->center()(1),1e-8);


	// 			for (std::map<double, typename DoFHandler<dim>::cell_iterator>::iterator it = base_solver.xminus1_set.begin();
	// 					 it != base_solver.xminus1_set.end() ;it++)
	// 				EXPECT_NEAR(it->first,it->second->center()(1),1e-8);


	// 		// check the ycoordinates of the periodic pairs
	// 		for (int j = 0 ; j < (int)base_solver.periodicity_vector.size(); j++)
	// 		{
	// 			const typename DoFHandler<dim>::cell_iterator cell_1 = base_solver.periodicity_vector[i].cell[0];
 //        		const typename DoFHandler<dim>::cell_iterator cell_2 = base_solver.periodicity_vector[i].cell[1];

 //        		// the ycoordinates of the cells should be the same
 //        		EXPECT_NEAR(cell_1->center()(1),cell_2->center()(1),1e-5);
	// 		}

	// 		if (i < constants.constants.refine_cycles - 1)
	// 			base_solver.mesh_internal_periodic_square(constants.constants.part_x,constants.constants.part_y + 50 * (i + 1));

	// 	}
	// }



	// TEST(PoissonHeatG20,HandlesPoissonHeatG20)
	// {
	// 	const unsigned int dim = 2;
	// 	ASSERT_EQ(dim,2) << "3D not implemented" << std::endl;

	// 	std::string folder_name = "../system_matrices/";
	// 	Constants::Base_Constants constants(input_file);
	// 	G20::G20<dim> G20(constants.constants,folder_name);

	// 	ExactSolution::G20_PoissonHeat<dim>  G20_PoissonHeat(constants.constants,G20.base_tensorinfo.S_half);

	// 	FEM_Solver::Base_Solver<dim> base_solver("mesh_file_name",
	// 										 "grid",
	// 										 constants.constants,
	// 										 &G20,
	// 										 &G20_PoissonHeat);


	// 	base_solver.run_periodic();

	// }

	TEST(SquareCavityG20,SquareCavityG20)
	{
		const unsigned int dim = 2;
		ASSERT_EQ(dim,2) << "3D not implemented" << std::endl;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		G20::G20<dim> G20(constants.constants,folder_name);

		ExactSolution::G20_PoissonHeat<dim>  G20_PoissonHeat(constants.constants,G20.base_tensorinfo.S_half);

		FEM_Solver::Base_Solver<dim> base_solver("mesh_file_name",
											 "grid",
											 constants.constants,
											 &G20,
											 &G20_PoissonHeat);


		base_solver.run_square();

	}
}
