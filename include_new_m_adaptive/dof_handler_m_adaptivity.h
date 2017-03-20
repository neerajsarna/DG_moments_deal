template<int dim>
void 
Base_Solver<dim>::construct_block_structure(std::vector<int> &block_structure,const std::vector<int> &nEqn)
{
		Assert(nEqn.size() != 0, ExcNotInitialized());
		Assert(std::is_sorted(std::begin(nEqn),std::end(nEqn)),ExcMessage("number of equations not sorted"));

		// the very first entry should be the number of equations in the first system
		block_structure.push_back(nEqn[0]);

		for (unsigned long int i = 1 ; i < nEqn.size() ; i++)
			block_structure.push_back(nEqn[i]-nEqn[i-1]);

		AssertDimension(block_structure.size(),nEqn.size());
}

template<int dim>
void 
Base_Solver<dim>::construct_fe_collection()
{
		std::vector<int> block_structure;

		// first create the block structure for the finite element object which will then be used to construct the fe system
		construct_block_structure(block_structure,nEqn);

		if (block_structure.size() == 1)
			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0]));

		if (block_structure.size() == 2)
		{

			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
													FE_Nothing<dim>(),block_structure[1]));

			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
													fe_dg,block_structure[1]));
		}

		if (block_structure.size() == 3)
		{
			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
													FE_Nothing<dim>(),block_structure[1],
													FE_Nothing<dim>(),block_structure[2]));

			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
												   fe_dg,block_structure[1],
												   FE_Nothing<dim>(),block_structure[2]));

			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
												   fe_dg,block_structure[1],
												   fe_dg,block_structure[2]));
		}

		if (block_structure.size() == 4)
		{
			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
												   FE_Nothing<dim>(),block_structure[1],
							 						FE_Nothing<dim>(),block_structure[2],
							 						FE_Nothing<dim>(),block_structure[3]));

			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
				                                   fe_dg,block_structure[1],
							                       FE_Nothing<dim>(),block_structure[2],
							                       FE_Nothing<dim>(),block_structure[3]));

			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
												   fe_dg,block_structure[1],
							 					   fe_dg,block_structure[2],
							 					   FE_Nothing<dim>(),block_structure[3]));

			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
				                                   fe_dg,block_structure[1],
							 					   fe_dg,block_structure[2],
							 					   fe_dg,block_structure[3]));
		}
}

template<>
void 
Base_Solver<1>::allocate_fe_index()
{
	typename hp::DoFHandler<1>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();

	// fraction of the half of the domain which will receive lower order moments
	const double domain_adapt = 2.0;

	// has not been implemented for more than two systems
	AssertThrow(nEqn.size() == 2 || nEqn.size() == 1 || nEqn.size() == 4,ExcNotImplemented());

	switch(nEqn.size())
	{
		case 1:
		{
			for (; cell != endc ; cell++)
				cell->set_active_fe_index(0);		

			break;			
		}

		case 2:
		{
			for (; cell != endc ; cell++)
			{
				// // higher order moment theory near the boundary
				if (cell->center()(0) >= domain_adapt * 0.5 || cell->center()(0) <= -domain_adapt * 0.5)
					cell->set_active_fe_index(1);

				// lower order moment theory towards the interior
				if (cell->center()(0) > -domain_adapt * 0.5 && cell->center()(0) < domain_adapt * 0.5)					
					cell->set_active_fe_index(0);		

				// we now consider two cells. one of them gets
				// Assert(cell->index() <= 1 ,ExcNotImplemented());
				// if (cell->index() == 0 )
				// 	cell->set_active_fe_index(0);

				// if (cell->index() == 1)
				// 	cell->set_active_fe_index(1);

			}

			break;
		}


		case 4:
		{
			for (; cell != endc ; cell++)
			{
				// // higher order moment theory near the boundary
				// if (cell->center()(0) > domain_adapt * 0.5 || cell->center()(0) < -domain_adapt * 0.5)
				// 	cell->set_active_fe_index(1);

				// // lower order moment theory towards the interior
				// if (cell->center()(0) >= -domain_adapt * 0.5 && cell->center()(0) <= domain_adapt * 0.5)					
				// 	cell->set_active_fe_index(0);		

				// we now consider two cells. one of them gets
				Assert(cell->index() <= 3 ,ExcNotImplemented());
				if (cell->index() == 0 )
					cell->set_active_fe_index(0);

				if (cell->index() == 1)
					cell->set_active_fe_index(1);

				if (cell->index() == 2)
					cell->set_active_fe_index(2);

				if (cell->index() == 3)
					cell->set_active_fe_index(3);
			}

			break;
		}

		default:
		{
			AssertThrow(1 == 0,ExcMessage("Should not have reached here"));
			break;
		}

	}

}


// allocation of finite element objects for a 2D problem
template<>
void 
Base_Solver<2>::allocate_fe_index()
{
	typename hp::DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();


	// has not been implemented for more than 1 systems
	AssertDimension(nEqn.size(),1);

	for (; cell != endc ; cell++)
		cell->set_active_fe_index(0);

}