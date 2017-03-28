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

		switch (block_structure.size())
		{
			case 1:
			{
				finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0]));
				break;
			}

			case 2:
			{
			
				finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
													FE_Nothing<dim>(),block_structure[1]));

				finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
													fe_dg,block_structure[1]));

				break;
			}

			case 3:
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
				break;
			}

			case 4:
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
				break;
			}

			case 5:
			{
			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
												   FE_Nothing<dim>(),block_structure[1],
							 						FE_Nothing<dim>(),block_structure[2],
							 						FE_Nothing<dim>(),block_structure[3],
							 						FE_Nothing<dim>(),block_structure[4]));

			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
												   fe_dg,block_structure[1],
							 						FE_Nothing<dim>(),block_structure[2],
							 						FE_Nothing<dim>(),block_structure[3],
							 						FE_Nothing<dim>(),block_structure[4]));

			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
												   fe_dg,block_structure[1],
							 						fe_dg,block_structure[2],
							 						FE_Nothing<dim>(),block_structure[3],
							 						FE_Nothing<dim>(),block_structure[4]));

			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
												   fe_dg,block_structure[1],
							 						fe_dg,block_structure[2],
							 						fe_dg,block_structure[3],
							 						FE_Nothing<dim>(),block_structure[4]));

		finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
												   fe_dg,block_structure[1],
							 						fe_dg,block_structure[2],
							 						fe_dg,block_structure[3],
							 						fe_dg,block_structure[4]));
				break;
			}

			// case 6:
			// {
			// 			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
			// 									   FE_Nothing<dim>(),block_structure[1],
			// 				 						FE_Nothing<dim>(),block_structure[2],
			// 				 						FE_Nothing<dim>(),block_structure[3],
			// 				 						FE_Nothing<dim>(),block_structure[4],
			// 										FE_Nothing<dim>(),block_structure[5]));

			// 			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
			// 									   fe_dg,block_structure[1],
			// 				 						FE_Nothing<dim>(),block_structure[2],
			// 				 						FE_Nothing<dim>(),block_structure[3],
			// 				 						FE_Nothing<dim>(),block_structure[4],
			// 										FE_Nothing<dim>(),block_structure[5]));


			// 			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
			// 									   fe_dg,block_structure[1],
			// 				 						fe_dg,block_structure[2],
			// 				 						FE_Nothing<dim>(),block_structure[3],
			// 				 						FE_Nothing<dim>(),block_structure[4],
			// 										FE_Nothing<dim>(),block_structure[5]));


			// 			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
			// 									   fe_dg,block_structure[1],
			// 				 						fe_dg,block_structure[2],
			// 				 						fe_dg,block_structure[3],
			// 				 						FE_Nothing<dim>(),block_structure[4],
			// 										FE_Nothing<dim>(),block_structure[5]));


			// 			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
			// 									   fe_dg,block_structure[1],
			// 				 						fe_dg,block_structure[2],
			// 				 						fe_dg,block_structure[3],
			// 				 						fe_dg,block_structure[4],
			// 										FE_Nothing<dim>(),block_structure[5]));


			// 			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
			// 									   fe_dg,block_structure[1],
			// 				 						fe_dg,block_structure[2],
			// 				 						fe_dg,block_structure[3],
			// 				 						fe_dg,block_structure[4],
			// 										fe_dg,block_structure[5]));
			// 	break;
			// }


			// case 7:
			// {
			// 			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
			// 									   FE_Nothing<dim>(),block_structure[1],
			// 				 						FE_Nothing<dim>(),block_structure[2],
			// 				 						FE_Nothing<dim>(),block_structure[3],
			// 				 						FE_Nothing<dim>(),block_structure[4],
			// 										FE_Nothing<dim>(),block_structure[5],
			// 										FE_Nothing<dim>(),block_structure[6]));

			// 			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
			// 									   fe_dg,block_structure[1],
			// 				 						FE_Nothing<dim>(),block_structure[2],
			// 				 						FE_Nothing<dim>(),block_structure[3],
			// 				 						FE_Nothing<dim>(),block_structure[4],
			// 										FE_Nothing<dim>(),block_structure[5],
			// 										FE_Nothing<dim>(),block_structure[6]));



			// 			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
			// 									   fe_dg,block_structure[1],
			// 				 						fe_dg,block_structure[2],
			// 				 						FE_Nothing<dim>(),block_structure[3],
			// 				 						FE_Nothing<dim>(),block_structure[4],
			// 										FE_Nothing<dim>(),block_structure[5],
			// 										FE_Nothing<dim>(),block_structure[6]));

			// 			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
			// 									   fe_dg,block_structure[1],
			// 				 						fe_dg,block_structure[2],
			// 				 						fe_dg,block_structure[3],
			// 				 						FE_Nothing<dim>(),block_structure[4],
			// 										FE_Nothing<dim>(),block_structure[5],
			// 										FE_Nothing<dim>(),block_structure[6]));

			// 			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
			// 									   fe_dg,block_structure[1],
			// 				 						fe_dg,block_structure[2],
			// 				 						fe_dg,block_structure[3],
			// 				 						fe_dg,block_structure[4],
			// 										FE_Nothing<dim>(),block_structure[5],
			// 										FE_Nothing<dim>(),block_structure[6]));


			// 			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
			// 									   fe_dg,block_structure[1],
			// 				 						fe_dg,block_structure[2],
			// 				 						fe_dg,block_structure[3],
			// 				 						fe_dg,block_structure[4],
			// 										fe_dg,block_structure[5],
			// 										FE_Nothing<dim>(),block_structure[6]));


			// 			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
			// 									   fe_dg,block_structure[1],
			// 				 						fe_dg,block_structure[2],
			// 				 						fe_dg,block_structure[3],
			// 				 						fe_dg,block_structure[4],
			// 										fe_dg,block_structure[5],
			// 										fe_dg,block_structure[6]));


			// 	break;
			// }

			// case 8:
			// {

			// 	finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
			// 									   FE_Nothing<dim>(),block_structure[1],
			// 				 						FE_Nothing<dim>(),block_structure[2],
			// 				 						FE_Nothing<dim>(),block_structure[3],
			// 				 						FE_Nothing<dim>(),block_structure[4],
			// 										FE_Nothing<dim>(),block_structure[5],
			// 										FE_Nothing<dim>(),block_structure[6],
			// 										FE_Nothing<dim>(),block_structure[7]));
				
			// 	finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
			// 									   fe_dg,block_structure[1],
			// 				 						FE_Nothing<dim>(),block_structure[2],
			// 				 						FE_Nothing<dim>(),block_structure[3],
			// 				 						FE_Nothing<dim>(),block_structure[4],
			// 										FE_Nothing<dim>(),block_structure[5],
			// 										FE_Nothing<dim>(),block_structure[6],
			// 										FE_Nothing<dim>(),block_structure[7]));


			// 	finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
			// 									   fe_dg,block_structure[1],
			// 				 						fe_dg,block_structure[2],
			// 				 						FE_Nothing<dim>(),block_structure[3],
			// 				 						FE_Nothing<dim>(),block_structure[4],
			// 										FE_Nothing<dim>(),block_structure[5],
			// 										FE_Nothing<dim>(),block_structure[6],
			// 										FE_Nothing<dim>(),block_structure[7]));


			// 	finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
			// 									   fe_dg,block_structure[1],
			// 				 						fe_dg,block_structure[2],
			// 				 						fe_dg,block_structure[3],
			// 				 						FE_Nothing<dim>(),block_structure[4],
			// 										FE_Nothing<dim>(),block_structure[5],
			// 										FE_Nothing<dim>(),block_structure[6],
			// 										FE_Nothing<dim>(),block_structure[7]));

			// 	finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
			// 									   fe_dg,block_structure[1],
			// 				 						fe_dg,block_structure[2],
			// 				 						fe_dg,block_structure[3],
			// 				 						fe_dg,block_structure[4],
			// 										FE_Nothing<dim>(),block_structure[5],
			// 										FE_Nothing<dim>(),block_structure[6],
			// 										FE_Nothing<dim>(),block_structure[7]));


			// 	finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
			// 									   fe_dg,block_structure[1],
			// 				 						fe_dg,block_structure[2],
			// 				 						fe_dg,block_structure[3],
			// 				 						fe_dg,block_structure[4],
			// 										fe_dg,block_structure[5],
			// 										FE_Nothing<dim>(),block_structure[6],
			// 										FE_Nothing<dim>(),block_structure[7]));



			// 	finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
			// 									   fe_dg,block_structure[1],
			// 				 						fe_dg,block_structure[2],
			// 				 						fe_dg,block_structure[3],
			// 				 						fe_dg,block_structure[4],
			// 										fe_dg,block_structure[5],
			// 										fe_dg,block_structure[6],
			// 										FE_Nothing<dim>(),block_structure[7]));

			// 	finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
			// 										   fe_dg,block_structure[1],
			// 				 						   fe_dg,block_structure[2],
			// 				 						   fe_dg,block_structure[3],
			// 				 						fe_dg,block_structure[4],
			// 										fe_dg,block_structure[5],
			// 										fe_dg,block_structure[6],
			// 										fe_dg,block_structure[7]));
			// 	break;
			// }
			default:
			{
				AssertThrow(1 == 0, ExcMessage("Should not have reached here"));
				break;
			}
		}

			
}


template<int dim>
void 
Base_Solver<dim>::allocate_fe_index(const unsigned int present_cycle,const unsigned int total_cycles)
{
	typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
	unsigned int counter = 0;
	Vector<double> refined(max_fe_index+1);
	refined = 0;

	// the very first cycle of m refinement
	 if (present_cycle == 0)
	 {
		for (; cell != endc ; cell++)
			cell->set_active_fe_index(0);

		cell = dof_handler.begin_active();

	}


	else
	{

		for (; cell != endc ; cell++)
		{
			const double error_value = VelocitySpace_error_per_cell(counter);

			for (unsigned int i = 0 ; i < VelocitySpace_error_tolerance.size()-1 ; i++)
				if (error_value >= VelocitySpace_error_tolerance[i] && error_value <= VelocitySpace_error_tolerance[i+1])
				{
					refined(i)++;
					cell->set_active_fe_index(i);
				}
 
			counter++;

		}		
	}


	for (unsigned long int i = 0 ; i < refined.size() ; i++)
		std::cout << "Fe Index: " << i << " Times: " << refined(i) << std::endl;



}
