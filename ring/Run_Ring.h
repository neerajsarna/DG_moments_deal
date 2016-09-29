template<int dim>
void
Base_Solver<dim>::distribute_dof_allocate_matrix_ring()
{
    dof_handler.distribute_dofs(finite_element);

    DynamicSparsityPattern dsp(dof_handler.n_dofs(),dof_handler.n_dofs());


    DoFTools::make_flux_sparsity_pattern (dof_handler, dsp);

    sparsity_pattern.copy_from(dsp);
 
    global_matrix.reinit(sparsity_pattern);   
 
    solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());

}

template<int dim>
void
Base_Solver<dim>::run_ring(const unsigned int refine_cycles)
{
	for (unsigned int i = 0 ; i < refine_cycles ; i ++)
	{
		if (i == 0)
		{
			// first we generate the grid
			switch(constants.mesh_options)
			{
				case generate_internal:
				{
			// generate internally
					this->mesh_internal_ring();
					break;
				}
				case read_msh:
				{
				// read the gmsh generated ring
					this->mesh_gmsh_ring();
					break;
				}
				default:
				{
					Assert(1 == 0,ExcMessage("Should not have reached here"));
					break;
				}
			}
		}

		//now we distribute the dofs and allocate the memory for the global matrix
		distribute_dof_allocate_matrix_ring();

		// the following routine assembles
		switch(constants.assembly_type)
		{
			case meshworker:
			{
				assemble_system_meshworker();
				break;
			}

			// case manuel:
			// {
			// 	assemble_system();
			// 	break;
			// }
			default :
			{
				Assert(1 == 0,ExcMessage("Should not have reached here"));
				break;
			}
		}

		//LinearSolve(global_matrix,system_rhs,solution,Pardiso);

	}

}