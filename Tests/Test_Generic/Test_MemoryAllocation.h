using namespace dealii;

TEST(Memory_Allocation,HandlesMemoryAllocation)
{
		const unsigned int dim = 2;

		Triangulation<dim> triangulation;
		GridGenerator::subdivided_hyper_cube(triangulation,50);

		FESystem<dim> finite_element(FE_DGQ<dim>(1),100);
        DoFHandler<dim> dof_handler(triangulation);


		TrilinosWrappers::SparseMatrix global_matrix;

		dof_handler.distribute_dofs(finite_element);

        printf("sparsity pattern\n"); 
        fflush(stdout);

		DynamicSparsityPattern dsp(dof_handler.n_dofs(),dof_handler.n_dofs());

	//std::cout << "making flux sparsity pattern " << std::endl;
		DoFTools::make_flux_sparsity_pattern (dof_handler, dsp);

        printf("global_matrix\n"); 
        fflush(stdout);

		global_matrix.reinit(dsp);  

}