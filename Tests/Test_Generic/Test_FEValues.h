namespace Test_FEValues
{
	using namespace dealii;

	// we would like to know whether all the cells have the same fe_values and the fe_shapevalues
	TEST(FEValues,HandlesFEValues)
	{
                const unsigned int dim = 2;
		const unsigned int p = 2;
                ASSERT_EQ(dim,2) << "3D not implemented" << std::endl;

		SphericalManifold<dim> boundary;
		Triangulation<dim> triangulation;
	        
		Point<dim> center;
                const double inner_radius = 0.5;
                const double outer_radius = 2.0;
                const unsigned int parts = 10;

                center(0) = 0.0;
                center(1) = 0.0;

                GridGenerator::hyper_shell (triangulation,
                                                                center, inner_radius, outer_radius,
                                                        parts);

       		triangulation.set_all_manifold_ids_on_boundary(0);
        	triangulation.set_manifold(0,boundary);

		
		FE_DGQ<dim> finite_element(p);
		DoFHandler<dim> dof_handler(triangulation);
		
		dof_handler.distribute_dofs(finite_element);
		MappingQ<dim> mapping(p+1);

		UpdateFlags update_flags = update_quadrature_points |
                	           update_values            |
                        	   update_gradients;
	

		QGauss<dim> quadrature(p+1);
	
		FEValues<dim> fe_values(mapping,finite_element,quadrature,update_flags);

		typename Triangulation<dim>::active_cell_iterator cell=triangulation.begin_active(), endc = triangulation.end();


		
	}
}
