namespace Test_FEValues
{
	using namespace dealii;

	// we would like to know whether all the cells have the same fe_values and the fe_shapevalues
	TEST(FEValues,HandlesFEValues)
	{
                const unsigned int dim = 2;
                ASSERT_EQ(dim,2) << "3D not implemented" << std::endl;

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
	}
}
