	// generate a ring shaped mesh internally
	template<int dim>
	void
	Base_MeshGenerator<dim>::mesh_internal_ring()
	{
		
		Point<dim> center;
		const double inner_radius = 0.5;	// default values 0.5
		const double outer_radius = 2.0;	// default values 2.0
		const unsigned int parts = 10;
		const unsigned int refinement_level = initial_refinement; // default value is 1

		center(0) = 0.0;
		center(1) = 0.0;

		std::cout << "Inner radius " << inner_radius << " Outer radius " << outer_radius << std::endl;
		
		GridGenerator::hyper_shell (triangulation,
        							center, inner_radius, outer_radius,
              						parts);


		switch(problem_type)
		{
			// heat conduction between the outer and the inner ring
			case heat_conduction:
			{

				triangulation.set_all_manifold_ids_on_boundary(0);
				triangulation.set_manifold(0,boundary);
				triangulation.refine_global(refinement_level);
				break;
			}

			// in case of inflow and outflow the boundary ids will be prescribed by gmsh
			case inflow_outflow:
			case periodic:
			case lid_driven_cavity:
			{
				AssertThrow(1 == 0,ExcNotImplemented());
				break;
			}

			default:
			{
				AssertThrow( 1 == 0, ExcMessage("Should not have reached here"));
				break;
			}
		}


	}


