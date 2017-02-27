	// generate a ring shaped mesh internally
	template<int dim>
	void
	Base_MeshGenerator<dim>::mesh_internal_ring()
	{
		
		Point<dim> center;
		const double inner_radius = constants.inner_radius;	// default values 0.5
		const double outer_radius = constants.outer_radius;	// default values 2.0
		const unsigned int parts = 10;
		const unsigned int refinement_level = constants.initial_refinement; // default value is 1

		center(0) = 0.0;
		center(1) = 0.0;

		GridGenerator::hyper_shell (triangulation,
        							center, inner_radius, outer_radius,
              						parts);


		switch(constants.problem_type)
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
			{
				break;
			}

			case periodic:
			{
				AssertThrow(1 == 0,ExcNotImplemented());
				break;
			}

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


