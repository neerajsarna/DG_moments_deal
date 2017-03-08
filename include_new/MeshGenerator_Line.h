	template<int dim>
	void
	Base_MeshGenerator<dim>::mesh_internal_line()
	{
		GridGenerator::hyper_cube(triangulation,-0.5,0.5);

		switch(constants.problem_type)
		{
			case inflow_outflow:
			{
				set_bid_line_inflow_outflow();
				break;
			}
			
			case heat_conduction:
			{
				set_bid_line_heat_conduction();
				break;
			}

			case periodic:
			case lid_driven_cavity:
			{
				AssertThrow(1 == 0, ExcMessage("Should not have reached here"));
				break;
			}

			default:
			{
				AssertThrow(1 == 0 ,ExcMessage("Should not have reached here"));
				break;
			}
		}

		triangulation.refine_global(3);
	}

	template<int dim>
	void 
	Base_MeshGenerator<dim>::set_bid_line_inflow_outflow()
	{
		typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();

		for (; cell != endc ; cell++)
			for (unsigned int face = 0 ; face < GeometryInfo<dim>::faces_per_cell ; face++)
				{
					// inflow of the line
					if (fabs(cell->face(face)->center()(0)) < 1e-10)
						cell->face(face)->set_boundary_id(101);

					// outflow of the line
					if (fabs(cell->face(face)->center()(0) - 1) < 1e-10)
						cell->face(face)->set_boundary_id(102);
				}

	}

	template<int dim>
	void 
	Base_MeshGenerator<dim>::set_bid_line_heat_conduction()
	{
		typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();

		for (; cell != endc ; cell++)
			for (unsigned int face = 0 ; face < GeometryInfo<dim>::faces_per_cell ; face++)
				{
					// left edge of the line
					if (fabs(cell->face(face)->center()(0)) < 1e-10)
						cell->face(face)->set_boundary_id(0);

					// right edge of the line
					if (fabs(cell->face(face)->center()(0) - 1) < 1e-10)
						cell->face(face)->set_boundary_id(1);
				}		
	}