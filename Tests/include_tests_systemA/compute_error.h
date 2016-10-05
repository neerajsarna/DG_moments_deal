// the following routine return the value of the error corresponding to a particualr refine cycle
	double compute_error(Vector<double> &solution, 
						const MappingQ<dim> &mapping,
						ExactSolution::Base_ExactSolution<dim> *exact_solution,
						const unsigned int error_component,
						const unsigned int poly_order,
						const unsigned int active_cells,
						const unsigned int nEqn,
						DoFHandler<dim> *dof_handler)
	{	
        unsigned int component = error_component;
        
        constants.variable_map.find(constants.error_variable)->second;

        const unsigned int ngp = poly_order + 1;

        // error per cell of the domain
        Vector<double> error_per_cell(active_cells);      

        ComponentSelectFunction<dim> weight(component,nEqn);                              // used to compute only the error in theta

        // computation of L2 error
        VectorTools::integrate_difference (*mapping,*dof_handler,solution,
          									*base_exactsolution,
          									error_per_cell,
          									QGauss<dim>(ngp),
          									VectorTools::L2_norm,
          									&weight);  


        const double L2_error = error_per_cell.l2_norm();

        return(L2_error);
	}