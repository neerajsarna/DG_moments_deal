template<int num_flux,int dim>
void
Solver_DG<num_flux,dim>
::output_solution_details(const Triangulation<dim> &triangulation,
							const string file_solution,
						const string file_exact,
						 const string file_error)const 
{
	typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();

	QMidpoint<dim> quadrature;
	const UpdateFlags update_flags  = update_values | update_JxW_values | update_quadrature_points;
	FEValues<dim>  fe_v(mapping,finite_element,quadrature, update_flags);

	FILE *fp_solution, *fp_exact, *fp_error;

	fp_solution = fopen(file_solution.c_str(),"w+");
	fp_exact = fopen(file_exact.c_str(),"w+");
	fp_error = fopen(file_error.c_str(), "w+");

	AssertThrow(fp_solution != NULL,ExcMessage("file not open"));
	AssertThrow(fp_exact != NULL,ExcMessage("file not open"));
	AssertThrow(fp_error != NULL,ExcMessage("file not open"));

	fprintf(fp_solution, "#%s\n","x y at the midpoint of each cell all the solution components");
	fprintf(fp_exact, "#%s\n","x y at the midpoint of each cell all the exact solution components");
	fprintf(fp_error, "#%s\n","x y at the midpoint of each cell all the error components");

	for (; cell != endc ; cell++)
	{
		fe_v.reinit(cell);

		for (unsigned int vertex = 0 ; vertex < GeometryInfo<dim>::vertices_per_cell ; vertex ++)
		{

		Vector<double> solution_value(this->nEqn);
		Vector<double> exact_solution_value(this->nEqn);
		Vector<double> error_value(this->nEqn);

		VectorTools::point_value(dof_handler, solution, cell->vertex(vertex),solution_value);	
		exact_solution->vector_value(cell->vertex(vertex),exact_solution_value);

		for (unsigned int i = 0 ; i < this->nEqn ; i++)
			error_value(i) = fabs(solution_value(i)-exact_solution_value(i));

		
		switch(equation_system_data->system_type)
		{
			case symmetric:
			{
				// we would like to convert back to our original variables
				Sparse_matrix_dot_Vector(equation_system_data->system_data[solve_system].S_half_inv,solution_value);
				break;
			}

			// don't do anything in case of using an un_symmetric system
			case un_symmetric:
				break;
			
		}

		

		fprintf(fp_solution, "%f %f ",cell->vertex(vertex)[0],cell->vertex(vertex)[1]);
		fprintf(fp_exact, "%f %f ",cell->vertex(vertex)[0],cell->vertex(vertex)[1]);
		fprintf(fp_error, "%f %f ",cell->vertex(vertex)[0],cell->vertex(vertex)[1]);

		for (unsigned int i = 0 ; i < this->nEqn ; i++)
		{
			fprintf(fp_solution, "%f ",solution_value(i));
			fprintf(fp_exact, "%f ",exact_solution_value(i));
			fprintf(fp_error, "%f ",error_value(i));
		}

		fprintf(fp_solution, "\n");
		fprintf(fp_exact, "\n");
		fprintf(fp_error, "\n");
			
		}
	}

	fclose(fp_solution);
	fclose(fp_exact);
	fclose(fp_error);
}


// same as the above function, but does not compute the error values and simply outputs the solution
template<int num_flux,int dim>
void
Solver_DG<num_flux,dim>
::output_solution_details(const Triangulation<dim> &triangulation,
							const string file_solution)const 
{
	typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();

	QMidpoint<dim> quadrature;
	const UpdateFlags update_flags  = update_values | update_JxW_values | update_quadrature_points;
	FEValues<dim>  fe_v(mapping,finite_element,quadrature, update_flags);

	FILE *fp_solution;

	fp_solution = fopen(file_solution.c_str(),"w+");
	AssertThrow(fp_solution != NULL,ExcMessage("file not open"));
	fprintf(fp_solution, "#%s\n","x y at the midpoint of each cell all the solution components");
	
	for (; cell != endc ; cell++)
	{
		fe_v.reinit(cell);

		for (unsigned int vertex = 0 ; vertex < GeometryInfo<dim>::vertices_per_cell ; vertex ++)
		{

		Vector<double> solution_value(nEqn);


		VectorTools::point_value(dof_handler, solution, cell->vertex(vertex),solution_value);	
		
		switch(equation_system_data->system_type)
		{
			case symmetric:
			{
				// we would like to convert back to our original variables
				Sparse_matrix_dot_Vector(equation_system_data->system_data[solve_system].S_half_inv,solution_value);
				break;
			}

			// don't do anything in case of using an un_symmetric system
			case un_symmetric:
				break;
			
		}

		fprintf(fp_solution, "%f %f ",cell->vertex(vertex)[0],cell->vertex(vertex)[1]);
		

		for (unsigned int i = 0 ; i < this->nEqn ; i++)
			fprintf(fp_solution, "%f ",solution_value(i));
		

		fprintf(fp_solution, "\n");
		}
	}

	fclose(fp_solution);
	
}