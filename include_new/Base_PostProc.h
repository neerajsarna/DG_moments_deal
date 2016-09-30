namespace PostProc
{
	using namespace dealii;
	struct stat st = { .st_dev = 0 };

	template<int dim>
	class
	Base_PostProc
	{
		public:
			Base_PostProc(const constant_data &constants,
						  ExactSolution::Base_ExactSolution<dim> *exact_solution,
						  const DoFHandler<dim> *dof_handle,
						  const MappingQ<dim> *mapping_obj);

			// we need the following information from the calling routine for this class to work
			ExactSolution::Base_ExactSolution<dim> *base_exactsolution;
			const constant_data constants;
			const DoFHandler<dim> *dof_handler;
			const MappingQ<dim> *mapping;


			struct output_files
      		{
        			std::string file_for_convergence_tables;
        			std::string file_for_grid;
        			std::string file_for_num_solution;
        			std::string file_for_exact_solution;
        			std::string file_for_error;
      		};

      		output_files output_file_names;

      		bool used_midpoint;
      		bool used_qgauss;

      		// we prescribe the filename for the output routines
      		void prescribe_file_names();

      		// we make the required directories
      		void make_directories();

      		// error evaluation based upon gauss quadrature
			void error_evaluation_QGauss(const Vector<double> &solution,
										const unsigned int active_cells,
										const double hMax,
									   ConvergenceTable &convergence_table);

			// error evaluation based upon midpoint rule
			void error_evaluation_QMidpoint(const Vector<double> &solution,
											const unsigned int active_cells,
											const double hMax,
									   		ConvergenceTable &convergence_table);

			// print convergence table to a file
			void print_convergence_table(ConvergenceTable &convergence_table);

			// print the solution and the error to a file
			// With the help of the two booleans we decide whether we want to print the 
			// solution or we want to print the error.
			void print_solution(
									  const Triangulation<dim> &triangulation,
									  const Vector<double> &solution);

			void print_error(const Triangulation<dim> &triangulation,
							 const Vector<double> &solution);

			void print_exactsolution(const Triangulation<dim> &triangulation);

	};

	template<int dim>
	Base_PostProc<dim>::Base_PostProc(const constant_data &constants,
									  ExactSolution::Base_ExactSolution<dim> *exact_solution,
									  const DoFHandler<dim> *dof_handle,
									  const MappingQ<dim> *mapping_obj)
	:
	constants(constants),
	base_exactsolution(base_exactsolution),
	dof_handler(dof_handle),
	mapping(mapping_obj)
	{
		// we prescribe the file names where we wish to write
		prescribe_file_names();

		used_midpoint = false;
		used_qgauss = false;

		make_directories();
	}

	template<int dim>
	void
	Base_PostProc<dim>::prescribe_file_names()
	{
		  const unsigned int poly_degree = dof_handler->get_fe()->degree();

		  output_file_names.file_for_convergence_tables = constants.sub_directory_names[2] + "/convergence_table_global_degree_"
                                                          + std::to_string(poly_degree);

          output_file_names.file_for_num_solution = constants.sub_directory_names[1] + "/numerical_solution_global_degree_"
                                                          + std::to_string(poly_degree)+"_DOF_"+std::to_string(dof_handler->n_dofs());


          output_file_names.file_for_exact_solution = constants.sub_directory_names[1] + "/exact_solution_global_degree_"
                                                          + std::to_string(poly_degree)+"_DOF_"+std::to_string(dof_handler->n_dofs());


          output_file_names.file_for_error = constants.sub_directory_names[1] + "/error_global_degree_"
                                                          + std::to_string(poly_degree)+"_DOF_"+std::to_string(dof_handler->n_dofs());

	}

	template<int dim>
	void 
	Base_PostProc<dim>::make_directories()
	{
		const unsigned int num_outputs = 5;

		if (stat(constants.main_output_dir.c_str(),&st) == -1)
			mkdir(constants.main_output_dir.c_str(),0777);

		for (unsigned int i = 0 ; i < num_outputs ; i ++)
			if (stat(constants.sub_directory_names[i].c_str(),&st) == -1)
				mkdir(constants.sub_directory_names[i].c_str(),0777);

	}

	// evaluate the error using gaussian qaudrature
	template<int dim>
	void
	Base_PostProc<dim>::error_evaluation_QGauss(const Vector<double> &solution,
											   const unsigned int active_cells,
											   const double hMax,
											   ConvergenceTable &convergence_table)
	{
		used_qgauss = true;

		// error variable comes from Basics. The following is the component in which
		// we wish to find the error.
        unsigned int component = constants.variable_map.find(constants.error_variable)->second;

        const unsigned int ngp = constants.p + 1;
        // error per cell of the domain
        Vector<double> error_per_cell(active_cells);      

        ComponentSelectFunction<dim> weight(component,constants.nEqn);                              // used to compute only the error in theta

        // computation of L2 error
        VectorTools::integrate_difference (*mapping,*dof_handler,solution,
          									base_exactsolution,
          									error_per_cell,
          									QGauss<dim>(ngp),
          									VectorTools::L2_norm,
          									&weight);  


        const double L2_error = error_per_cell.l2_norm();

        // computation of L_inifinity error
        error_per_cell = 0;
        VectorTools::integrate_difference (*mapping,*dof_handler,solution,
          									base_exactsolution,
          									error_per_cell,
          									QGauss<dim>(ngp),
          									VectorTools::Linfty_norm,
          									&weight);  
        

        const double Linfty_error = error_per_cell.linfty_norm();
                
        std::string column_name_L2;
        std::string column_name_Linfty;

        column_name_L2 = "#L2 in u(Using QGauss)" + std::to_string(component);
        column_name_Linfty = "#Linfty in u(Using QGauss)" + std::to_string(component);

        std::string error_details;
        error_details = " L2_error: " + std::to_string(L2_error) +" Linfty_error: " + std::to_string(Linfty_error) + 
        				" #DOF: " + std::to_string(dof_handler->n_dofs()) +
                        " #Cells " + std::to_string(active_cells); 

        std::cout << "<<<<<<<<<<<<<<<<<<<<<Error detail>>>>>>>>>>>>>>>>>>>>>"<< std::endl ;
        std::cout << error_details << std::endl;
        std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl ;

        convergence_table.add_value(column_name_L2,L2_error);
        convergence_table.add_value(column_name_Linfty,Linfty_error);
        convergence_table.add_value("#degree of freedom",dof_handler->n_dofs());
        convergence_table.add_value("#number of cells",active_cells);
        convergence_table.add_value("#hMax",hMax);

        convergence_table.set_scientific(column_name_L2,true);
        convergence_table.set_scientific(column_name_Linfty,true);
        convergence_table.set_scientific("#hMax",true);
	}


	// evaluate the error from the solution using midpoint quadrature rule
	template<int dim>
	void
	Base_PostProc<dim>::error_evaluation_QMidpoint(const Vector<double> &solution,
											   		const unsigned int active_cells,
											   		const double hMax,
											   		ConvergenceTable &convergence_table)
	{

		used_midpoint = false;

		// error variable comes from Basics. The following is the component in which
		// we wish to find the error.
        unsigned int component = constants.variable_map.find(constants.error_variable)->second;

        // error per cell of the domain
        Vector<double> error_per_cell(active_cells);      

        ComponentSelectFunction<dim> weight(component,constants.nEqn);                              // used to compute only the error in theta

        // computation of L2 error
        VectorTools::integrate_difference (*mapping,*dof_handler,solution,
          									base_exactsolution,
          									error_per_cell,
          									QMidpoint<dim>(),
          									VectorTools::L2_norm,
          									&weight);  


        const double L2_error = error_per_cell.l2_norm();

        // computation of L_inifinity error
        error_per_cell = 0;
        VectorTools::integrate_difference (*mapping,*dof_handler,solution,
          									base_exactsolution,
          									error_per_cell,
          									QMidpoint<dim>(),
          									VectorTools::Linfty_norm,
          									&weight);  
        

        const double Linfty_error = error_per_cell.linfty_norm();
                
        std::string column_name_L2;
        std::string column_name_Linfty;

        column_name_L2 = "#L2 in u(Using QMidpoint)" + std::to_string(component);
        column_name_Linfty = "#Linfty in u(Using QMidpoint)" + std::to_string(component);

        std::string error_details;
        error_details = " L2_error: " + std::to_string(L2_error) +" Linfty_error: " + std::to_string(Linfty_error) + " #DOF: " + std::to_string(dof_handler->n_dofs()) +
                        " #Cells " + std::to_string(active_cells); 

        std::cout << "<<<<<<<<<<<<<<<<<<<<<Error detail>>>>>>>>>>>>>>>>>>>>>"<< std::endl ;
        std::cout << error_details << std::endl;
        std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl ;

        convergence_table.add_value(column_name_L2,L2_error);
        convergence_table.add_value(column_name_Linfty,Linfty_error);
        convergence_table.add_value("#degree of freedom",dof_handler->n_dofs());
        convergence_table.add_value("#number of cells",active_cells);
        convergence_table.add_value("#hMax",hMax);

        convergence_table.set_scientific(column_name_L2,true);
        convergence_table.set_scientific(column_name_Linfty,true);
        convergence_table.set_scientific("#hMax",true);
	}

	// print the convergence table to a file
	template<int dim> 
	void 
	Base_PostProc<dim>::print_convergence_table(ConvergenceTable &convergence_table)
        {
           unsigned int component = constants.variable_map.find(constants.error_variable)->second;
            std::string column_name_L2;
        	std::string column_name_Linfty;

        	if (used_midpoint)
        	{
        		column_name_L2 = "#L2 in u(Using QMidpoint)" + std::to_string(component);
        		column_name_Linfty = "#Linfty in u(Using QMidpoint)" + std::to_string(component);        		
        	}

        	if (used_qgauss)
        	{
        		column_name_L2 = "#L2 in u(Using QGauss)" + std::to_string(component);
        		column_name_Linfty = "#Linfty in u(Using QGauss)" + std::to_string(component);        		
        	}
		   convergence_table.evaluate_convergence_rates(column_name_L2,"#number of cells", ConvergenceTable::reduction_rate_log2);
           convergence_table.evaluate_convergence_rates(column_name_Linfty,"#number of cells", ConvergenceTable::reduction_rate_log2);

           std::ofstream output_convergence(output_file_names.file_for_convergence_tables);
           convergence_table.write_text(output_convergence);

        }

     // compute the solution at all the vertices and print it to a file
    template<int dim>
   	void 
   	Base_PostProc<dim>::
   	print_solution(
   					    const Triangulation<dim> &triangulation,
   					    const Vector<double> &solution)
   	{
   	typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();

	FILE *fp_solution;

	fp_solution = fopen(output_file_names.file_for_num_solution.c_str(),"w+");

	AssertThrow(fp_solution != NULL,ExcMessage("file not open"));

	fprintf(fp_solution, "#%s\n","x y at the midpoint of each cell all the solution components");

	for (; cell != endc ; cell++)
		for (unsigned int vertex = 0 ; vertex < GeometryInfo<dim>::vertices_per_cell ; vertex ++)
		{

		Vector<double> solution_value(constants.nEqn);
		VectorTools::point_value(&dof_handler, solution, cell->vertex(vertex),solution_value);	
		

		for (unsigned int space = 0 ; space < dim ; space ++)
			fprintf(fp_solution, "%f ",cell->vertex(vertex)[space]);
		

		for (unsigned int i = 0 ; i < this->nEqn ; i++)
			fprintf(fp_solution, "%f ",solution_value(i));
		

		fprintf(fp_solution, "\n");
			
		}

	fclose(fp_solution);
   	}

   	// we compute the error at all the vertices and then print it to a file
    template<int dim>
   	void 
   	Base_PostProc<dim>::
   	print_error(
   					    const Triangulation<dim> &triangulation,
   					    const Vector<double> &solution)
   	{
   	typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();

	FILE *fp_error;

	fp_error = fopen(output_file_names.file_for_error.c_str(),"w+");

	AssertThrow(fp_error != NULL,ExcMessage("file not open"));

	fprintf(fp_error, "#%s\n","x y at the midpoint of each cell and the error in every component");

	for (; cell != endc ; cell++)
		for (unsigned int vertex = 0 ; vertex < GeometryInfo<dim>::vertices_per_cell ; vertex ++)
		{

		Vector<double> solution_value(constants.nEqn);
		Vector<double> exact_solution_value(constants.nEqn);
		Vector<double> error_value(constants.nEqn);

		VectorTools::point_value(*dof_handler, solution, cell->vertex(vertex),solution_value);	
		base_exactsolution->vector_value(cell->vertex(vertex),exact_solution_value);

		for (unsigned int i = 0 ; i < constants.nEqn ; i++)
			error_value(i) = fabs(solution_value(i)-exact_solution_value(i));

		for (unsigned int space = 0 ; space < dim ; space ++)
			fprintf(fp_error, "%f ",cell->vertex(vertex)[space]);

		for (unsigned int i = 0 ; i < this->nEqn ; i++)
			fprintf(fp_error, "%f ",error_value(i));

		fprintf(fp_error, "\n");
			
		}


	fclose(fp_error);
   	}


    template<int dim>
   	void 
   	Base_PostProc<dim>::
   	print_exactsolution(const Triangulation<dim> &triangulation)
   	{
   	typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();

	FILE *fp_exact;

	fp_exact = fopen(output_file_names.file_for_exact_solution.c_str(),"w+");

	AssertThrow(fp_exact != NULL,ExcMessage("file not open"));

	fprintf(fp_exact, "#%s\n","cell coordinates and the values of the exact solution ");

	for (; cell != endc ; cell++)
		for (unsigned int vertex = 0 ; vertex < GeometryInfo<dim>::vertices_per_cell ; vertex ++)
		{

		Vector<double> exact_solution_value(constants.nEqn);

		base_exactsolution->vector_value(cell->vertex(vertex),exact_solution_value);

		for (unsigned int space = 0 ; space < dim ; space ++)
			fprintf(fp_exact, "%f ",cell->vertex(vertex)[space]);

		for (unsigned int i = 0 ; i < this->nEqn ; i++)
			fprintf(fp_exact, "%f ",exact_solution_value(i));

		fprintf(fp_exact, "\n");
			
		}

	fclose(fp_exact);
   	}
}