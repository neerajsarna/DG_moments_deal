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

      		double compute_residual(Vector<double> &solution,FESystem<dim> &finite_element,
										  EquationGenerator::Base_EquationGenerator<dim> *system_info,
										  const unsigned int active_cells);

      		// error evaluation based upon gauss quadrature
			void error_evaluation_QGauss(const Vector<double> &solution,
										const unsigned int active_cells,
										double &error_value,
										const double hMax,
									   ConvergenceTable &convergence_table,
									   const double residual);

			double compute_L2_norm(const Vector<double> &solution,const unsigned int active_cells);

			// return the value of the L2_erro
			double L2_error_QGauss(const Vector<double> &solution,const unsigned int active_cells);
			double Linfty_error_QGauss(const Vector<double> &solution,const unsigned int active_cells);


			// error evaluation based upon midpoint rule
			void error_evaluation_QMidpoint(const Vector<double> &solution,
											const unsigned int active_cells,
											double &error_value,
											const double hMax,
									   		ConvergenceTable &convergence_table,
							const double residual);

			double L2_error_QMidpoint(const Vector<double> &solution,const unsigned int active_cells);
			double Linfty_error_QMidpoint(const Vector<double> &solution,const unsigned int active_cells);

			// print convergence table to a file
			void print_convergence_table_to_file(ConvergenceTable &convergence_table);


			void print_convergence_table_to_file(const double L2_error,const double Linfty_error,
												const double hMax,const unsigned int active_cells);

			// print the solution and the error to a file
			// With the help of the two booleans we decide whether we want to print the 
			// solution or we want to print the error.
			void print_solution_to_file(const Triangulation<dim> &triangulation,
										const Vector<double> &solution);

			void print_error_to_file(const Triangulation<dim> &triangulation,
							 		 const Vector<double> &solution);

			void print_exactsolution_to_file(const Triangulation<dim> &triangulation);

			// print the solution depending upon printing options
			void print_options(const Triangulation<dim> &triangulation,
							   const Vector<double> &solution,
							   const unsigned int present_cycle,
							   const unsigned int refine_cycle,
							   ConvergenceTable &convergence_table);

			// same as above but prints the error to the file manually
			void print_options(const Triangulation<dim> &triangulation,
							   const Vector<double> &solution,
							   const unsigned int present_cycle,
							   const unsigned int refine_cycle,
							   const double L2_error,
							   const double Linfty_error,
							   const unsigned int active_cells,
							   const double hMax);

			// write the values of computational constants to a file
			void create_stamp();

			// we compute the error from our code and print it manually to file
			// is useful when using external grid
			void compute_error_print_manuel(const Vector<double> &solution,
											const unsigned int active_cells,
											const double hMax);

			// compute error manually. Sometimes we need to compute the error of the 
			// orginal system. The following routine is for that purpose
			void compute_error_unsymmetric(const Triangulation<dim> &triangulation,
										   const Vector<double> &solution,
										   const Sparse_matrix &S_half_inv,
										   double &error_value);

	};

	template<int dim>
	Base_PostProc<dim>::Base_PostProc(const constant_data &constants,
									  ExactSolution::Base_ExactSolution<dim> *exact_solution,
									  const DoFHandler<dim> *dof_handle,
									  const MappingQ<dim,dim> *mapping_obj)
	:
	base_exactsolution(exact_solution),
	constants(constants),
	dof_handler(dof_handle),
	mapping(mapping_obj)
	{
		prescribe_file_names();

		used_midpoint = false;
		used_qgauss = false;

		make_directories();

		// we leave the details of the computational parameters in the main directory
		create_stamp();

	}

	template<int dim>
	void
	Base_PostProc<dim>::prescribe_file_names()
	{
		  const unsigned int poly_degree = constants.p;

		  Assert(constants.sub_directory_names.size() != 0,ExcMessage("Not initialized"));
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

	template<int dim>
	void
	Base_PostProc<dim>::create_stamp()
	{
		FILE *fp;
		std::string filename = constants.main_output_dir + "/computational_parameters.txt";

		fp = fopen(filename.c_str(),"w+");
		AssertThrow(fp != NULL , ExcMessage("Could not open file for stamping"));

		std::string parameters = "A0 " + std::to_string(constants.A0) + "\n" +
								 "A1 " + std::to_string(constants.A1) + "\n" +
								 "A2 " + std::to_string(constants.A2) + "\n" +
								 "poly degree " + std::to_string(constants.p) + "\n" +
								 "mapping order" + std::to_string(constants.mapping_order) + "\n" +
								 "uW " + std::to_string(constants.uW) + "\n" +
								 "tau " + std::to_string(constants.tau)+ "\n" +
								 "nEqn " + std::to_string(constants.nEqn) + "\n" +
								 "nBC" + std::to_string(constants.nBC) + "\n" +
								 "bc_type" + std::to_string(constants.bc_type)+ "\n" +
								 "force_type " + std::to_string(constants.force_type);

		fprintf(fp, "%s\n",parameters.c_str());
		fclose(fp);
	}

	template<int dim>
	double 
	Base_PostProc<dim>::compute_residual(Vector<double> &solution,FESystem<dim> &finite_element,
										  EquationGenerator::Base_EquationGenerator<dim> *system_info,
										  const unsigned int active_cells)
	{
			typename DoFHandler<dim>::active_cell_iterator cell = dof_handler->begin_active(), endc = dof_handler->end();

			// we dont need much precision with the residual computation so we will simply use the midpoint quadrature rule 
			// for the following computation

			QGauss<dim> quadrature(constants.p + 2);
			const int total_ngp = quadrature.size();

			// an array for quadrature points
			std::vector<Point<dim>> q_points(quadrature.size());


			// we now have the update flags
			const UpdateFlags update_flags               =  update_gradients
                                                     | update_q_points
                                                     | update_JxW_values;

            // the fe_v object that will be used to provide us with the jacobian values and the shape values
			FEValues<dim>  fe_v(*mapping,finite_element,quadrature, update_flags);

			// the results returned from the point_gradient function
			std::vector<Tensor<1,dim,double>> value_gradient(constants.nEqn);
			std::vector<Vector<double>> value_per_quad(dim,Vector<double>(constants.nEqn));

			// l2 norm of the residual per cell
			Vector<double> residual_per_cell(active_cells);

			// the residual vector per cell
			Vector<double> difference_per_quad(constants.nEqn);

			// value of the source terms at the quadrature points
			std::vector<Vector<double>> source_term_value(total_ngp,Vector<double>(constants.nEqn));

			// value of the jacobians at all the quadrature points
			std::vector<double> Jacobian(total_ngp);
			unsigned int counter = 0;

			// vector containing the value of the solution
			Vector<double> value(constants.nEqn);
			Vector<double> exact_solution(constants.nEqn);

			residual_per_cell = 0;
			MatrixOpt::Base_MatrixOpt matrix_opt;

			for (; cell != endc ; cell++)
			{
				fe_v.reinit(cell);
				q_points = fe_v.get_quadrature_points();
				system_info->source_term(fe_v.get_quadrature_points(),source_term_value);
				Jacobian = fe_v.get_JxW_values();


				// we now compute the l2 norm of the solution per cell
				for (int q = 0 ;q < total_ngp ; q++)
				{
					difference_per_quad = 0;

					// the gradient value has to be multiplied by the system matrices
					VectorTools::point_gradient	(*dof_handler,
												 solution,
												q_points[q],
												value_gradient);	

					// the point value has to be multiplied by the Production term
					VectorTools::point_value(*dof_handler,
											solution,
											q_points[q],
											value);

					//base_exactsolution->vector_value(q_points[q],exact_solution);

					// now we take the transpose of values so that we can use the sparse vector product already developed
					
		
					for (unsigned int eq = 0 ; eq < constants.nEqn ; eq++)
						for (unsigned int space = 0 ; space < dim ; space ++)
							value_per_quad[space](eq) = value_gradient[eq][space];

					//the residual from the convective term
					for (unsigned int space = 0 ; space < dim ; space++)
						difference_per_quad += matrix_opt.Sparse_matrix_dot_Vector(system_info->system_data.A[space].matrix,value_per_quad[space]);

					// the residual from the right hand side of the equation
					difference_per_quad -= matrix_opt.Sparse_matrix_dot_Vector(system_info->system_data.P.matrix,value);
										   
					// computation of the l2 norm of the solution in this cell
					for (unsigned int eq = 0 ; eq < 1 ; eq++)
					{
						// contribution from the source term
						difference_per_quad(eq) -= source_term_value[q](eq);

						// l2-norm of the solutoin
						residual_per_cell(counter) += pow(difference_per_quad(eq),2) * Jacobian[q];
					}						

				}

				residual_per_cell(counter) = sqrt(residual_per_cell(counter));
				counter ++;
			}

			return(residual_per_cell.l2_norm());


	}


	template<int dim>
	double
	Base_PostProc<dim>::compute_L2_norm(const Vector<double> &solution,const unsigned int active_cells)
	{
		unsigned int component = constants.variable_map.find(constants.error_variable)->second;
		Vector<double> norm_per_cell(active_cells); 
		const unsigned int ngp = constants.p + 1;

		ComponentSelectFunction<dim> weight(component,constants.nEqn);                              // used to compute only the error in theta

		ZeroFunction<dim> zero_function(constants.nEqn);

        // computation of L2 error
        VectorTools::integrate_difference (*mapping,*dof_handler,solution,
          									zero_function,
          									norm_per_cell,
          									QGauss<dim>(ngp),
          									VectorTools::L2_norm,
          									&weight); 


        return(norm_per_cell.l2_norm());

	}

	// evaluate the error using gaussian qaudrature
	template<int dim>
	void
	Base_PostProc<dim>::error_evaluation_QGauss(const Vector<double> &solution,
											   const unsigned int active_cells,
											   double &error_value,
											   const double hMax,
											   ConvergenceTable &convergence_table,
						    const double residual)
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
          									*base_exactsolution,
          									error_per_cell,
          									QGauss<dim>(ngp),
          									VectorTools::L2_norm,
          									&weight);  


        const double L2_error = error_per_cell.l2_norm();

        error_value = L2_error;
        // computation of L_inifinity error
        error_per_cell = 0;
        VectorTools::integrate_difference (*mapping,*dof_handler,solution,
          									*base_exactsolution,
          									error_per_cell,
          									QGauss<dim>(ngp),
          									VectorTools::Linfty_norm,
          									&weight);  
        

        const double Linfty_error = error_per_cell.linfty_norm();
                
        std::string column_name_L2;
        std::string column_name_Linfty;
	std::string column_name_residual = "#Residual" ;

        column_name_L2 = "#L2 in u(Using QGauss)" + std::to_string(component);
        column_name_Linfty = "#Linfty in u(Using QGauss)" + std::to_string(component);

        std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Error Details>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"  << std::endl;
       	printf("L2_error: %e, Linf_error: %e, #DOF: %llu, #Cells %llu, #Residual %e \n",L2_error,Linfty_error,dof_handler->n_dofs(),active_cells,residual);
       	std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Error Details>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

        convergence_table.add_value(column_name_L2,L2_error);
        convergence_table.add_value(column_name_Linfty,Linfty_error);
        convergence_table.add_value("#degree of freedom",dof_handler->n_dofs());
        convergence_table.add_value("#number of cells",active_cells);
        convergence_table.add_value("#hMax",hMax);
	convergence_table.add_value("#Residual",residual);

        convergence_table.set_scientific(column_name_L2,true);
        convergence_table.set_scientific(column_name_Linfty,true);
        convergence_table.set_scientific("#hMax",true);
	convergence_table.set_scientific(column_name_residual,true);

	}


	// returns the L2 error using the gauss quadrature
	template<int dim>
	double
	Base_PostProc<dim>::L2_error_QGauss(const Vector<double> &solution,
										const unsigned int active_cells)
	{
		// error variable comes from Basics. The following is the component in which
		// we wish to find the error.
        unsigned int component = constants.variable_map.find(constants.error_variable)->second;

        const unsigned int ngp = constants.p + 1;
        // error per cell of the domain
        Vector<double> error_per_cell(active_cells);      

        ComponentSelectFunction<dim> weight(component,constants.nEqn);                              // used to compute only the error in theta

        // computation of L2 error
        VectorTools::integrate_difference (*mapping,*dof_handler,solution,
          									*base_exactsolution,
          									error_per_cell,
          									QGauss<dim>(ngp),
          									VectorTools::L2_norm,
          									&weight);  


        return (error_per_cell.l2_norm());
                
	}

	template<int dim>
	double
	Base_PostProc<dim>::Linfty_error_QGauss(const Vector<double> &solution,
											const unsigned int active_cells)
	{
		// error variable comes from Basics. The following is the component in which
		// we wish to find the error.
        unsigned int component = constants.variable_map.find(constants.error_variable)->second;

        const unsigned int ngp = constants.p + 1;
        // error per cell of the domain
        Vector<double> error_per_cell(active_cells);      

        ComponentSelectFunction<dim> weight(component,constants.nEqn);                              // used to compute only the error in theta

        // computation of L_inifinity error
        error_per_cell = 0;
        VectorTools::integrate_difference (*mapping,*dof_handler,solution,
          									*base_exactsolution,
          									error_per_cell,
          									QGauss<dim>(ngp),
          									VectorTools::Linfty_norm,
          									&weight);  
        

        return(error_per_cell.linfty_norm());
                
	}

	// evaluate the error from the solution using midpoint quadrature rule
	template<int dim>
	void
	Base_PostProc<dim>::error_evaluation_QMidpoint(const Vector<double> &solution,
											   		const unsigned int active_cells,
											   		double &error_value,
											   		const double hMax,
											   		ConvergenceTable &convergence_table,
						     const double residual)
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
          									*base_exactsolution,
          									error_per_cell,
          									QMidpoint<dim>(),
          									VectorTools::L2_norm,
          									&weight);  


        const double L2_error = error_per_cell.l2_norm();

        error_value = L2_error;

        // computation of L_inifinity error
        error_per_cell = 0;
        VectorTools::integrate_difference (*mapping,*dof_handler,solution,
          									*base_exactsolution,
          									error_per_cell,
          									QMidpoint<dim>(),
          									VectorTools::Linfty_norm,
          									&weight);  
        

        const double Linfty_error = error_per_cell.linfty_norm();
                
        std::string column_name_L2;
        std::string column_name_Linfty;

        column_name_L2 = "#L2 in u(Using QMidpoint)" + std::to_string(component);
        column_name_Linfty = "#Linfty in u(Using QMidpoint)" + std::to_string(component);

       	
       	std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Error Details>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"  <<  std::endl;
       	printf("L2_error: %e, Linf_error: %e, #DOF: %u, #Cells %u\n",L2_error,Linfty_error,dof_handler->n_dofs(),active_cells);
       	std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Error Details>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" <<  std::endl;


        convergence_table.add_value(column_name_L2,L2_error);
        convergence_table.add_value(column_name_Linfty,Linfty_error);
        convergence_table.add_value("#degree of freedom",dof_handler->n_dofs());
        convergence_table.add_value("#number of cells",active_cells);
        convergence_table.add_value("#hMax",hMax);

        convergence_table.set_scientific(column_name_L2,true);
        convergence_table.set_scientific(column_name_Linfty,true);
        convergence_table.set_scientific("#hMax",true);
	}


	// returns the L2 error using the gauss quadrature
	template<int dim>
	double
	Base_PostProc<dim>::L2_error_QMidpoint(const Vector<double> &solution,
										  const unsigned int active_cells)
	{
		// error variable comes from Basics. The following is the component in which
		// we wish to find the error.
        unsigned int component = constants.variable_map.find(constants.error_variable)->second;

        const unsigned int ngp = constants.p + 1;
        // error per cell of the domain
        Vector<double> error_per_cell(active_cells);      

        ComponentSelectFunction<dim> weight(component,constants.nEqn);                              // used to compute only the error in theta

        // computation of L2 error
        VectorTools::integrate_difference (*mapping,*dof_handler,solution,
          									*base_exactsolution,
          									error_per_cell,
          									QGauss<dim>(ngp),
          									VectorTools::L2_norm,
          									&weight);  


        return (error_per_cell.l2_norm());

                
	}

	template<int dim>
	double
	Base_PostProc<dim>::Linfty_error_QMidpoint(const Vector<double> &solution,
												const unsigned int active_cells)
	{
		// error variable comes from Basics. The following is the component in which
		// we wish to find the error.
        unsigned int component = constants.variable_map.find(constants.error_variable)->second;

        const unsigned int ngp = constants.p + 1;
        // error per cell of the domain
        Vector<double> error_per_cell(active_cells);      

        ComponentSelectFunction<dim> weight(component,constants.nEqn);                              // used to compute only the error in theta

        // computation of L_inifinity error
        error_per_cell = 0;
        VectorTools::integrate_difference (*mapping,*dof_handler,solution,
          									*base_exactsolution,
          									error_per_cell,
          									QGauss<dim>(ngp),
          									VectorTools::Linfty_norm,
          									&weight);  
        

        return(error_per_cell.linfty_norm());
                
	}


	// print the convergence table to a file
	template<int dim> 
	void 
	Base_PostProc<dim>::print_convergence_table_to_file(ConvergenceTable &convergence_table)
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

	// print convergence table manually
	template<int dim> 
	void 
	Base_PostProc<dim>::print_convergence_table_to_file(const double L2_error,const double Linfty_error,
														const double hMax,const unsigned int active_cells)
        {
           unsigned int component = constants.variable_map.find(constants.error_variable)->second;
           
           FILE *fp;
           fp = fopen(output_file_names.file_for_convergence_tables.c_str(),"a+");

           fprintf(fp, "L2_error Linf_error hMax active_cells (u%u)\n",component);

           fprintf(fp, "%f %f %f %u\n",L2_error,Linfty_error,hMax,active_cells);

           fclose(fp);

        }

     // compute the solution at all the vertices and print it to a file
    template<int dim>
   	void 
   	Base_PostProc<dim>::
   	print_solution_to_file(
   					    const Triangulation<dim> &triangulation,
   					    const Vector<double> &solution)
   	{
   	typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();

	FILE *fp_solution;

	fp_solution = fopen(output_file_names.file_for_num_solution.c_str(),"w+");

	AssertThrow(fp_solution != NULL,ExcMessage("file not open"));

	fprintf(fp_solution, "#%s\n","x y at the midpoint of each cell all the solution components");
	Vector<double> solution_value(constants.nEqn);
	

	for (; cell != endc ; cell++)
	{
		solution_value = 0;
		VectorTools::point_value(*dof_handler, solution, cell->center(),solution_value);	

                for (unsigned int space = 0 ; space < dim ; space ++)
                        fprintf(fp_solution, "%f ",cell->center()[space]);

		// we only print variables uptill theta
                for (int i = 0 ; i < 4 ; i++)
                        fprintf(fp_solution, "%f ",solution_value(i));

		
		fprintf(fp_solution,"\n");
	}


	fclose(fp_solution);
   	}

   	// we compute the error at all the vertices and then print it to a file
    template<int dim>
   	void 
   	Base_PostProc<dim>::
   	print_error_to_file(
   					    const Triangulation<dim> &triangulation,
   					    const Vector<double> &solution)
   	{
   	typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();

	FILE *fp_error;

	fp_error = fopen(output_file_names.file_for_error.c_str(),"w+");

	AssertThrow(fp_error != NULL,ExcMessage("file not open"));

	fprintf(fp_error, "#%s\n","x y at the midpoint of each cell and the error in every component");

	// we can print  
	for (; cell != endc ; cell++)
		{

		Vector<double> solution_value(constants.nEqn);
		Vector<double> exact_solution_value(constants.nEqn);
		Vector<double> error_value(constants.nEqn);

		VectorTools::point_value(*dof_handler, solution, cell->center(),solution_value);	
		base_exactsolution->vector_value(cell->center(),exact_solution_value);

		for ( int i = 0 ; i < constants.nEqn ; i++)
			error_value(i) = fabs(solution_value(i)-exact_solution_value(i));

		for (int space = 0 ; space < dim ; space ++)
			fprintf(fp_error, "%f ",cell->center()[space]);

		for (int i = 0 ; i < constants.nEqn ; i++)
			fprintf(fp_error, "%f ",error_value(i));

		fprintf(fp_error, "\n");
			
		}


	fclose(fp_error);
   	}


    template<int dim>
   	void 
   	Base_PostProc<dim>::
   	print_exactsolution_to_file(const Triangulation<dim> &triangulation)
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

		for (int i = 0 ; i < constants.nEqn ; i++)
			fprintf(fp_exact, "%f ",exact_solution_value(i));

		fprintf(fp_exact, "\n");
			
		}

	fclose(fp_exact);
   	}

   	template<int dim>
   	void 
   	Base_PostProc<dim>::
   	print_options(const Triangulation<dim> &triangulation,
   				  const Vector<double> &solution,
   				  const unsigned int present_cycle,
   				  const unsigned int total_cycles,
   				  ConvergenceTable &convergence_table)
   	{

   		if (constants.print_all)
   		{
   			if (constants.print_solution)
   				print_solution_to_file(triangulation,solution);

   			if(constants.print_error)
   				print_error_to_file(triangulation,solution);

   			if(constants.print_exactsolution)
   				print_exactsolution_to_file(triangulation);

   		}

   		else
   		{
   			// only print in the final cycle
   			if (present_cycle == total_cycles - 1)
   			{
   				if(constants.print_solution)
   					print_solution_to_file(triangulation,solution);

   				if(constants.print_error)
   					print_error_to_file(triangulation,solution);

   				if(constants.print_exactsolution)
   					print_exactsolution_to_file(triangulation);
   			}
   		}

   		if (constants.print_convergence_table)
 			print_convergence_table_to_file(convergence_table);

   	}

   	template<int dim>
   	void 
   	Base_PostProc<dim>::
   	print_options(const Triangulation<dim> &triangulation,
   				  const Vector<double> &solution,
   				  const unsigned int present_cycle,
   				  const unsigned int total_cycles,
   				  const double L2_error,
   				  const double Linfty_error,
   				  const unsigned int active_cells,
   				  const double hMax)
   	{

   		// if we wish to print for all the refinement cycles
   		if (constants.print_all)
   		{
   			if (constants.print_solution)
   				print_solution_to_file(triangulation,solution);

   			if(constants.print_error)
   				print_error_to_file(triangulation,solution);

   			if(constants.print_exactsolution)
   				print_exactsolution_to_file(triangulation);

   		}

   		// if we wish to print only in the end of the computation
   		else
   		{
   			// only print in the final cycle
   			if (present_cycle == total_cycles - 1)
   			{
   				if(constants.print_solution)
   					print_solution_to_file(triangulation,solution);

   				if(constants.print_error)
   					print_error_to_file(triangulation,solution);

   				if(constants.print_exactsolution)
   					print_exactsolution_to_file(triangulation);
   			}
   		}

   		if (constants.print_convergence_table)
 			print_convergence_table_to_file(L2_error,Linfty_error,hMax,active_cells);

   	}

   	template<int dim>
   	void
   	Base_PostProc<dim>::compute_error_print_manuel(const Vector<double> &solution,
											   		const unsigned int active_cells,
											   		const double hMax)
   	{
 		// error variable comes from Basics. The following is the component in which
		// we wish to find the error.
        unsigned int component = constants.variable_map.find(constants.error_variable)->second;

        // error per cell of the domain
        Vector<double> error_per_cell(active_cells);      

        ComponentSelectFunction<dim> weight(component,constants.nEqn);                              // used to compute only the error in theta

        // computation of L2 error
        VectorTools::integrate_difference (*mapping,*dof_handler,solution,
          									*base_exactsolution,
          									error_per_cell,
          									QMidpoint<dim>(),
          									VectorTools::L2_norm,
          									&weight);  


        const double L2_error = error_per_cell.l2_norm();

        // computation of L_inifinity error
        error_per_cell = 0;
        VectorTools::integrate_difference (*mapping,*dof_handler,solution,
          									*base_exactsolution,
          									error_per_cell,
          									QMidpoint<dim>(),
          									VectorTools::Linfty_norm,
          									&weight);  
        

        const double Linfty_error = error_per_cell.linfty_norm(); 
        
        FILE *fp;

        // append the existing file and write values to it
        fp = fopen(output_file_names.file_for_convergence_tables.c_str(),"a+"); 		
        AssertThrow(fp != NULL,ExcMessage("Cant open file for convergence table writting"));

        fprintf(fp, "L2_error Linfity_error active_cells hMax\n");
        fprintf(fp, "%f %f %u %f \n",L2_error, Linfty_error, active_cells,hMax);

        fclose(fp);

   	}

   	template<int dim>
   	void
   	Base_PostProc<dim>::compute_error_unsymmetric(const Triangulation<dim> &triangulation,
												   const Vector<double> &solution,
												   const Sparse_matrix &S_half_inv,
												   double &error)
   	{
   		typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(),
   														  endc = triangulation.end();

   		Assert(triangulation.n_active_cells() != 0, ExcMessage("No grid available"));
   		Assert(solution.size() !=0 ,ExcMessage("solution not available"));

   		Vector<double> exact_solution_value(constants.nEqn);
   		Vector<double> solution_value(constants.nEqn);
   		Vector<double> error_value(triangulation.n_active_cells());
   		const unsigned int component = constants.variable_map.find(constants.error_variable)->second;

   		MatrixOpt::Base_MatrixOpt matrix_opt;


   		for (; cell != endc ; cell++)
   		{
 				
					VectorTools::point_value(*dof_handler, solution, cell->center(),solution_value);// compute the exact solution first
					base_exactsolution->vector_value(cell->center(),exact_solution_value);  	

					exact_solution_value = 	matrix_opt.Sparse_matrix_dot_Vector(S_half_inv,exact_solution_value);	
					solution_value = 	matrix_opt.Sparse_matrix_dot_Vector(S_half_inv,solution_value);	

					// compute the Linf error at the cell center
					error_value(cell->index()) = fabs(exact_solution_value(component) - solution_value(component));
   		}

   		const double Linf_error = error_value.linfty_norm();
   		error = Linf_error;
   		std::cout << "Linf error from unsymmetric system: " << Linf_error << " #Cells: "<< triangulation.n_active_cells() << std::endl;
   	}
}
