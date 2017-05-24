namespace PostProc
{
	using namespace dealii;
	struct stat st = { .st_dev = 0 };

	template<int dim>
	class
	Base_PostProc
	{
	public:
		Base_PostProc(const constant_numerics &constants,
					  ExactSolution::Base_ExactSolution<dim> *exact_solution,
					  const std::vector<int> &nEqn,
				      const std::vector<int> &nBC);

			// we need the following information from the calling routine for this class to work
		ExactSolution::Base_ExactSolution<dim> *base_exactsolution;
		const constant_numerics constants;
		const std::vector<int> nEqn;
		const std::vector<int> nBC;

		MatrixOpt::Base_MatrixOpt matrix_opt;

		struct output_files
		{
			std::string file_for_convergence_tables;
			std::string file_for_grid;
			std::string file_for_num_solution;
			std::string file_for_exact_solution;
			std::string file_for_error;
			std::string file_for_velocity_space_error;
			std::string file_for_fe_index;
		};

		output_files output_file_names;

		bool used_midpoint;
		bool used_qgauss;
		bool class_initialized;

		void reinit(const hp::DoFHandler<dim> &dof_handler);
		void reinit(const DoFHandler<dim> &dof_handler);

      		// we prescribe the filename for the output routines
		void prescribe_file_names(const hp::DoFHandler<dim> &dof_handler);
		void prescribe_file_names(const DoFHandler<dim> &dof_handler);

      		// we make the required directories
		void make_directories();

		double compute_residual(Vector<double> &solution,const int n_active_cells,
											const MappingQ<dim> &mapping, const DoFHandler<dim> &dof_handler);


		double compute_residual(Vector<double> &solution,const int n_active_cells,
											const hp::MappingCollection<dim> &mapping, const hp::DoFHandler<dim> &dof_handler);


		double compute_energy(Vector<double> &solution,const int n_active_cells,
											const MappingQ<dim> &mapping, const DoFHandler<dim> &dof_handler);

		double compute_energy_bound();

      		// error evaluation based upon gauss quadrature
		void error_evaluation_QGauss(const Vector<double> &solution,
									const unsigned int active_cells,
									double &error_value,
									const double hMax,
									ConvergenceTable &convergence_table,
									const double residual,
									const hp::MappingCollection<dim> &mapping,
									const hp::DoFHandler<dim> &dof_handler);


      		// error evaluation based upon gauss quadrature
		Vector<double> return_error_per_cell(const Vector<double> &solution,
											 const unsigned int active_cells,
											 const hp::MappingCollection<dim> &mapping,
											 const hp::DoFHandler<dim> &dof_handler);


		void error_evaluation_QGauss(const Vector<double> &solution,
									const unsigned int active_cells,
									double &error_value,
									const double hMax,
									ConvergenceTable &convergence_table,
									const double residual,
									const MappingQ<dim> &mapping,
									const DoFHandler<dim> &dof_handler);




			// print convergence table to a file
		void print_convergence_table_to_file(ConvergenceTable &convergence_table);


		void print_convergence_table_to_file(const double L2_error,const double Linfty_error,
									const double hMax,const unsigned int active_cells);

			// print the solution and the error to a file
			// With the help of the two booleans we decide whether we want to print the 
			// solution or we want to print the error.
		void print_solution_to_file(const Triangulation<dim> &triangulation,
									const Vector<double> &solution,
									const Sparse_matrix &S_half_inv,
									const hp::DoFHandler<dim> &dof_handler);

		void print_solution_to_file(const Triangulation<dim> &triangulation,
									const Vector<double> &solution,
									const Sparse_matrix &S_half_inv,
									const DoFHandler<dim> &dof_handler);

		void print_solution_to_file(const Triangulation<dim> &triangulation,
									const Vector<double> &solution,
									const Sparse_matrix &S_half_inv,
									const DoFHandler<dim> &dof_handler,
									std::string &filename);



		void print_error_to_file(const Triangulation<dim> &triangulation,
								 const Vector<double> &solution,
								 const hp::DoFHandler<dim> &dof_handler);	

		void print_error_to_file(const Triangulation<dim> &triangulation,
								 const Vector<double> &solution,
								 const DoFHandler<dim> &dof_handler);	

		void print_exactsolution_to_file(const Triangulation<dim> &triangulation,const Sparse_matrix &S_half_inv);

		void print_VelocitySpace_error_to_file(const Triangulation<dim> &triangulation,const Vector<double> &error);

		// print the fe index of every cell
		void print_fe_index(const hp::DoFHandler<dim> &dof_handler);


			// print the solution depending upon printing options
		void print_options(const Triangulation<dim> &triangulation,
			const Vector<double> &solution,
			const unsigned int present_cycle,
			const unsigned int refine_cycle,
			ConvergenceTable &convergence_table,
			const Sparse_matrix &S_half_inv,
			const DoFHandler<dim> &dof_handler);

		void print_options(const Triangulation<dim> &triangulation,
						   const Vector<double> &solution,
						   const unsigned int present_cycle,
				           const unsigned int refine_cycle,
			  			   ConvergenceTable &convergence_table,
						   const Sparse_matrix &S_half_inv,
						   const hp::DoFHandler<dim> &dof_handler,
						   const Vector<double> &VelocitySpace_error_per_cell);

			// write the values of computational constants to a file
		void create_stamp(const hp::DoFHandler<dim> &dof_handler);
		void create_stamp(const DoFHandler<dim> &dof_handler);
		int max_equations;

		 Vector<double>
    	compute_lift_drag(const MappingQ<dim> &mapping,
    									  const FESystem<dim> &finite_element,
    									  const DoFHandler<dim> &dof_handler,
    									  const Sparse_matrix &S_half_inv,
    									  const Vector<double> &solution,
    									  ConvergenceTable &convergence_table,
    									  const unsigned int b_id_surface);

	};

	template<int dim>
	Base_PostProc<dim>::Base_PostProc(const constant_numerics &constants,
		ExactSolution::Base_ExactSolution<dim> *exact_solution,
		const std::vector<int> &nEqn,
		const std::vector<int> &nBC)
	:
	base_exactsolution(exact_solution),
	constants(constants),
	nEqn(nEqn),
	nBC(nBC)
	{
		// we take the last entry of the vector since it corresponds to the moment systems which has the 
		// maximum number of equations
			max_equations = nEqn[nEqn.size()-1];


	}

	template<int dim>
	void
	Base_PostProc<dim>::reinit(const hp::DoFHandler<dim> &dof_handler)
	{
		class_initialized = true;
		prescribe_file_names(dof_handler);
		

		used_midpoint = false;
		used_qgauss = false;

		make_directories();

		// we leave the details of the computational parameters in the main directory
		create_stamp(dof_handler);


	}

	template<int dim>
	void
	Base_PostProc<dim>::reinit(const DoFHandler<dim> &dof_handler)
	{
		class_initialized = true;
		prescribe_file_names(dof_handler);
		


		used_midpoint = false;
		used_qgauss = false;

		make_directories();

		// we leave the details of the computational parameters in the main directory
		create_stamp(dof_handler);

	}


	template<int dim>
	void
	Base_PostProc<dim>::prescribe_file_names(const hp::DoFHandler<dim> &dof_handler)
	{
		
		Assert(class_initialized,ExcMessage("Please initialize the post proc class"));
		const unsigned int poly_degree = constants.p;

		Assert(constants.sub_directory_names.size() != 0,ExcMessage("Not initialized"));
		output_file_names.file_for_convergence_tables = constants.sub_directory_names[2] + "/convergence_table_global_degree_"
		+ std::to_string(poly_degree);

		output_file_names.file_for_num_solution = constants.sub_directory_names[1] + "/numerical_solution_global_degree_"
		+ std::to_string(poly_degree)+"_DOF_"+std::to_string(dof_handler.n_dofs());

		output_file_names.file_for_exact_solution = constants.sub_directory_names[1] + "/exact_solution_global_degree_"
		+ std::to_string(poly_degree)+"_DOF_"+std::to_string(dof_handler.n_dofs());

		output_file_names.file_for_error = constants.sub_directory_names[1] + "/error_global_degree_"
		+ std::to_string(poly_degree)+"_DOF_"+std::to_string(dof_handler.n_dofs());

		output_file_names.file_for_velocity_space_error = constants.sub_directory_names[1] + "/error_velocity_space_global_degree_"
		+ std::to_string(poly_degree)+"_DOF_"+std::to_string(dof_handler.n_dofs());

		output_file_names.file_for_fe_index = constants.sub_directory_names[1] + "/fe_index_global_degree_"
		+ std::to_string(poly_degree)+"_DOF_"+std::to_string(dof_handler.n_dofs());
	}

	template<int dim>
	void
	Base_PostProc<dim>::prescribe_file_names(const DoFHandler<dim> &dof_handler)
	{
		Assert(class_initialized == true,ExcMessage("Please initialize the post proc class"));
		const unsigned int poly_degree = constants.p;

		Assert(constants.sub_directory_names.size() != 0,ExcMessage("Not initialized"));
		output_file_names.file_for_convergence_tables = constants.sub_directory_names[2] + "/convergence_table_global_degree_"
		+ std::to_string(poly_degree);

		output_file_names.file_for_num_solution = constants.sub_directory_names[1] + "/numerical_solution_global_degree_"
		+ std::to_string(poly_degree)+"_DOF_"+std::to_string(dof_handler.n_dofs());

		output_file_names.file_for_exact_solution = constants.sub_directory_names[1] + "/exact_solution_global_degree_"
		+ std::to_string(poly_degree)+"_DOF_"+std::to_string(dof_handler.n_dofs());

		output_file_names.file_for_error = constants.sub_directory_names[1] + "/error_global_degree_"
		+ std::to_string(poly_degree)+"_DOF_"+std::to_string(dof_handler.n_dofs());

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
		Base_PostProc<dim>::create_stamp(const hp::DoFHandler<dim> &dof_handler)
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
			"bc_type" + std::to_string(constants.bc_type)+ "\n" +
			"force_type " + std::to_string(constants.force_type) + 
			"#Dofs " + std::to_string(dof_handler.n_dofs()) + "\n"
			"#Cells " + std::to_string(dof_handler.get_triangulation().n_active_cells());

			for (unsigned long int eq = 0 ; eq < nEqn.size(); eq++)
			{
				parameters += "nEqn" + std::to_string(nEqn[eq]) + "\n";
				parameters += "nBC" + std::to_string(nEqn[eq]) + "\n";
			}


			fprintf(fp, "%s\n",parameters.c_str());
			fclose(fp);
		}

	template<int dim>
		void
		Base_PostProc<dim>::create_stamp(const DoFHandler<dim> &dof_handler)
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
			"bc_type" + std::to_string(constants.bc_type)+ "\n" +
			"force_type " + std::to_string(constants.force_type) + 
			"#Dofs " + std::to_string(dof_handler.n_dofs()) + "\n"
			"#Cells " + std::to_string(dof_handler.get_triangulation().n_active_cells());

			for (unsigned long int eq = 0 ; eq < nEqn.size(); eq++)
			{
				parameters += "nEqn" + std::to_string(nEqn[eq]) + "\n";
				parameters += "nBC" + std::to_string(nEqn[eq]) + "\n";
			}


			fprintf(fp, "%s\n",parameters.c_str());
			fclose(fp);
		}


	template<int dim>
		double 
		Base_PostProc<dim>::compute_residual(Vector<double> &solution,const int n_active_cells,
											const hp::MappingCollection<dim> &mapping, const hp::DoFHandler<dim> &dof_handler)
		{
			Assert(class_initialized == true,ExcMessage("Please initialize the post proc class"));

			const unsigned int ngp = constants.p + 1;
			Vector<double> error_per_cell(n_active_cells);
			QGauss<dim> quadrature_basic(ngp);			
			hp::QCollection<dim> hp_quadrature;

			for (unsigned long int i = 0 ; i < nEqn.size() ; i++)
				hp_quadrature.push_back(quadrature_basic);
				


			VectorTools::integrate_difference (mapping,dof_handler,solution,
				ZeroFunction<dim>(max_equations),
				error_per_cell,
				hp_quadrature,
				VectorTools::L2_norm); 


			return(error_per_cell.l2_norm()) ;


		}


	template<int dim>
		double 
		Base_PostProc<dim>::compute_residual(Vector<double> &solution,const int n_active_cells,
											const MappingQ<dim> &mapping, const DoFHandler<dim> &dof_handler)
		{
			Assert(class_initialized == true,ExcMessage("Please initialize the post proc class"));

			const unsigned int ngp = constants.p + 1;
			Vector<double> error_per_cell(n_active_cells);

			VectorTools::integrate_difference (mapping,dof_handler,solution,
				ZeroFunction<dim>(max_equations),
				error_per_cell,
				QGauss<dim>(ngp),
				VectorTools::L2_norm); 


			return(error_per_cell.l2_norm()) ;


		}

	// compute energy or the entropy function corresponding to the 1D moment system
	template<>
		double 
		Base_PostProc<1>::compute_energy(Vector<double> &solution,const int n_active_cells,
											const MappingQ<1> &mapping, const DoFHandler<1> &dof_handler)
		{
			Assert(class_initialized == true,ExcMessage("Please initialize the post proc class"));

			const unsigned int ngp = constants.p + 1;
			Vector<double> energy_per_cell(n_active_cells);

			VectorTools::integrate_difference (mapping,dof_handler,solution,
				ZeroFunction<1>(max_equations),
				energy_per_cell,
				QGauss<1>(ngp),
				VectorTools::L2_norm); 


			// divide by 0.5 due to gauss theorem, sqaure because of the entropy
			return(pow(energy_per_cell.l2_norm(),2) * 0.5);


		}

	// evaluate the error using gaussian qaudrature
	template<int dim>
	void
	Base_PostProc<dim>::error_evaluation_QGauss(const Vector<double> &solution,
		const unsigned int active_cells,
		double &error_value,
		const double hMax,
		ConvergenceTable &convergence_table,
		const double residual,
		const hp::MappingCollection<dim> &mapping,
		const hp::DoFHandler<dim> &dof_handler)
	{
		const unsigned int ngp = constants.p + 1;
		QGauss<dim> quadrature_basic(ngp);			
		hp::QCollection<dim> hp_quadrature;
		
		for (unsigned long int i = 0 ; i < nEqn.size() ; i++)
				hp_quadrature.push_back(quadrature_basic);


		used_qgauss = true;
		Assert(class_initialized == true,ExcMessage("Please initialize the post proc class"));

		// error variable comes from Basics. The following is the component in which
		// we wish to find the error.
        		// error variable comes from Basics. The following is the component in which
		// we wish to find the error.
		Assert(constants.variable_map.size() != 0,ExcMessage("Variable map not initialized"));
		Assert(constants.variable_map_1D.size() != 0,ExcMessage("Variable map not initialized"));

		unsigned int component;
		
		if (dim == 1)
			component = constants.variable_map_1D.find(constants.error_variable)->second;

		if (dim == 2)
			component = constants.variable_map.find(constants.error_variable)->second;

        // error per cell of the domain
		Vector<double> error_per_cell(active_cells);      

        ComponentSelectFunction<dim> weight(component,max_equations);                              // used to compute only the error in theta

        // computation of L2 error
        VectorTools::integrate_difference (mapping,dof_handler,solution,
        								  *base_exactsolution,
        								   error_per_cell,
        								   hp_quadrature,
        								   VectorTools::L2_norm,
        								   &weight);  


        const double L2_error = error_per_cell.l2_norm();

        error_value = L2_error;
        // computation of L_inifinity error
        error_per_cell = 0;
        VectorTools::integrate_difference (mapping,dof_handler,solution,
        									*base_exactsolution,
        									error_per_cell,
        									hp_quadrature,
        									VectorTools::Linfty_norm,
        									&weight);  
        

        const double Linfty_error = error_per_cell.linfty_norm();
        
        std::string column_name_L2;
        std::string column_name_Linfty;
        std::string column_name_residual = "#Residual" ;

        column_name_L2 = "#L2 in u(Using QGauss)" + std::to_string(component);
        column_name_Linfty = "#Linfty in u(Using QGauss)" + std::to_string(component);

        std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Error Details>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"  << std::endl;
        printf("L2_error: %e, Linf_error: %e, #DOF: %u, #Cells %u, #Residual %e \n",L2_error,Linfty_error,dof_handler.n_dofs(),active_cells,residual);
        std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Error Details>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

        convergence_table.add_value(column_name_L2,L2_error);
        convergence_table.add_value(column_name_Linfty,Linfty_error);
        convergence_table.add_value("#degree of freedom",dof_handler.n_dofs());
        convergence_table.add_value("#number of cells",active_cells);
        convergence_table.add_value("#hMax",hMax);
        convergence_table.add_value("#Residual",residual);

        convergence_table.set_scientific(column_name_L2,true);
        convergence_table.set_scientific(column_name_Linfty,true);
        convergence_table.set_scientific("#hMax",true);
        convergence_table.set_scientific(column_name_residual,true);

    }


	// evaluate the error using gaussian qaudrature
	template<int dim>
	void
	Base_PostProc<dim>::error_evaluation_QGauss(const Vector<double> &solution,
		const unsigned int active_cells,
		double &error_value,
		const double hMax,
		ConvergenceTable &convergence_table,
		const double residual,
		const MappingQ<dim> &mapping,
		const DoFHandler<dim> &dof_handler)
	{
		used_qgauss = true;
		Assert(class_initialized == true,ExcMessage("Please initialize the post proc class"));

		// error variable comes from Basics. The following is the component in which
		// we wish to find the error.
        		// error variable comes from Basics. The following is the component in which
		// we wish to find the error.
		Assert(constants.variable_map.size() != 0,ExcMessage("Variable map not initialized"));
		Assert(constants.variable_map_1D.size() != 0,ExcMessage("Variable map not initialized"));

		unsigned int component;
		
		if (dim == 1)
			component = constants.variable_map_1D.find(constants.error_variable)->second;

		if (dim == 2)
			component = constants.variable_map.find(constants.error_variable)->second;

		const unsigned int ngp = constants.p + 1;
        // error per cell of the domain
		Vector<double> error_per_cell(active_cells);      

        ComponentSelectFunction<dim> weight(component,max_equations);                              // used to compute only the error in theta

        // computation of L2 error
        VectorTools::integrate_difference (mapping,dof_handler,solution,
        	*base_exactsolution,
        	error_per_cell,
        	QGauss<dim>(ngp),
        	VectorTools::L2_norm,
        	&weight);  


        const double L2_error = error_per_cell.l2_norm();

        error_value = L2_error;
        // computation of L_inifinity error
        error_per_cell = 0;
        VectorTools::integrate_difference (mapping,dof_handler,solution,
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
        printf("L2_error: %e, Linf_error: %e, #DOF: %u, #Cells %u, #Residual %e  \n",
        		L2_error,Linfty_error,dof_handler.n_dofs(),active_cells,residual);
        std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Error Details>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

        convergence_table.add_value(column_name_L2,L2_error);
        convergence_table.add_value(column_name_Linfty,Linfty_error);
        convergence_table.add_value("#degree of freedom",dof_handler.n_dofs());
        convergence_table.add_value("#number of cells",active_cells);
        convergence_table.add_value("#hMax",hMax);
        convergence_table.add_value("#Residual",residual);

        convergence_table.set_scientific(column_name_L2,true);
        convergence_table.set_scientific(column_name_Linfty,true);
        convergence_table.set_scientific("#hMax",true);
        convergence_table.set_scientific(column_name_residual,true);

    }




	// print the convergence table to a file
	template<int dim> 
    void 
    Base_PostProc<dim>::print_convergence_table_to_file(ConvergenceTable &convergence_table)
    {
           		// error variable comes from Basics. The following is the component in which
		// we wish to find the error.
    	Assert(constants.variable_map.size() != 0,ExcMessage("Variable map not initialized"));
    	Assert(constants.variable_map_1D.size() != 0,ExcMessage("Variable map not initialized"));
		Assert(class_initialized == true,ExcMessage("Please initialize the post proc class"));

    	unsigned int component;
    	
    	if (dim == 1)
    		component = constants.variable_map_1D.find(constants.error_variable)->second;

    	if (dim == 2)
    		component = constants.variable_map.find(constants.error_variable)->second;
    	
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


    template<int dim>
    void 
    Base_PostProc<dim>::
    print_solution_to_file(
    					const Triangulation<dim> &triangulation,
    					const Vector<double> &solution,
    					const Sparse_matrix &S_half_inv,
    					const hp::DoFHandler<dim> &dof_handler)
    {
    	typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();

		Assert(class_initialized == true,ExcMessage("Please initialize the post proc class"));

    	FILE *fp_solution;

    	fp_solution = fopen(output_file_names.file_for_num_solution.c_str(),"w+");

    	AssertThrow(fp_solution != NULL,ExcMessage("file not open"));

    	fprintf(fp_solution, "#%s\n","x y at the midpoint of each cell all the solution components");
    	Vector<double> solution_value(max_equations);
    	    	
    	int variables_to_print;

    	if (dim == 2)
		variables_to_print = 9;	// we print all the moments till the heat flux

	if (dim == 1)
		variables_to_print = 5; // we print all the moments till the heat flux

	for (; cell != endc ; cell++)
	{
		for (unsigned int vertex = 0 ; vertex < GeometryInfo<dim>::vertices_per_cell ; vertex++)
		{
			solution_value = 0;
			VectorTools::point_value(dof_handler, solution, cell->vertex(vertex),solution_value);	

			// we now convert back to the conventional variables. That is the variables in unsymmetric system	
			solution_value = matrix_opt.Sparse_matrix_dot_Vector(S_half_inv,solution_value);

			for (unsigned int space = 0 ; space < dim ; space ++)
				fprintf(fp_solution, "%f ",cell->vertex(vertex)(space));

		// we only print variables uptill heat flux
			for (int i = 0 ; i < variables_to_print ; i++)
				fprintf(fp_solution, "%f ",solution_value(i));

			fprintf(fp_solution, "\n");
		}
	}


	fclose(fp_solution);
}

    template<int dim>
    void 
    Base_PostProc<dim>::
    print_solution_to_file(
    					const Triangulation<dim> &triangulation,
    					const Vector<double> &solution,
    					const Sparse_matrix &S_half_inv,
    					const DoFHandler<dim> &dof_handler)
    {
    	typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();
    	Assert(class_initialized == true,ExcMessage("Please initialize the post proc class"));

    	FILE *fp_solution;

    	fp_solution = fopen(output_file_names.file_for_num_solution.c_str(),"w+");

    	AssertThrow(fp_solution != NULL,ExcMessage("file not open"));

    	fprintf(fp_solution, "#%s\n","x y at the midpoint of each cell all the solution components");
    	Vector<double> solution_value(max_equations);
    	
    	int variables_to_print;

    	if (dim == 2)
		variables_to_print = 9;	// we print all the moments till the heat flux

	if (dim == 1)
		variables_to_print = 5; // we print all the moments till the heat flux

    	//variables_to_print = 6;

    	Assert(variables_to_print<=max_equations,ExcMessage("to many variables to print"));
	for (; cell != endc ; cell++)
	{
		for (unsigned int vertex = 0 ; vertex < GeometryInfo<dim>::vertices_per_cell ; vertex++)
		{
			solution_value = 0;
			VectorTools::point_value(dof_handler, solution, cell->vertex(vertex),solution_value);	

			// we now convert back to the conventional variables. That is the variables in unsymmetric system	
			solution_value = matrix_opt.Sparse_matrix_dot_Vector(S_half_inv,solution_value);

			for (unsigned int space = 0 ; space < dim ; space ++)
				fprintf(fp_solution, "%f ",cell->vertex(vertex)(space));

		// we only print variables uptill heat flux
			for (int i = 0 ; i < variables_to_print ; i++)
				fprintf(fp_solution, "%f ",solution_value(i));

			fprintf(fp_solution, "\n");
		}
	}


	fclose(fp_solution);
}


	// same as above but can take a user defined input filename
    template<int dim>
    void 
    Base_PostProc<dim>::
    print_solution_to_file(
    					const Triangulation<dim> &triangulation,
    					const Vector<double> &solution,
    					const Sparse_matrix &S_half_inv,
    					const DoFHandler<dim> &dof_handler,
    					std::string &filename)
    {
    	typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();
    	Assert(class_initialized == true,ExcMessage("Please initialize the post proc class"));

    	FILE *fp_solution;

    	fp_solution = fopen(filename.c_str(),"w+");

    	AssertThrow(fp_solution != NULL,ExcMessage("file not open"));

    	fprintf(fp_solution, "#%s\n","x y at the midpoint of each cell all the solution components");
    	Vector<double> solution_value(max_equations);
    	
    	int variables_to_print;

    	if (dim == 2)
		variables_to_print = 9;	// we print all the moments till the heat flux

	if (dim == 1)
		variables_to_print = 5; // we print all the moments till the heat flux

	for (; cell != endc ; cell++)
	{
		for (unsigned int vertex = 0 ; vertex < GeometryInfo<dim>::vertices_per_cell ; vertex++)
		{
			solution_value = 0;
			VectorTools::point_value(dof_handler, solution, cell->vertex(vertex),solution_value);	

			// we now convert back to the conventional variables. That is the variables in unsymmetric system	
			solution_value = matrix_opt.Sparse_matrix_dot_Vector(S_half_inv,solution_value);

			for (unsigned int space = 0 ; space < dim ; space ++)
				fprintf(fp_solution, "%f ",cell->vertex(vertex)(space));

		// we only print variables uptill heat flux
			for (int i = 0 ; i < variables_to_print ; i++)
				fprintf(fp_solution, "%f ",solution_value(i));

			fprintf(fp_solution, "\n");
		}
	}


	fclose(fp_solution);
}



   	// we compute the error at all the vertices and then print it to a file
    template<int dim>
void 
Base_PostProc<dim>::
print_error_to_file(
	const Triangulation<dim> &triangulation,
	const Vector<double> &solution,
	const hp::DoFHandler<dim> &dof_handler)
{
	Assert(class_initialized == true,ExcMessage("Please initialize the post proc class"));
	typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();

	FILE *fp_error;

	fp_error = fopen(output_file_names.file_for_error.c_str(),"w+");

	AssertThrow(fp_error != NULL,ExcMessage("file not open"));

	fprintf(fp_error, "#%s\n","x y at the midpoint of each cell and the error in every component");

	// we can print  
	for (; cell != endc ; cell++)
	{

		Vector<double> solution_value(max_equations);
		Vector<double> exact_solution_value(max_equations);
		Vector<double> error_value(max_equations);

		VectorTools::point_value(dof_handler, solution, cell->center(),solution_value);	
		base_exactsolution->vector_value(cell->center(),exact_solution_value);

		for ( int i = 0 ; i < max_equations ; i++)
			error_value(i) = fabs(solution_value(i)-exact_solution_value(i));

		for (int space = 0 ; space < dim ; space ++)
			fprintf(fp_error, "%f ",cell->center()[space]);

		for (int i = 0 ; i < max_equations ; i++)
			fprintf(fp_error, "%f ",error_value(i));

		fprintf(fp_error, "\n");
		
	}


	fclose(fp_error);
}

    template<int dim>
void 
Base_PostProc<dim>::
print_error_to_file(
	const Triangulation<dim> &triangulation,
	const Vector<double> &solution,
	const DoFHandler<dim> &dof_handler)
{
	Assert(class_initialized == true,ExcMessage("Please initialize the post proc class"));
	typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();

	FILE *fp_error;

	fp_error = fopen(output_file_names.file_for_error.c_str(),"w+");

	AssertThrow(fp_error != NULL,ExcMessage("file not open"));

	fprintf(fp_error, "#%s\n","x y at the midpoint of each cell and the error in every component");

	// we can print  
	for (; cell != endc ; cell++)
	{

		Vector<double> solution_value(max_equations);
		Vector<double> exact_solution_value(max_equations);
		Vector<double> error_value(max_equations);

		VectorTools::point_value(dof_handler, solution, cell->center(),solution_value);	
		base_exactsolution->vector_value(cell->center(),exact_solution_value);

		for ( int i = 0 ; i < max_equations ; i++)
			error_value(i) = fabs(solution_value(i)-exact_solution_value(i));

		for (int space = 0 ; space < dim ; space ++)
			fprintf(fp_error, "%f ",cell->center()[space]);

		for (int i = 0 ; i < max_equations ; i++)
			fprintf(fp_error, "%f ",error_value(i));

		fprintf(fp_error, "\n");
		
	}


	fclose(fp_error);
}

template<int dim>
void 
Base_PostProc<dim>::
print_exactsolution_to_file(const Triangulation<dim> &triangulation,const Sparse_matrix &S_half_inv)
{
	Assert(class_initialized == true,ExcMessage("Please initialize the post proc class"));
	typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();

	FILE *fp_exact;

	fp_exact = fopen(output_file_names.file_for_exact_solution.c_str(),"w+");

	AssertThrow(fp_exact != NULL,ExcMessage("file not open"));

	fprintf(fp_exact, "#%s\n","cell coordinates and the values of the exact solution ");

	for (; cell != endc ; cell++)
		for (unsigned int vertex = 0 ; vertex < GeometryInfo<dim>::vertices_per_cell ; vertex ++)
		{

			Vector<double> exact_solution_value(max_equations);

			base_exactsolution->vector_value(cell->vertex(vertex),exact_solution_value);

		// we now convert back to the original variables
			exact_solution_value = matrix_opt.Sparse_matrix_dot_Vector(S_half_inv,exact_solution_value);

			for (unsigned int space = 0 ; space < dim ; space ++)
				fprintf(fp_exact, "%f ",cell->vertex(vertex)[space]);

			for (int i = 0 ; i < max_equations ; i++)
				fprintf(fp_exact, "%f ",exact_solution_value(i));

			fprintf(fp_exact, "\n");
			
		}

		fclose(fp_exact);
}

template<int dim>
void 
Base_PostProc<dim>::
print_VelocitySpace_error_to_file(const Triangulation<dim> &triangulation,const Vector<double> &error)
{
	Assert(class_initialized == true,ExcMessage("Please initialize the post proc class"));
	Assert(error.size()!=0,ExcNotInitialized());
	typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();

	FILE *fp_exact;

	fp_exact = fopen(output_file_names.file_for_velocity_space_error.c_str(),"w+");

	AssertThrow(fp_exact != NULL,ExcMessage("file not open"));

	fprintf(fp_exact, "#%s\n","cell coordinates and the corresponding deviation from the equilibrium ");

	unsigned int counter = 0;

	for (; cell != endc ; cell++)
	{
		for (unsigned int space = 0 ; space < dim ; space ++)
				fprintf(fp_exact, "%f ",cell->center()[space]);

		fprintf(fp_exact, "%f ",error(counter));

		fprintf(fp_exact, "\n");

		counter++;

	}

		fclose(fp_exact);
}



template<int dim>
void 
Base_PostProc<dim>::
print_fe_index(const hp::DoFHandler<dim> &dof_handler)
{
	Assert(class_initialized == true,ExcMessage("Please initialize the post proc class"));
	typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
														 endc = dof_handler.end();

	FILE *fp_exact;

	fp_exact = fopen(output_file_names.file_for_fe_index.c_str(),"w+");

	AssertThrow(fp_exact != NULL,ExcMessage("file not open"));

	fprintf(fp_exact, "#%s\n","cell coordinates and the corresponding fe_index");

	unsigned int counter = 0;

	for (; cell != endc ; cell++)
	{
		for (unsigned int space = 0 ; space < dim ; space ++)
				fprintf(fp_exact, "%f ",cell->center()[space]);

	
		fprintf(fp_exact, "%d ",cell->active_fe_index());

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
		ConvergenceTable &convergence_table,
		const Sparse_matrix &S_half_inv,
		const DoFHandler<dim> &dof_handler)
	{
		Assert(class_initialized == true,ExcMessage("Please initialize the post proc class"));
		// if we would like to print for all the refinement cycles
		if (constants.print_all)
		{
			if (constants.print_solution)
				print_solution_to_file(triangulation,solution,S_half_inv,dof_handler);

			if(constants.print_error)
				print_error_to_file(triangulation,solution,dof_handler);

			if(constants.print_exactsolution)
				print_exactsolution_to_file(triangulation,S_half_inv);

		}

		else
		{
   			// only print in the final cycle
			if (present_cycle == total_cycles - 1)
			{
				if(constants.print_solution)
					print_solution_to_file(triangulation,solution,S_half_inv,dof_handler);

				if(constants.print_error)
					print_error_to_file(triangulation,solution,dof_handler);

				if(constants.print_exactsolution)
					print_exactsolution_to_file(triangulation,S_half_inv);

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
		ConvergenceTable &convergence_table,
		const Sparse_matrix &S_half_inv,
		const hp::DoFHandler<dim> &dof_handler,
		const Vector<double> &VelocitySpace_error_per_cell)
	{
		Assert(class_initialized == true,ExcMessage("Please initialize the post proc class"));
		if (constants.print_all)
		{
			if (constants.print_solution)
				print_solution_to_file(triangulation,solution,S_half_inv,dof_handler);

			if(constants.print_error)
				print_error_to_file(triangulation,solution,dof_handler);

			if(constants.print_exactsolution)
				print_exactsolution_to_file(triangulation,S_half_inv);

			if(constants.print_velocity_space_error)
				print_VelocitySpace_error_to_file(triangulation,VelocitySpace_error_per_cell);

			if(constants.print_fe_index)
				print_fe_index(dof_handler);
		}

		else
		{
   			// only print in the final cycle
			if (present_cycle == total_cycles - 1)
			{
				if(constants.print_solution)
					print_solution_to_file(triangulation,solution,S_half_inv,dof_handler);

				if(constants.print_error)
					print_error_to_file(triangulation,solution,dof_handler);

				if(constants.print_exactsolution)
					print_exactsolution_to_file(triangulation,S_half_inv);

				if(constants.print_velocity_space_error)
					print_VelocitySpace_error_to_file(triangulation,VelocitySpace_error_per_cell);

				if(constants.print_fe_index)
					print_fe_index(dof_handler);
			}
		}

		if (constants.print_convergence_table)
			print_convergence_table_to_file(convergence_table);

	}

	// evaluates the error per cell using gausian quadrature
	template<int dim>
	Vector<double>
	Base_PostProc<dim>::return_error_per_cell(const Vector<double> &solution,
		const unsigned int active_cells,
		const hp::MappingCollection<dim> &mapping,
		const hp::DoFHandler<dim> &dof_handler)
	{
		const unsigned int ngp = constants.p + 1;
		QGauss<dim> quadrature_basic(ngp);			
		hp::QCollection<dim> hp_quadrature;
		
		for (unsigned long int i = 0 ; i < nEqn.size() ; i++)
				hp_quadrature.push_back(quadrature_basic);


		used_qgauss = true;
		Assert(class_initialized == true,ExcMessage("Please initialize the post proc class"));

		// error variable comes from Basics. The following is the component in which
		// we wish to find the error.
        		// error variable comes from Basics. The following is the component in which
		// we wish to find the error.
		Assert(constants.variable_map.size() != 0,ExcMessage("Variable map not initialized"));
		Assert(constants.variable_map_1D.size() != 0,ExcMessage("Variable map not initialized"));

		unsigned int component;
		
		if (dim == 1)
			component = constants.variable_map_1D.find(constants.error_variable)->second;

		if (dim == 2)
			component = constants.variable_map.find(constants.error_variable)->second;

        // error per cell of the domain
		Vector<double> error_per_cell(active_cells);      

        ComponentSelectFunction<dim> weight(component,max_equations);                              // used to compute only the error in theta

        // computation of L2 error
        VectorTools::integrate_difference (mapping,dof_handler,solution,
        								  *base_exactsolution,
        								   error_per_cell,
        								   hp_quadrature,
        								   VectorTools::L2_norm,
        								   &weight);  

        return(error_per_cell);

    }

    // in the following routine we compute the lift and drag on an object
    template<int dim>
    Vector<double>
    Base_PostProc<dim>::compute_lift_drag(const MappingQ<dim> &mapping,
    									  const FESystem<dim> &finite_element,
    									  const DoFHandler<dim> &dof_handler,
    									  const Sparse_matrix &S_half_inv,
    									  const Vector<double> &solution,
    									  ConvergenceTable &convergence_table,
    									  const unsigned int b_id_surface)
    {
    	AssertDimension(dim,2);
    	Assert(class_initialized == true,ExcMessage("Please initialize the post proc class"));

    	// location of the needed variables and the corresponding conversion factors
    	const unsigned int id_rho = this->constants.variable_map.find("rho")->second; 
    	const unsigned int id_theta = this->constants.variable_map.find("theta")->second;
    	const unsigned int id_sigmaxx = this->constants.variable_map.find("sigmaxx")->second;
    	const unsigned int id_sigmayy = this->constants.variable_map.find("sigmayy")->second;
    	const unsigned int id_sigmaxy = this->constants.variable_map.find("sigmaxy")->second;

    	double rho;					// value of density
    	double theta;				// value of temperature
    	double sigmaxx;				// value of sigma
    	double sigmaxy;				
    	double sigmayy;

    	double Tx; 					// traction force in the x-direction
    	double Ty;					// traction force in the y-direction
    	double pressure_x;			// force due to pressure in the x-direction
    	double pressure_y;			// force due to pressure in the y-direction

    	// conversion factors for different moments
    	const double fac_rho = 1.0;
    	const double fac_theta = -sqrt(2.0/3.0);
    	const double fac_sigma = sqrt(2);

    	// values of the lift and drag coefficients
    	Vector<double> lift_drag(2);

      	const UpdateFlags face_update_flags  = update_values
      											| update_q_points
      											| update_JxW_values
      											| update_normal_vectors;


      	lift_drag = 0;
      	QGauss<dim-1> face_quadrature(constants.p+1);
      	const unsigned int total_ngp = face_quadrature.size();
    	FEFaceValues<dim> fe_v_face(mapping,finite_element, face_quadrature, face_update_flags);

    	Vector<double> solution_value(max_equations);
    	typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();

    	for (; cell != endc ; cell++)
    	{
    		for (unsigned int face = 0 ; face < GeometryInfo<dim>::faces_per_cell ; face++)
    			// if the face is at the boundary at which we wish to compute the forces
    			if (cell->face(face)->at_boundary() && cell->face(face)->boundary_id() == b_id_surface)
    			{
    				fe_v_face.reinit(cell,face);
    				 const std::vector<double> &Jacobian_face = fe_v_face.get_JxW_values();
 					 const std::vector<Point<dim>> &quad_points = fe_v_face.get_quadrature_points();

    				for (unsigned int q = 0 ; q < total_ngp ; q++)
    				{
    					Tensor<1,dim> normal_vector = fe_v_face.normal_vector(q);
    					const double nx = normal_vector[0];
    					const double ny = normal_vector[1];

    					// value of the solution at the quadrature point
    					VectorTools::point_value(dof_handler, solution, quad_points[q],solution_value);	

    					// converting back to unsymmetric variables
    					solution_value = matrix_opt.Sparse_matrix_dot_Vector(S_half_inv,solution_value);

    					// pressure force in the x-direction
    					pressure_x = nx * (fac_rho * solution_value(id_rho) + fac_theta * solution_value(id_theta));

    					// pressure force in the y-direction
    					pressure_y = ny * (fac_rho * solution_value(id_rho) + fac_theta * solution_value(id_theta));


    					// traction in the x-direction
    					Tx = fac_sigma * (solution_value(id_sigmaxx) * nx + solution_value(id_sigmaxy) * ny);

    					// traction in the y-direction
    					Ty = fac_sigma * (solution_value(id_sigmaxy) * nx + solution_value(id_sigmayy) * ny);

    					// std::cout << "Solution value: " << solution_value << std::endl;
    					// std::cout << "Qudrature point: " << quad_points[q] << std::endl;
    					// std::cout << "nx: " << nx << " ny: " << ny << std::endl;
    					// std::cout << "pressure_x: " << pressure_x << " pressure_y: " << pressure_y << std::endl;
    					// std::cout << "Tx: " << Tx << " Ty: " << Ty << std::endl;

    					// we now integrate over the surface
    					// drag on the force, i.e the force in the x-direction
    					lift_drag(0) += (pressure_x + Tx) * Jacobian_face[q];

    					// lift on the surface, i.e. the force in the y-direction
    					lift_drag(1) += (pressure_y + Ty) * Jacobian_face[q];

    				}
    			}
    	}

    	std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<Lift and Drag>>>>>>>>>>>>>>>>>>>>>>>>"<<std::endl;
    	std::cout << "Drag: " << lift_drag(0) << " Lift: " << lift_drag(1) << std::endl;
    	std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

    	convergence_table.add_value("#Drag",lift_drag(0));
    	convergence_table.add_value("#Lift",lift_drag(1));

        convergence_table.set_scientific("#Drag",true);
        convergence_table.set_scientific("#Lift",true);

    	return(lift_drag);


    }

}
