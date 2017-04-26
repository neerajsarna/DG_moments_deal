namespace FEM_Solver
{
	// base class for running a system, We just use it for storing the eact solution and also for defining certain
	// basic routines
	using namespace dealii;

	template<int dim>
	class
	Run_Problem
	{
	public:
		Run_Problem(ExactSolution::Base_ExactSolution<dim> *exact_solution);

            // the exact solution will correspond to a particular system
		ExactSolution::Base_ExactSolution<dim> *base_exactsolution;

           	// x in final Ax = b
		Vector<double> solution;

			 // Data for post proc
		ConvergenceTable convergence_table;
		std::vector<double> error_per_itr;

         // the following quantity is a measure of how good our FE solution satisfies the weak form of the equations
		Vector<double> residual;

		void distribute_dof_allocate_matrix(DoFHandler<dim> &dof_handler,FESystem<dim> &finite_element,TrilinosWrappers::SparseMatrix &global_matrix);
		void distribute_dof_allocate_matrix(hp::DoFHandler<dim> &dof_handler,hp::FECollection<dim> &finite_element,TrilinosWrappers::SparseMatrix &global_matrix);

		void allocate_vectors(DoFHandler<dim> &dof_handler,Vector<double> &solution,Vector<double> &system_rhs,
							Vector<double> &residual);

		void allocate_vectors(hp::DoFHandler<dim> &dof_handler,Vector<double> &solution,Vector<double> &system_rhs,
							  Vector<double> &residual);



	};

	template<int dim>
	Run_Problem<dim>::Run_Problem(ExactSolution::Base_ExactSolution<dim> *exact_solution)
	:
	base_exactsolution(exact_solution)
	{}

	template<int dim>
	void 
	Run_Problem<dim>::distribute_dof_allocate_matrix(DoFHandler<dim> &dof_handler,FESystem<dim> &finite_element,
													 TrilinosWrappers::SparseMatrix &global_matrix)
	{
		dof_handler.distribute_dofs(finite_element);

        // the vector which stores the residual 
		DynamicSparsityPattern dsp(dof_handler.n_dofs(),dof_handler.n_dofs());


	//std::cout << "making flux sparsity pattern " << std::endl;
		DoFTools::make_flux_sparsity_pattern (dof_handler, dsp);

		global_matrix.reinit(dsp);  
	}

	template<int dim>
	void 
	Run_Problem<dim>::distribute_dof_allocate_matrix(hp::DoFHandler<dim> &dof_handler,hp::FECollection<dim> &finite_element,
													TrilinosWrappers::SparseMatrix &global_matrix)
	{
		dof_handler.distribute_dofs(finite_element);

        // the vector which stores the residual 
		DynamicSparsityPattern dsp(dof_handler.n_dofs(),dof_handler.n_dofs());


	//std::cout << "making flux sparsity pattern " << std::endl;
		DoFTools::make_flux_sparsity_pattern (dof_handler, dsp);

		global_matrix.reinit(dsp);  
	}



template<int dim>
	void 
	Run_Problem<dim>::allocate_vectors(DoFHandler<dim> &dof_handler,Vector<double> &solution,Vector<double> &system_rhs,
										Vector<double> &residual)
	{
		solution.reinit(dof_handler.n_dofs());
		system_rhs.reinit(dof_handler.n_dofs());
		residual.reinit(dof_handler.n_dofs());
	}

template<int dim>
	void 
	Run_Problem<dim>::allocate_vectors(hp::DoFHandler<dim> &dof_handler,Vector<double> &solution,Vector<double> &system_rhs,
										Vector<double> &residual)
	{
		solution.reinit(dof_handler.n_dofs());
		system_rhs.reinit(dof_handler.n_dofs());
		residual.reinit(dof_handler.n_dofs());
	}

template<int dim>
	class
	Run_Problem_FE:public Assembly_Manager_FE<dim>,
					public Run_Problem<dim>
	{
		public:
			Run_Problem_FE(const std::string &output_file_name,
						const constant_numerics &constants,
						std::vector<Develop_System::System<dim>> &equation_info,
                        ExactSolution::Base_ExactSolution<dim> *exact_solution,
                        const std::vector<int> &nEqn,
                        const std::vector<int> &nBC);

			void run();

	};

	template<int dim>
	Run_Problem_FE<dim>::Run_Problem_FE(const std::string &output_file_name,
						const constant_numerics &constants,
						std::vector<Develop_System::System<dim>> &equation_info,
                        ExactSolution::Base_ExactSolution<dim> *exact_solution,
                        const std::vector<int> &nEqn,
                        const std::vector<int> &nBC)
	:
	Assembly_Manager_FE<dim>(output_file_name,
					constants,equation_info,
                	nEqn,nBC),
	Run_Problem<dim>(exact_solution)
	{}

	template<int dim>
	void 
	Run_Problem_FE<dim>::run()
	{
			//
	const unsigned int refine_cycles = this->constants.refine_cycles;

	this->error_per_itr.resize(refine_cycles);

    TimerOutput timer (std::cout, TimerOutput::summary,
                	   TimerOutput::wall_times);

	for (unsigned int cycle = 0 ; cycle < refine_cycles ; cycle ++)
	{
		this->fe_data_structure.print_mesh_info();


	
		AssertDimension(this->error_per_itr.size(),refine_cycles);

		//now we distribute the dofs and allocate the memory for the global matrix. We have 
		// already created the mesh so we can directly distribute the degrees of freedom now.

		std::cout << "Distributing dof " << std::endl;
		timer.enter_subsection("Dof Distribution");
		this->distribute_dof_allocate_matrix(this->fe_data_structure.dof_handler,this->fe_data_structure.finite_element,this->global_matrix);
		this->allocate_vectors(this->fe_data_structure.dof_handler,this->solution,this->system_rhs,this->residual);
		timer.leave_subsection();

		std::cout << "#CELLS " << this->fe_data_structure.triangulation.n_active_cells() << std::endl;
		std::cout << "#DOFS " << this->fe_data_structure.dof_handler.n_dofs() << std::endl;
		std::cout << "Memory by dof handler(Gb) " << 
					this->fe_data_structure.dof_handler.memory_consumption()/pow(10,9)<< std::endl;	

		// the following routine assembles
		std::cout << "Assembling" << std::endl;
		timer.enter_subsection("Assembly");
		Assert(this->constants.assembly_type!= manuel,ExcMessage("Manuel Assembly not supported for this problem."));
		this->assemble_system_meshworker();
		timer.leave_subsection();

		
		// we initialize the object which will solve our system
		// We do int the following way so as to keep the solver independent of all the other implementations.
		// This makes the code highly reusable. So one can directly copy the following class and use it somewhere
		// else if one wants to.

		timer.enter_subsection("Linear Solver");
		LinearSolver::LinearSolver linear_solver;
		std::cout << "Preparing data for pardiso " << std::endl;
		linear_solver.develop_pardiso_data(this->global_matrix);
		this->residual = linear_solver.solve_with_pardiso(this->system_rhs,this->solution);
		timer.leave_subsection();

		timer.enter_subsection("Post Processing");



		PostProc::Base_PostProc<dim> postproc(this->constants,this->base_exactsolution,
											 this->nEqn,this->nBC);

		postproc.reinit(this->fe_data_structure.dof_handler);

		// // now we compute the error due to computation
		postproc.error_evaluation_QGauss(this->solution,
										 this->fe_data_structure.triangulation.n_active_cells(),
										  this->error_per_itr[cycle],
										GridTools::maximal_cell_diameter(this->fe_data_structure.triangulation),
										this->convergence_table,
										this->residual.l2_norm(),
										this->fe_data_structure.mapping,
										this->fe_data_structure.dof_handler);


		// postproc.compute_lift_drag(this->fe_data_structure.mapping,
  //   							   this->fe_data_structure.finite_element,
  //   							   this->fe_data_structure.dof_handler,
  //   							   this->system_info[0].base_tensorinfo.S_half_inv,
  //   							   this->solution,
  //   							   this->convergence_table,
  //   								2);

		// we sent the symmetrizer corresponding to the maximum moment system which we are solving for
		// In this case, since we are only considering a single moment system therefore this value corresponds to zero.
		
		postproc.print_options(this->fe_data_structure.triangulation,this->solution,cycle,refine_cycles,
							  this->convergence_table,
							this->system_info[0].base_tensorinfo.S_half_inv,
								this->fe_data_structure.dof_handler);		
		
		timer.leave_subsection();

		// Grid refinement should be done in the end.
		this->fe_data_structure.refinement_handling(cycle,refine_cycles);

	}
	}


	template<int dim>
	class
	Run_Problem_hp_FE:public Assembly_Manager_hp_FE<dim>,
					  public Run_Problem<dim>
	{
		public:
			Run_Problem_hp_FE(const std::string &output_file_name,
						const constant_numerics &constants,
						std::vector<Develop_System::System<dim>> &equation_info,
                        ExactSolution::Base_ExactSolution<dim> *exact_solution,
                        const std::vector<int> &nEqn,
                        const std::vector<int> &nBC);

			void run();

	};

	template<int dim>
	Run_Problem_hp_FE<dim>::Run_Problem_hp_FE(const std::string &output_file_name,
						const constant_numerics &constants,
						std::vector<Develop_System::System<dim>> &equation_info,
                        ExactSolution::Base_ExactSolution<dim> *exact_solution,
                        const std::vector<int> &nEqn,
                        const std::vector<int> &nBC)
	:
	Assembly_Manager_hp_FE<dim>(output_file_name,
					constants,equation_info,
                	nEqn,nBC),
	Run_Problem<dim>(exact_solution)
	{}

	template<int dim>
	void 
	Run_Problem_hp_FE<dim>::run()
	{
					// total number of refinement cycles in the physical space
			const int refine_cycles_h = this->constants.refine_cycles;

			// total number of refinement cycles in the velocity space
			const int refine_cycles_c = this->constants.refine_cycles_c;

			this->error_per_itr.resize(refine_cycles_c * refine_cycles_h);

			Assert(this->constants.problem_type != periodic , ExcMessage("Use different routine for periodic boudnary conditions"));

			TimerOutput timer (std::cout, TimerOutput::summary,
                	   TimerOutput::wall_times);

			this->hp_fe_data_structure.print_mesh_info();
			// we first do h-refinement
			for (int cycle_h = 0 ; cycle_h < refine_cycles_h ; cycle_h ++)
			{
				for (int cycle_c = 0 ; cycle_c < refine_cycles_c ; cycle_c ++)
				{
					// allocate the index for every cell
					timer.enter_subsection("Dof Distribution");
					std::cout << "Dof distirubtion" << std::endl;
					this->hp_fe_data_structure.allocate_fe_index_distance_center(cycle_c,refine_cycles_c);

					// distribute the degrees of freedom for the different fe indices which have been distributed
					this->distribute_dof_allocate_matrix(this->hp_fe_data_structure.dof_handler,
														 this->hp_fe_data_structure.finite_element,
														 this->global_matrix);


					this->allocate_vectors(this->hp_fe_data_structure.dof_handler,this->solution,this->system_rhs,
										   this->residual);

					timer.leave_subsection();
					Assert(this->constants.assembly_type == manuel,ExcMessage("Only manuel assembly allowed"));
					
					timer.enter_subsection("Assembly");
					std::cout << "Assembly " << std::endl;
					this->assemble_system_manuel();
					timer.leave_subsection();

					LinearSolver::LinearSolver linear_solver;
					std::cout << "Linear Solver" << std::endl;
					timer.enter_subsection("Linear Solver");
					linear_solver.develop_pardiso_data(this->global_matrix);
					this->residual = linear_solver.solve_with_pardiso(this->system_rhs,this->solution);

					timer.leave_subsection();

					std::cout << "post processing " << std::endl;
					timer.enter_subsection("post processing");


					this->hp_fe_data_structure.compute_equilibrium_deviation(this->ngp,
                                                				this->nEqn,
                                                				this->hp_fe_data_structure.triangulation,
                                                				this->solution,
                                                				cycle_h);


					PostProc::Base_PostProc<dim> postproc(this->constants,
														 this->base_exactsolution,
														this->nEqn,this->nBC);

					postproc.reinit(this->hp_fe_data_structure.dof_handler);


					// // now we compute the error due to computation
					postproc.error_evaluation_QGauss(this->solution,
										 			this->hp_fe_data_structure.triangulation.n_active_cells(),
										   			this->error_per_itr[cycle_h *(refine_cycles_c) + cycle_c],
													GridTools::maximal_cell_diameter(this->hp_fe_data_structure.triangulation),
													this->convergence_table,
													this->residual.l2_norm(),
													this->hp_fe_data_structure.mapping,
													this->hp_fe_data_structure.dof_handler);


					postproc.print_options(this->hp_fe_data_structure.triangulation,this->solution,cycle_c,refine_cycles_c,
										   this->convergence_table,
										  this->system_info[this->hp_fe_data_structure.max_fe_index].base_tensorinfo.S_half_inv,
										  this->hp_fe_data_structure.dof_handler,
										  this->hp_fe_data_structure.VelocitySpace_error_per_cell);




				  timer.leave_subsection();

				}	
				
				// Grid refinement should be done in the end.
				//this->hp_fe_data_structure.refinement_handling(cycle_h,refine_cycles_h);			
			}

	}

}
