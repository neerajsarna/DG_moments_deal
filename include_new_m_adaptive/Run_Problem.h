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

		MatrixOpt::Base_MatrixOpt matrix_opt;


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
	Run_Problem<dim>::allocate_vectors(DoFHandler<dim> &dof_handler,
										Vector<double> &solution,
										Vector<double> &system_rhs,
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
                        const std::vector<int> &nBC,
                        const unsigned int system_to_solve,
                        Triangulation<dim> &triangulation);


			void run(MeshGenerator::Base_MeshGenerator<dim> &Mesh_Info);

	};

	template<int dim>
	Run_Problem_FE<dim>::Run_Problem_FE(const std::string &output_file_name,
						const constant_numerics &constants,
						std::vector<Develop_System::System<dim>> &equation_info,
                        ExactSolution::Base_ExactSolution<dim> *exact_solution,
                        const std::vector<int> &nEqn,
                        const std::vector<int> &nBC,
                        const unsigned int system_to_solve,
                        Triangulation<dim> &triangulation)
	:
	Assembly_Manager_FE<dim>(output_file_name,
					constants,equation_info,
                	nEqn,nBC,system_to_solve,triangulation),
	Run_Problem<dim>(exact_solution)
	{}

	template<int dim>
	void 
	Run_Problem_FE<dim>::run(MeshGenerator::Base_MeshGenerator<dim> &Mesh_Info)
	{
	
	const unsigned int refine_cycles = this->constants.refine_cycles;

	this->error_per_itr.resize(refine_cycles);

    TimerOutput timer (std::cout, TimerOutput::summary,
                	   TimerOutput::wall_times);

	for (unsigned int cycle = 0 ; cycle < refine_cycles ; cycle ++)
	{
		Mesh_Info.print_mesh_info();


	
		AssertDimension(this->error_per_itr.size(),refine_cycles);

		//now we distribute the dofs and allocate the memory for the global matrix. We have 
		// already created the mesh so we can directly distribute the degrees of freedom now.

		std::cout << "Distributing dof " << std::endl;
		timer.enter_subsection("Dof Distribution");
		this->distribute_dof_allocate_matrix(this->fe_data_structure.dof_handler,
											this->fe_data_structure.finite_element,this->global_matrix);

		std::cout << "#Dofs " << this->fe_data_structure.dof_handler.n_dofs() << std::endl;
		this->allocate_vectors(this->fe_data_structure.dof_handler,this->solution,this->system_rhs,this->residual);
		timer.leave_subsection();


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



		PostProc::Base_PostProc<dim> postproc(this->constants,this->base_exactsolution);

		postproc.reinit(this->fe_data_structure.dof_handler);

		// // now we compute the error due to computation
		postproc.error_evaluation_QGauss(this->solution,
										 Mesh_Info.triangulation.n_active_cells(),
										  this->error_per_itr[cycle],
										GridTools::maximal_cell_diameter(Mesh_Info.triangulation),
										this->convergence_table,
										this->residual.l2_norm(),
										this->fe_data_structure.mapping,
										this->fe_data_structure.dof_handler,
										this->nEqn);


		// postproc.compute_lift_drag(this->fe_data_structure.mapping,
  //   							   this->fe_data_structure.finite_element,
  //   							   this->fe_data_structure.dof_handler,
  //   							   this->system_info[0].base_tensorinfo.S_half_inv,
  //   							   this->solution,
  //   							   this->convergence_table,
  //   								2);

		// we sent the symmetrizer corresponding to the maximum moment system which we are solving for
		// In this case, since we are only considering a single moment system therefore this value corresponds to zero.
		
		postproc.print_options(Mesh_Info.triangulation,this->solution,cycle,refine_cycles,
							  this->convergence_table,
							this->system_info.base_tensorinfo.S_half_inv,
								this->fe_data_structure.dof_handler,this->nEqn);		
		
		timer.leave_subsection();

		// Grid refinement should be done in the end.
		Mesh_Info.refinement_handling(cycle,refine_cycles);

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
                        const std::vector<int> &nBC,
                        Triangulation<dim> &triangulation);

			// different ways to adapt
			// in the following routine we compute the deviation of the distribution function from
			// the distribution function of the previous moment theory. For e.g. for the lowest moment theory
			// we compute the deviation from the the maxwellian, for higher moment theory we compute the deviation from
			// the previous one.
			void run_distribution_deviation(MeshGenerator::Base_MeshGenerator<dim> &Mesh_Info);

			// compare the solution obtained with a higher order reference solution and then refine
			void run_higher_order_reference(DoFHandler<dim> &dof_handler_reference,
											const Vector<double> &solution_reference,
											MeshGenerator::Base_MeshGenerator<dim> &Mesh_Info);

	};

	template<int dim>
	Run_Problem_hp_FE<dim>::Run_Problem_hp_FE(const std::string &output_file_name,
						const constant_numerics &constants,
						std::vector<Develop_System::System<dim>> &equation_info,
                        ExactSolution::Base_ExactSolution<dim> *exact_solution,
                        const std::vector<int> &nEqn,
                        const std::vector<int> &nBC,
                        Triangulation<dim> &triangulation)
	:
	Assembly_Manager_hp_FE<dim>(output_file_name,
					constants,equation_info,
                	nEqn,nBC,triangulation),
	Run_Problem<dim>(exact_solution)
	{}

	template<int dim>
	void 
	Run_Problem_hp_FE<dim>::run_distribution_deviation(MeshGenerator::Base_MeshGenerator<dim> &Mesh_Info)
	{
					// total number of refinement cycles in the physical space
			const int refine_cycles_h = this->constants.refine_cycles;

			// total number of refinement cycles in the velocity space
			const int refine_cycles_c = this->constants.refine_cycles_c;

			this->error_per_itr.resize(refine_cycles_c * refine_cycles_h);

			Assert(this->constants.problem_type != periodic , ExcMessage("Use different routine for periodic boudnary conditions"));

			TimerOutput timer (std::cout, TimerOutput::summary,
                	   TimerOutput::wall_times);

			Mesh_Info.print_mesh_info();

			// we only allow for one refinement cycle in the space dimension
			AssertDimension(refine_cycles_h,1);
			AssertDimension(refine_cycles_c,(int)this->nEqn.size());

			// we first do h-refinement
			for (int cycle_h = 0 ; cycle_h < refine_cycles_h ; cycle_h ++)
			{
				for (int cycle_c = 0 ; cycle_c < refine_cycles_c ; cycle_c ++)
				{
					// allocate the index for every cell
					timer.enter_subsection("Dof Distribution");
					std::cout << "Dof distirubtion" << std::endl;
					this->hp_fe_data_structure.allocate_fe_index_distribution_deviation(cycle_c);

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


					this->hp_fe_data_structure.compute_distribution_deviation(this->ngp,
                                                				this->nEqn,
                                                				Mesh_Info.triangulation,
                                                				this->solution,
                                                				cycle_c);


					PostProc::Base_PostProc<dim> postproc(this->constants,
														 this->base_exactsolution);

					postproc.reinit(this->hp_fe_data_structure.dof_handler);


					// // now we compute the error due to computation
					postproc.error_evaluation_QGauss(this->solution,
										 			Mesh_Info.triangulation.n_active_cells(),
										   			this->error_per_itr[cycle_h *(refine_cycles_c) + cycle_c],
													GridTools::maximal_cell_diameter(Mesh_Info.triangulation),
													this->convergence_table,
													this->residual.l2_norm(),
													this->hp_fe_data_structure.mapping,
													this->hp_fe_data_structure.dof_handler,
													this->nEqn);


					postproc.print_options(Mesh_Info.triangulation,this->solution,cycle_c,refine_cycles_c,
										   this->convergence_table,
										  this->system_info[this->hp_fe_data_structure.max_fe_index].base_tensorinfo.S_half_inv,
										  this->hp_fe_data_structure.dof_handler,
										  this->hp_fe_data_structure.VelocitySpace_error_per_cell,
										  this->nEqn);




				  timer.leave_subsection();

				}	
				
				// Grid refinement should be done in the end.
				//this->hp_fe_data_structure.refinement_handling(cycle_h,refine_cycles_h);			
			}

	}

	template<int dim>
	void 
	Run_Problem_hp_FE<dim>::run_higher_order_reference(DoFHandler<dim> &dof_handler_reference,
													   const Vector<double> &solution_reference,
													   MeshGenerator::Base_MeshGenerator<dim> &Mesh_Info)
	{
			// total number of refinement cycles in the physical space
			const int refine_cycles_h = this->constants.refine_cycles;
			const int refine_cycles_c = this->constants.refine_cycles_c;


			this->error_per_itr.resize(refine_cycles_c * refine_cycles_h);

			TimerOutput timer (std::cout, TimerOutput::summary,
                	   TimerOutput::wall_times);

			Mesh_Info.print_mesh_info();

			// we only allow for one refinement cycle in the space dimension
			AssertDimension(refine_cycles_h,1);

			

			// we first do h-refinement
			for (int cycle_h = 0 ; cycle_h < refine_cycles_h ; cycle_h ++)
			{
				for (int cycle_c = 0 ; cycle_c < refine_cycles_c ; cycle_c ++)
				{
					std::cout << "M refinement Cycle.." << cycle_c <<  std::endl;

					// allocate the index for every cell
					timer.enter_subsection("Dof Distribution");
					std::cout << "Dof distirubtion" << std::endl;
					this->hp_fe_data_structure.allocate_fe_index_distribution_deviation(cycle_c);

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


					this->hp_fe_data_structure.compute_error_comparitive(this->solution,
                                          								 this->constants,
                                          								dof_handler_reference,
                                          								solution_reference);


					PostProc::Base_PostProc<dim> postproc(this->constants,
														 this->base_exactsolution);

					postproc.reinit(this->hp_fe_data_structure.dof_handler);


					// // now we compute the error due to computation
					postproc.error_evaluation_QGauss(this->solution,
										 			Mesh_Info.triangulation.n_active_cells(),
										   			this->error_per_itr[cycle_h *(refine_cycles_c) + cycle_c],
													GridTools::maximal_cell_diameter(Mesh_Info.triangulation),
													this->convergence_table,
													this->residual.l2_norm(),
													this->hp_fe_data_structure.mapping,
													this->hp_fe_data_structure.dof_handler,
													this->nEqn);


					postproc.print_options(Mesh_Info.triangulation,this->solution,cycle_c,refine_cycles_c,
										   this->convergence_table,
										  this->system_info[this->hp_fe_data_structure.max_fe_index].base_tensorinfo.S_half_inv,
										  this->hp_fe_data_structure.dof_handler,
										  this->hp_fe_data_structure.VelocitySpace_error_per_cell,
										  this->nEqn);




				  timer.leave_subsection();

				}	
				
			}
	}

	// time stepping using meshworker
	template<int dim>
	class
	Run_Problem_FE_Time_Stepping:public Assembly_Manager_FE<dim>,
								 public Run_Problem<dim>
	{
		public:
			Run_Problem_FE_Time_Stepping(const std::string &output_file_name,
						const constant_numerics &constants,
						std::vector<Develop_System::System<dim>> &equation_info,
                        ExactSolution::Base_ExactSolution<dim> *exact_solution,
                        const std::vector<int> &nEqn,
                        const std::vector<int> &nBC,
                        const unsigned int system_to_solve,
                        Triangulation<dim> &triangulation);


			// the following function prescribes the initial conditions
			class
			initial_conditions:public Function<dim>
			{
				public:
					initial_conditions(const unsigned int  &nEqn);
					virtual void vector_value(const Point<dim> &p,Vector<double> &value) const ;
			};

			// in the following routine we check mass conservation
			double mean_value(const unsigned int active_cells,const unsigned int comp);
			const double delta_t = 0.1;
			void run(MeshGenerator::Base_MeshGenerator<dim> &Mesh_Info);
			void compute_energy_bound();

	};

	template<int dim>
	Run_Problem_FE_Time_Stepping<dim>::Run_Problem_FE_Time_Stepping(const std::string &output_file_name,
						const constant_numerics &constants,
						std::vector<Develop_System::System<dim>> &equation_info,
                        ExactSolution::Base_ExactSolution<dim> *exact_solution,
                        const std::vector<int> &nEqn,
                        const std::vector<int> &nBC,
                        const unsigned int system_to_solve,
                        Triangulation<dim> &triangulation)
	:
	Assembly_Manager_FE<dim>(output_file_name,
					constants,equation_info,
                	nEqn,nBC,system_to_solve,triangulation),
	Run_Problem<dim>(exact_solution)
	{
		// only works for finite volume scheme, does not have higher order time stepping scheme
		//AssertDimension(constants.p,0);

		// has not been implemented for more than one system
		AssertDimension(nEqn.size(),1);
	}

	template<int dim>
	Run_Problem_FE_Time_Stepping<dim>::initial_conditions::initial_conditions(const unsigned int &nEqn)
	:
	Function<dim>(nEqn)
	{;};

	//returns the mean of the comp component of the solution
	template<int dim>
	double 
	Run_Problem_FE_Time_Stepping<dim>::mean_value(const unsigned int active_cells,const unsigned int comp)
	{

		AssertIndexRange(comp,(unsigned int)this->nEqn);
		// we select the component corresponding to density
		unsigned int component= comp;
		

		const unsigned int ngp = this->constants.p + 1;
        // error per cell of the domain
		Vector<double> density_per_cell(active_cells);      


        ComponentSelectFunction<dim> weight(component,this->nEqn);       
        // used to compute only the error in theta

        // computation of L2 error
        VectorTools::integrate_difference (this->fe_data_structure.mapping,this->fe_data_structure.dof_handler,
        									this->solution,
        								   ZeroFunction<dim>(this->nEqn),
        									density_per_cell,
        									QGauss<dim>(ngp),
        									VectorTools::mean,
        									&weight); 

        // the above routine computes the negative of the mean therefore we multiply by a negative
        density_per_cell *= -1;


        return(density_per_cell.mean_value() * density_per_cell.size());
	}

	template<int dim>
	void
	Run_Problem_FE_Time_Stepping<dim>::initial_conditions::vector_value(const Point<dim> &p,Vector<double> &value) const
	{
		const double x = p(0);

		// zero initial conditions to all the variables
		unsigned int ID_rho;
		unsigned int ID_vx;
		unsigned int ID_vy;
		unsigned int ID_theta;

		value = 0;

	}


	template<int dim>
	void 
	Run_Problem_FE_Time_Stepping<dim>::run(MeshGenerator::Base_MeshGenerator<dim> &Mesh_Info)
	{
	// the total number of refinement cycles 
	const unsigned int refine_cycles = this->constants.refine_cycles;

	// does not support grid refinement right now
	AssertDimension(refine_cycles,1);

	this->error_per_itr.resize(refine_cycles);

    TimerOutput timer (std::cout, TimerOutput::summary,
                	   TimerOutput::wall_times);


   	// we first assemble the differential operator
   	Mesh_Info.print_mesh_info();



   	AssertDimension(this->error_per_itr.size(),refine_cycles);

		//now we distribute the dofs and allocate the memory for the global matrix. We have 
		// already created the mesh so we can directly distribute the degrees of freedom now.

   	std::cout << "Distributing dof " << std::endl;
   	timer.enter_subsection("Dof Distribution");
   	this->distribute_dof_allocate_matrix(this->fe_data_structure.dof_handler,this->fe_data_structure.finite_element,this->global_matrix);
   	this->allocate_vectors(this->fe_data_structure.dof_handler,this->solution,this->system_rhs,this->residual);
   	timer.leave_subsection();


		// the following routine assembles
   	std::cout << "Assembling" << std::endl;
   	timer.enter_subsection("Assembly");
   	Assert(this->constants.assembly_type!= manuel,ExcMessage("Manuel Assembly not supported for this problem."));
   	this->assemble_system_meshworker();
   	timer.leave_subsection();

   	std::cout << "Computing intial conditions " << std::endl;
   	initial_conditions initial_data(this->nEqn);

   	VectorTools::interpolate(this->fe_data_structure.mapping,
						     this->fe_data_structure.dof_handler,
							  initial_data,
							 this->solution);

   	// we compute the initial density
   	const double initial_density = mean_value(Mesh_Info.triangulation.n_active_cells(),0);
   	// deviation from the steady state value
   	Vector<double> residual_steady_state(this->solution.size());
   	Vector<double> new_solution(this->solution.size());
   	Vector<double> new_solution2(this->solution.size());
   	residual_steady_state = 100;

   	std::cout << "Time stepping " << std::endl;
   	fflush(stdout);

   	// counts the number of time steps
   	unsigned int counter = 0;

   	PostProc::Base_PostProc<dim> postproc(this->constants,this->base_exactsolution);

	postproc.reinit(this->fe_data_structure.dof_handler);

	FILE *energy_growth;
	energy_growth = fopen("energy_growth","w+");

	double present_time = 0;

	compute_energy_bound();

   	while (residual_steady_state.l2_norm() > 1e-8)
   	{
   		counter++;

   		double initial_energy;
   		double final_energy;

   		// forward euler step
   		// new_solution = delta_t * system_rhs
   		// if (counter%100 == 0 )
   		// 	initial_energy = postproc.compute_energy(this->solution,
   		// 											 Mesh_Info.triangulation.n_active_cells(),
					// 								 this->fe_data_structure.mapping, 
					// 								 this->fe_data_structure.dof_handler,this->nEqn);

   		// update with the solution rhs
   		new_solution.equ(delta_t,this->system_rhs);
   		// update with the spatial derivative
   		new_solution.add(-delta_t,
   			Run_Problem<dim>::matrix_opt.Sparse_matrix_dot_Vector(this->global_matrix,this->solution));
   		// update with the previous solution
   		new_solution += this->solution;

   		// second update
   		// update with the previous solution
   		// THIS is incorrect
   		new_solution2.equ(3.0/4.0,this->solution);
   		new_solution2.add(1.0/4.0,new_solution);
   		new_solution2.add(delta_t * 1.0/4.0,this->system_rhs);
   		new_solution2.add(-delta_t * 1.0/4,Run_Problem<dim>::matrix_opt.Sparse_matrix_dot_Vector(this->global_matrix,new_solution));
   

   		// if (counter %100 == 0)
   		// 	final_energy = postproc.compute_energy(new_solution2,
   		// 											Mesh_Info.triangulation.n_active_cells(),
					// 										this->fe_data_structure.mapping, 
					// 										this->fe_data_structure.dof_handler,
					// 										this->nEqn);

   		residual_steady_state = Run_Problem<dim>::matrix_opt.Sparse_matrix_dot_Vector(this->global_matrix,new_solution2);
   		residual_steady_state -= this->system_rhs;

   		// update the old solution
   		this->solution = new_solution2;

   		present_time +=delta_t;

   		if (counter % 100 == 0)
   		{
   			const double current_density = mean_value(Mesh_Info.triangulation.n_active_cells(),0);

   			std::cout << "Residual " << residual_steady_state.l2_norm() << 
   					" Solution Norm " << new_solution2.l2_norm() <<
   					" Deviation in mass " << fabs(initial_density-current_density) <<  std::endl;

   			// printf("energy growth rate: %e \n",(final_energy-initial_energy)/delta_t);
   			// fprintf(energy_growth, "%f %0.15f\n",present_time,(final_energy-initial_energy)/delta_t);

   		}

   	}



		// // now we compute the error due to computation
		postproc.error_evaluation_QGauss(this->solution,
										 Mesh_Info.triangulation.n_active_cells(),
										  this->error_per_itr[0],
										GridTools::maximal_cell_diameter(Mesh_Info.triangulation),
										this->convergence_table,
										this->residual.l2_norm(),
										this->fe_data_structure.mapping,
										this->fe_data_structure.dof_handler,
										this->nEqn);

		postproc.print_options(Mesh_Info.triangulation,this->solution,0,1,
							  this->convergence_table,
							this->system_info.base_tensorinfo.S_half_inv,
								this->fe_data_structure.dof_handler,this->nEqn);		


		fclose(energy_growth);

	}

	// the following routines computes the energy bound 
	template<>
	void 
	Run_Problem_FE_Time_Stepping<1>::compute_energy_bound()
	{

		Tensor<1,1> normal_vector;
		Tensor<1,1> boundary_point;
		unsigned int b_id;
		Vector<double> bc_rhs(this->nBC);

		// boundary at the left hand side
		b_id = 101;
		normal_vector[0] = -1.0;
		boundary_point[0] = -0.5;

		this->system_info.bcrhs_inflow.BCrhs(boundary_point,
						  			normal_vector,
						  			bc_rhs,
						  			b_id);


		// boundary at the right hand side
		b_id = 102;
		normal_vector[0] = 1.0;
		boundary_point[0] = 0.5;

		this->system_info.bcrhs_inflow.BCrhs(boundary_point,
						  			normal_vector,
						  			bc_rhs,
						  			b_id);		


	}


	// through the following class we develop a reference solution 
template<int dim>
	class
	Develop_Reference:public Assembly_Manager_FE<dim>,
					public Run_Problem<dim>
	{
		public:
			Develop_Reference(const std::string &output_file_name,
						const constant_numerics &constants,
						std::vector<Develop_System::System<dim>> &equation_info,
                        ExactSolution::Base_ExactSolution<dim> *exact_solution,
                        const std::vector<int> &nEqn,
                        const std::vector<int> &nBC,
                        const unsigned int system_to_solve,
                        Triangulation<dim> &triangulation);


			void run(MeshGenerator::Base_MeshGenerator<dim> &Mesh_Info);

	};

	template<int dim>
	Develop_Reference<dim>::Develop_Reference(const std::string &output_file_name,
						const constant_numerics &constants,
						std::vector<Develop_System::System<dim>> &equation_info,
                        ExactSolution::Base_ExactSolution<dim> *exact_solution,
                        const std::vector<int> &nEqn,
                        const std::vector<int> &nBC,
                        const unsigned int system_to_solve,
                        Triangulation<dim> &triangulation)
	:
	Assembly_Manager_FE<dim>(output_file_name,
					constants,equation_info,
                	nEqn,nBC,system_to_solve,triangulation),
	Run_Problem<dim>(exact_solution)
	{}

	template<int dim>
	void 
	Develop_Reference<dim>::run(MeshGenerator::Base_MeshGenerator<dim> &Mesh_Info)
	{
		TimerOutput timer (std::cout, TimerOutput::summary,
                	   TimerOutput::wall_times);


		std::cout << "Developing reference solution........" << std::endl;
		Mesh_Info.print_mesh_info();

		std::cout << "Distributing dof " << std::endl;
		timer.enter_subsection("Dof Distribution");
		this->distribute_dof_allocate_matrix(this->fe_data_structure.dof_handler,this->fe_data_structure.finite_element,this->global_matrix);
		this->allocate_vectors(this->fe_data_structure.dof_handler,this->solution,this->system_rhs,this->residual);
		timer.leave_subsection();

		std::cout << "#Dofs " << this->fe_data_structure.dof_handler.n_dofs() << std::endl;

		// the following routine assembles
		std::cout << "Assembling" << std::endl;
		timer.enter_subsection("Assembly");

		//The main assembly type can be any but we only use meshworker
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
	}


	template<int dim>
	class
	Run_Problem_Periodic:public Assembly_Manager_Periodic<dim>,
						 public Run_Problem<dim>

	{
		public:
			Run_Problem_Periodic(const std::string &output_file_name,
						const constant_numerics &constants,
						std::vector<Develop_System::System<dim>> &equation_info,
                        ExactSolution::Base_ExactSolution<dim> *exact_solution,
                        const std::vector<int> &nEqn,
                        const std::vector<int> &nBC,
                        const unsigned int system_to_solve,
                        Triangulation<dim> &triangulation);


			void run(MeshGenerator::Base_MeshGenerator<dim> &Mesh_Info);
			void run_time_step(MeshGenerator::Base_MeshGenerator<dim> &Mesh_Info);
			// in the following routine we check mass conservation
			double mean_value(const unsigned int active_cells,const unsigned int comp);
			void distribute_dof_allocate_matrix_periodic_box();
			const double delta_t = 0.1;

	};

	template<int dim>
	Run_Problem_Periodic<dim>::Run_Problem_Periodic(const std::string &output_file_name,
						const constant_numerics &constants,
						std::vector<Develop_System::System<dim>> &equation_info,
                        ExactSolution::Base_ExactSolution<dim> *exact_solution,
                        const std::vector<int> &nEqn,
                        const std::vector<int> &nBC,
                        const unsigned int system_to_solve,
                        Triangulation<dim> &triangulation)
	:
	Assembly_Manager_Periodic<dim>(output_file_name,
					constants,equation_info,
                	nEqn,nBC,system_to_solve,triangulation),
	Run_Problem<dim>(exact_solution)
	{}

	template<int dim>
	void 
	Run_Problem_Periodic<dim>::distribute_dof_allocate_matrix_periodic_box()
	{
		this->fe_data_structure.dof_handler.distribute_dofs(this->fe_data_structure.finite_element);

		DynamicSparsityPattern dsp(this->fe_data_structure.dof_handler.n_dofs(),
								  this->fe_data_structure.dof_handler.n_dofs());


		DoFTools::make_flux_sparsity_pattern (this->fe_data_structure.dof_handler, dsp);

		this->add_periodic_sparsity(dsp);
		this->global_matrix.reinit(dsp);  
	}

	template<int dim>
	void 
	Run_Problem_Periodic<dim>::run(MeshGenerator::Base_MeshGenerator<dim> &Mesh_Info)
	{

	std::cout << "Solving for " << this->nEqn << " Equations" << std::endl;
	

	this->error_per_itr.resize(this->constants.refine_cycles);

	TimerOutput timer (std::cout, TimerOutput::summary,
                   TimerOutput::wall_times);
	
	Assert(this->constants.mesh_type == square_domain 
			&& this->constants.problem_type == periodic,ExcMessage("Incorrect mesh"));
	
	for (int cycle = 0 ; cycle < this->constants.refine_cycles ; cycle++)
	{

		AssertDimension((int)this->error_per_itr.size(),this->constants.refine_cycles);

		timer.enter_subsection("Develop Periodicity");
		fflush(stdout);

		// First we develop the periodic faces using internal functions of dealii.
		this->develop_periodic_faces(this->fe_data_structure.dof_handler);

		// now we construct the required data structure
		this->divide_periodicity();

		timer.leave_subsection();

		//now we distribute the dofs and allocate the memory for the global matrix. We have 
		// already created the mesh so we can directly distribute the degrees of freedom now.
		timer.enter_subsection("Distribute Dof");
		fflush(stdout);
		this->distribute_dof_allocate_matrix_periodic_box();
		this->allocate_vectors(this->fe_data_structure.dof_handler,
								this->solution,this->system_rhs,
								this->residual);

		timer.leave_subsection();

		std::cout << "#Cells " << Mesh_Info.triangulation.n_active_cells() << std::endl;
		std::cout << "#Dofs " << this->fe_data_structure.dof_handler.n_dofs() << std::endl;
		fflush(stdout);

		// cannot use meshworker for the present problem
		timer.enter_subsection("Assemble");
		this->assemble_system_manuel();

		timer.leave_subsection();
        
		// we initialize the object which will solve our system
		// We do int the following way so as to keep the solver independent of all the other implementations.
		// This makes the code highly reusable. So one can directly copy the following class and use it somewhere
		// else is one wants to.
		
		timer.enter_subsection("Linear Solver");
		LinearSolver::LinearSolver linear_solver;
		linear_solver.develop_pardiso_data(this->global_matrix);
		this->residual = linear_solver.solve_with_pardiso(this->system_rhs,this->solution);
		timer.leave_subsection();

		timer.enter_subsection("Post Proc");
		PostProc::Base_PostProc<dim> postproc(this->constants,this->base_exactsolution);

		postproc.reinit(this->fe_data_structure.dof_handler);

		// // now we compute the error due to computation
		postproc.error_evaluation_QGauss(this->solution,
										 Mesh_Info.triangulation.n_active_cells(),
										  this->error_per_itr[cycle],
										GridTools::maximal_cell_diameter(Mesh_Info.triangulation),
										this->convergence_table,
										this->residual.l2_norm(),
										this->fe_data_structure.mapping,
										this->fe_data_structure.dof_handler,
										this->nEqn);


		// we sent the symmetrizer corresponding to the maximum moment system which we are solving for
		// In this case, since we are only considering a single moment system therefore this value corresponds to zero.
		
		postproc.print_options_quad_points(Mesh_Info.triangulation,this->solution,cycle,this->constants.refine_cycles,
							  this->convergence_table,
							  this->system_info.base_tensorinfo.S_half_inv,
							  this->fe_data_structure.dof_handler,
							  this->fe_data_structure.mapping,
							  this->nEqn);		
		
		timer.leave_subsection();

		// Grid refinement should be done in the end.
		Mesh_Info.refinement_handling(cycle,this->constants.refine_cycles);
		
	}
	}

	//returns the mean of the comp component of the solution
	template<int dim>
	double 
	Run_Problem_Periodic<dim>::mean_value(const unsigned int active_cells,const unsigned int comp)
	{

		AssertIndexRange(comp,(unsigned int)this->nEqn);
		// we select the component corresponding to density
		unsigned int component= comp;
		

		const unsigned int ngp = this->constants.p + 1;
        // error per cell of the domain
		Vector<double> density_per_cell(active_cells);      


        ComponentSelectFunction<dim> weight(component,this->nEqn);       
        // used to compute only the error in theta

        // computation of L2 error
        VectorTools::integrate_difference (this->fe_data_structure.mapping,this->fe_data_structure.dof_handler,
        									this->solution,
        								   ZeroFunction<dim>(this->nEqn),
        									density_per_cell,
        									QGauss<dim>(ngp),
        									VectorTools::mean,
        									&weight); 

        // the above routine computes the negative of the mean therefore we multiply by a negative
        density_per_cell *= -1;


        return(density_per_cell.mean_value() * density_per_cell.size());
	}


	template<int dim>
	void 
	Run_Problem_Periodic<dim>::run_time_step(MeshGenerator::Base_MeshGenerator<dim> &Mesh_Info)
	{

	this->error_per_itr.resize(this->constants.refine_cycles);

	TimerOutput timer (std::cout, TimerOutput::summary,
                   TimerOutput::wall_times);
	
	Assert(this->constants.mesh_type == square_domain 
			&& this->constants.problem_type == periodic,ExcMessage("Incorrect mesh"));
	
	for (int cycle = 0 ; cycle < 1 ; cycle++)
	{

		AssertDimension((int)this->error_per_itr.size(),this->constants.refine_cycles);

		timer.enter_subsection("Develop Periodicity");
		fflush(stdout);
		// first we develop the periodic faces using internal functions of dealii
		this->develop_periodic_faces(this->fe_data_structure.dof_handler);

		// now we construct the required data structure
		this->divide_periodicity();

		timer.leave_subsection();

		//now we distribute the dofs and allocate the memory for the global matrix. We have 
		// already created the mesh so we can directly distribute the degrees of freedom now.
		timer.enter_subsection("Distribute Dof");
		fflush(stdout);
		this->distribute_dof_allocate_matrix_periodic_box();
		this->allocate_vectors(this->fe_data_structure.dof_handler,
								this->solution,this->system_rhs,
								this->residual);

		timer.leave_subsection();

		std::cout << "#Cells " << Mesh_Info.triangulation.n_active_cells() << std::endl;
		std::cout << "#Dofs " << this->fe_data_structure.dof_handler.n_dofs() << std::endl;
		fflush(stdout);

		// cannot use meshworker for the present problem
		timer.enter_subsection("Assemble");
		this->assemble_system_manuel();

		timer.leave_subsection();

		// this->matrix_opt.print_dealii_sparse(this->global_matrix,"global_matrix_odd");
		// this->matrix_opt.print_dealii_vector(this->system_rhs,"rhs_odd");
        
		// we initialize the object which will solve our system
		// We do int the following way so as to keep the solver independent of all the other implementations.
		// This makes the code highly reusable. So one can directly copy the following class and use it somewhere
		// else is one wants to.
		

		timer.enter_subsection("Time Stepping");
		this->solution = 1;
		Vector<double> residual_steady_state(this->solution.size());
		Vector<double> new_solution(this->solution.size());
		Vector<double> new_solution2(this->solution.size());
		residual_steady_state = 100;

		const double initial_density = mean_value(Mesh_Info.triangulation.n_active_cells(),0);
		std::cout << "Time stepping " << std::endl;
		fflush(stdout);

   	// counts the number of time steps
		unsigned int counter = 0;


		while (residual_steady_state.l2_norm() > 1e-6)
		{
			counter++;

   		// update with the solution rhs
			new_solution.equ(delta_t,this->system_rhs);
   		// update with the spatial derivative
			new_solution.add(-delta_t,Run_Problem<dim>::matrix_opt.Sparse_matrix_dot_Vector(this->global_matrix,this->solution));
   		// update with the previous solution
			new_solution += this->solution;

   		// second update
   		// update with the previous solution
			// new_solution2.equ(3.0/4.0,this->solution);
			// new_solution2.add(1.0/4.0,new_solution);
			// new_solution2.add(delta_t * 1.0/4.0,this->system_rhs);
			// new_solution2.add(-delta_t * 1.0/4,Run_Problem<dim>::matrix_opt.Sparse_matrix_dot_Vector(this->global_matrix,new_solution));

			this->solution = new_solution;

			residual_steady_state = Run_Problem<dim>::matrix_opt.Sparse_matrix_dot_Vector(this->global_matrix,this->solution);
			residual_steady_state -= this->system_rhs;

   		// update the old solution
			


			if (counter % 10 == 0)
			{
				const double current_density = mean_value(Mesh_Info.triangulation.n_active_cells(),0);

				std::cout << "residual_steady_state " << residual_steady_state.l2_norm() << 
							 " Solution Norm " << this->solution.l2_norm() <<
							 " Deviation in mass" << fabs(initial_density-current_density) <<  std::endl;

			}

		}
		timer.leave_subsection();

		timer.enter_subsection("Post Proc");
		PostProc::Base_PostProc<dim> postproc(this->constants,this->base_exactsolution);

		postproc.reinit(this->fe_data_structure.dof_handler);

		// // now we compute the error due to computation
		postproc.error_evaluation_QGauss(this->solution,
			Mesh_Info.triangulation.n_active_cells(),
			this->error_per_itr[cycle],
			GridTools::maximal_cell_diameter(Mesh_Info.triangulation),
			this->convergence_table,
			this->residual.l2_norm(),
			this->fe_data_structure.mapping,
			this->fe_data_structure.dof_handler,
			this->nEqn);


		// we sent the symmetrizer corresponding to the maximum moment system which we are solving for
		// In this case, since we are only considering a single moment system therefore this value corresponds to zero.
		
		postproc.print_options_quad_points(Mesh_Info.triangulation,this->solution,cycle,this->constants.refine_cycles,
			this->convergence_table,
			this->system_info.base_tensorinfo.S_half_inv,
			this->fe_data_structure.dof_handler,
			this->fe_data_structure.mapping,
			this->nEqn);		
		
		timer.leave_subsection();

		// Grid refinement should be done in the end.
		Mesh_Info.refinement_handling(cycle,this->constants.refine_cycles);
		
	}
}



}
