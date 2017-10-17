namespace FEM_Solver
{
	using namespace dealii;

	template<int dim>
	class
	Run_Problem_hp_Periodic:public Assembly_Manager_hp_Periodic<dim>,
						 	public Run_Problem<dim>

	{
		public:
			Run_Problem_hp_Periodic(const std::string &output_file_name,
						const constant_numerics &constants,
						std::vector<Develop_System::System<dim>> &equation_info,
                        ExactSolution::Base_ExactSolution<dim> *exact_solution,
                        const std::vector<int> &nEqn,
                        const std::vector<int> &nBC,
                        Triangulation<dim> &triangulation);


			// we take the reference solution from a higher order moment method
			void run_higher_order_reference(DoFHandler<dim> &dof_handler_reference,
											const Vector<double> &solution_reference,
											MeshGenerator::Base_MeshGenerator<dim> &Mesh_Info);

			void distribute_dof_allocate_matrix_periodic_box();

	};

	template<int dim>
	Run_Problem_hp_Periodic<dim>::Run_Problem_hp_Periodic(const std::string &output_file_name,
												    	const constant_numerics &constants,
														std::vector<Develop_System::System<dim>> &equation_info,
                        							   ExactSolution::Base_ExactSolution<dim> *exact_solution,
                        							   const std::vector<int> &nEqn,
                        							   const std::vector<int> &nBC,
                        							   Triangulation<dim> &triangulation)
	:
	Assembly_Manager_hp_Periodic<dim>(output_file_name,
									 constants,equation_info,
                					nEqn,nBC,triangulation),
	Run_Problem<dim>(exact_solution)
	{}

	template<int dim>
	void 
	Run_Problem_hp_Periodic<dim>::distribute_dof_allocate_matrix_periodic_box()
	{
		this->hp_fe_data_structure.dof_handler.distribute_dofs(this->hp_fe_data_structure.finite_element);

		DynamicSparsityPattern dsp(this->hp_fe_data_structure.dof_handler.n_dofs(),
								  this->hp_fe_data_structure.dof_handler.n_dofs());


		DoFTools::make_flux_sparsity_pattern (this->hp_fe_data_structure.dof_handler, dsp);

		this->hp_add_periodic_sparsity(dsp);
		this->global_matrix.reinit(dsp);  
	}

	template<int dim>
	void 
	Run_Problem_hp_Periodic<dim>::run_higher_order_reference(DoFHandler<dim> &dof_handler_reference,
														  const Vector<double> &solution_reference,
									     				  MeshGenerator::Base_MeshGenerator<dim> &Mesh_Info)
	{
			const int refine_cycles_h = this->constants.refine_cycles;

			// total number of refinement cycles in the velocity space
			const int refine_cycles_c = this->constants.refine_cycles_c;

			this->error_per_itr.resize(refine_cycles_c * refine_cycles_h);

			TimerOutput timer (std::cout, TimerOutput::summary,
                	   TimerOutput::wall_times);

			Mesh_Info.print_mesh_info();

			// we first do h-refinement
			for (int cycle_h = 0 ; cycle_h < refine_cycles_h ; cycle_h ++)
			{
				for (int cycle_c = 0 ; cycle_c < refine_cycles_c ; cycle_c ++)
				{
					std::cout << "Velocity Space Refinement Cycle " << cycle_c << std::endl;

					// allocate the index for every cell
					timer.enter_subsection("Dof Distribution");
					std::cout << "Dof distirubtion" << std::endl;
					this->hp_fe_data_structure.allocate_fe_index_error_comparison(cycle_c);
					timer.leave_subsection();

					// // First we develop the periodic faces using internal functions of dealii.
					this->hp_develop_periodic_faces(this->hp_fe_data_structure.dof_handler);

					// // // now we construct the required data structure
					this->hp_divide_periodicity();

					// // distribute the degrees of freedom for the different fe indices which have been distributed
					this->distribute_dof_allocate_matrix_periodic_box();

					this->allocate_vectors(this->hp_fe_data_structure.dof_handler,this->solution,this->system_rhs,
										   this->residual);

					
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


					if (cycle_c != refine_cycles_c)
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


					postproc.print_options_quad_points(Mesh_Info.triangulation,
													   this->solution,
													   cycle_c,
														refine_cycles_c,
														this->convergence_table,
										  this->system_info[this->hp_fe_data_structure.max_fe_index].base_tensorinfo.S_half_inv,
										  this->hp_fe_data_structure.dof_handler,
										  this->hp_fe_data_structure.mapping,
										  this->nEqn);


				  timer.leave_subsection();

				}

			}

	}



}

	

