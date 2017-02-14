// in the following file we test our residual computation
namespace TestPostProc
{
	using namespace dealii;

	template<int dim> class dummy_function:public Function<dim>
	{
		public:
			dummy_function(const unsigned int nEqn):Function<dim>(nEqn),
													nEqn(nEqn){};

			virtual void vector_value(const Point<dim> &p,Vector<double> &value) const ;
			const unsigned int nEqn;
	};

	template<int dim> void dummy_function<dim>::vector_value(const Point<dim> &p,Vector<double> &value) const
	{
		
  		const double x = p[0];
  		const double y = p[1];
  		

  		Assert(dim == 2,ExcNotImplemented());

  		for(unsigned int eq = 0 ; eq < nEqn ; eq++)
  			value(eq) = sin(M_PI *x ) * sin(M_PI*y); 
  		
	}

	TEST(SolvingSystemA,HandlesSolvingSystemA)
	{
		const unsigned int dim = 2;
		ASSERT_EQ(dim,2) << "3D not implemented" << std::endl;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		SystemA::SystemA<dim> systemA(constants.constants,folder_name);

		ExactSolution::ExactSolution_SystemA_ring<dim>  exact_solution_systemA(constants.constants,systemA.base_tensorinfo.S_half);

		Triangulation<dim> triangulation;

		Point<dim> center;
		const double inner_radius = 0.5;
		const double outer_radius = 2.0;
		const unsigned int parts = 10;

		center(0) = 0.0;
		center(1) = 0.0;

		GridGenerator::hyper_shell (triangulation,
        							center, inner_radius, outer_radius,
              						parts);

		triangulation.refine_global(5);

			// triangulation.clear();
   //          Point<dim> p1;
   //          Point<dim> p2;
   //          std::vector<unsigned int > repetitions(dim);

   //          p1(0) = 0.5;
   //          p1(1) = 0.5;

   //          p2(0) = 1.0;
   //          p2(1) = 1.0;

   //          repetitions[0] = 100;
   //          repetitions[1] = 100;

   //          //The diagonal of the rectangle is the time joining p1 and p2
   //          GridGenerator::subdivided_hyper_rectangle(triangulation,
   //                                    	              repetitions,
   //                                      	            p1,
   //                                          	        p2);


			FESystem<dim> finite_element(FE_Q<dim>(1),constants.constants.nEqn);
			DoFHandler<dim> dof_handler(triangulation);
			MappingQ<dim,dim> mapping(1);
			dof_handler.distribute_dofs(finite_element);
			Vector<double> solution(dof_handler.n_dofs());
			ConstraintMatrix constraints;


			std::cout << "number of active cells " << triangulation.n_active_cells() << std::endl;
			std::cout << "number of Dofs " << dof_handler.n_dofs() << std::endl;

			constraints.close();
			VectorTools::project(dof_handler,
								constraints,
								QGauss<dim>(3),
								exact_solution_systemA,
								solution);

			typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
															end_c = dof_handler.end();

			const unsigned int dofs_per_cell = finite_element.dofs_per_cell;
			std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
			Vector<double> value(constants.constants.nEqn);
			std::vector<Tensor<1,dim,double>> value_gradient(constants.constants.nEqn);
			

			PostProc::Base_PostProc<dim> postproc(constants.constants,&exact_solution_systemA,&dof_handler,&mapping);				
			double residual_strong_form = postproc.compute_residual(solution,finite_element,&systemA,triangulation.n_active_cells()); 
			std::cout << "residual: " << residual_strong_form << std::endl;

	}

}