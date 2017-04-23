namespace Test_FEValues
{
	using namespace dealii;


	template<int dim>
	class
	initial_value:public Function<dim>
	{
		public:
			initial_value(const int nEqn);


		virtual void vector_value(const Point<dim> &p,Vector<double> &value) const ;

	};

	template<int dim>
	initial_value<dim>::initial_value(const int nEqn)
	:
	Function<dim>(nEqn)
	{
	}

	template<int dim>
	void 
	initial_value<dim>::vector_value(const Point<dim> &p,Vector<double> &value)const 
	{
		const double x = p(0);
		const double y = p(1);

		value = 0;

		// rho
		value(0) = pow(x,2);
		//velocity
		value(1) = 0;

		value(2) = 0;

		// theta
		value(3) = -sqrt(3.0/2.0) * exp(x) * exp(y);

		// sigmaxx
		value(4) = (sin(x) + sin(y))/sqrt(2);

		// sigmaxy
		value(5) = cos(x) * sin(y)/sqrt(2);

		//sigma yy
		value(6) = (pow(x,3) * 0.1 + pow(x,2) + pow(y,2))/sqrt(2);


	}

	TEST(DISABLED_LiftDrag,HandlesLiftDragComputation)
	{
		const unsigned int dim = 2;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);


	// 	// we construct a vector of all the system data being considered 
		std::vector<Develop_System::System<dim>> System;

		// initialize the vector containing all the systems
		// initialize the vector containing all the systems
		for(int i = 0 ; i < constants.constants_sys.total_systems ; i++)
		{
			System.push_back(Develop_System::System<dim>(constants.constants_num,constants.constants_sys.nEqn[i],
				constants.constants_sys.nBC[i],constants.constants_sys.Ntensors[i],folder_name));

		}
	
		for (int i = 0 ; i < constants.constants_sys.total_systems ; i++)
			System[i].initialize_system();


		// the exact solution can only be created for one of the systems
		ExactSolution::ExactSolution_Dummy<dim>  exactsolution_dummy(constants.constants_num,System[constants.constants_sys.total_systems-1].base_tensorinfo.S_half,
													constants.constants_sys.nEqn[constants.constants_sys.total_systems-1],constants.constants_sys.Ntensors[constants.constants_sys.total_systems-1]);


		// finite element solver for a single system
		FEM_Solver::Run_Problem_FE<dim> fe_solver("grid",
											constants.constants_num,
											System,
											&exactsolution_dummy,
											constants.constants_sys.nEqn,
											constants.constants_sys.nBC);

		fe_solver.distribute_dof_allocate_matrix(fe_solver.fe_data_structure.dof_handler,fe_solver.fe_data_structure.finite_element,
												fe_solver.global_matrix);



		fe_solver.allocate_vectors(fe_solver.fe_data_structure.dof_handler,fe_solver.solution,
									fe_solver.system_rhs,fe_solver.residual);

		initial_value<dim> solution_function(System[0].nEqn);

		VectorTools::interpolate(fe_solver.fe_data_structure.mapping,
								 fe_solver.fe_data_structure.dof_handler,
								 solution_function,
								 fe_solver.solution);


		PostProc::Base_PostProc<dim> postproc(fe_solver.constants,fe_solver.base_exactsolution,
											 	fe_solver.nEqn,fe_solver.nBC);


		 postproc.reinit(fe_solver.fe_data_structure.dof_handler);

		// // we now compute the lift and drag over a surface of a plate
		 Vector<double> lift_drag = postproc.compute_lift_drag(fe_solver.fe_data_structure.mapping,
    								fe_solver.fe_data_structure.finite_element,
    								fe_solver.fe_data_structure.dof_handler,
    								System[0].base_tensorinfo.S_half_inv,
    								fe_solver.solution,
    								fe_solver.convergence_table,
    								3);


		 std::cout << lift_drag << std::endl;
	}



	TEST(DISABLED_MassMatrix,HandlesMassMatrix)
	{
        const unsigned int dim = 2;
		const unsigned int p = 1;
        ASSERT_EQ(dim,2) << "3D not implemented" << std::endl;

		SphericalManifold<dim> boundary;
		Triangulation<dim> triangulation;
	        
		Point<dim> center;
		const double inner_radius = 0.5;
		const double outer_radius = 2.0;
		const unsigned int parts = 5;

		center(0) = 0.0;
		center(1) = 0.0;

		GridGenerator::hyper_shell (triangulation,
			center, inner_radius, outer_radius,
			parts);

		triangulation.set_all_manifold_ids_on_boundary(0);
		triangulation.set_manifold(0,boundary);

		std::cout << "#Cells " << triangulation.n_active_cells() << std::endl;
		
		FE_DGQ<dim> finite_element(p);
		DoFHandler<dim> dof_handler(triangulation);
		
		dof_handler.distribute_dofs(finite_element);
		MappingQ<dim> mapping(p+1);

		UpdateFlags update_flags = update_quadrature_points |
                	           		update_values            |
                        	   		update_gradients		  |
                        	   		update_JxW_values,
                face_update_flags = update_values
      								| update_q_points
      								| update_JxW_values
      								| update_normal_vectors;


	

		QGauss<dim> quadrature(p+1);
		QGauss<dim-1> quadrature_face(p+1);
	
		FEValues<dim> fe_values(mapping,finite_element,quadrature,update_flags);
		FEFaceValues<dim> fe_values_face(mapping,finite_element,quadrature_face,face_update_flags);

		typename DoFHandler<dim>::active_cell_iterator cell=dof_handler.begin_active(),
													 	endc = dof_handler.end();

		const unsigned int num_dof = finite_element.dofs_per_cell;
		const unsigned int num_gauss_points = quadrature.size();
		const unsigned int num_gauss_points_face = quadrature_face.size();

		std::vector<double> J(num_gauss_points);
		std::vector<double> J_face(num_gauss_points);

		
		std::vector<FullMatrix<double>> Mass_grad;
		FullMatrix<double> Mass_x_manuel(num_dof,num_dof);
		FullMatrix<double> Mass_y_manuel(num_dof,num_dof);
		FullMatrix<double> Mass_z_manuel(num_dof,num_dof);
		

		FullMatrix<double> Mass_manuel(num_dof,num_dof);
		FullMatrix<double> Mass(num_dof,num_dof);

		FullMatrix<double> Mass_face_manuel(num_dof,num_dof);
		FullMatrix<double> Mass_face(num_dof,num_dof);

		NumericalIntegration::Base_NumericalIntegration<dim> numerical_integration;

		// we first need to compute the shape value so that the shape_values in the class numerical_integration 
		// can be initialized 
		numerical_integration.Compute_Shape_Value(mapping,
												  p+1,
												  cell);


		for (; cell != endc ; cell++)
		{
			fe_values.reinit(cell);

			Mass_x_manuel = 0;
			Mass_y_manuel = 0;
			Mass_z_manuel = 0;
			Mass_manuel = 0;
			Mass_face_manuel = 0;

			J = fe_values.get_JxW_values();

			// we now prepare the mass matrix ourselves and compare  it with that returned by the function
			for (unsigned int i = 0 ; i < num_dof ; i++)
				for (unsigned int j = 0 ; j < num_dof ; j++)
					for (unsigned int q = 0 ;q < num_gauss_points ; q++)
						Mass_x_manuel(i,j) += fe_values.shape_value(i,q ) * fe_values.shape_grad(j,q)[0]
											  * J[q];

			Mass_grad = numerical_integration.Compute_Mass_shape_grad(fe_values, num_dof, J);


			for (unsigned int i = 0 ; i < num_dof ; i++)
				for (unsigned int j = 0 ; j < num_dof ; j++)
					for (unsigned int q = 0 ;q < num_gauss_points ; q++)
						Mass_y_manuel(i,j) += fe_values.shape_value(i,q) * fe_values.shape_grad(j,q)[1]
											  * J[q];



			for (unsigned int i = 0 ; i < num_dof ; i++)
				for (unsigned int j = 0 ; j < num_dof ; j++)
					for (unsigned int q = 0 ;q < num_gauss_points ; q++)
						Mass_manuel(i,j) += fe_values.shape_value(i,q) * fe_values.shape_value(j,q)
											  * J[q];


			Mass = numerical_integration.Compute_Mass_shape_value(fe_values,num_dof,J);

			// now we compute the Mass matrix for the boundary faces

			for (unsigned int i = 0 ; i < num_dof ; i++)
				for (unsigned int j = 0 ; j < num_dof ; j ++)
				{
					 EXPECT_NEAR(Mass_x_manuel(i,j),Mass_grad[0](i,j),1e-5);
					 EXPECT_NEAR(Mass_y_manuel(i,j),Mass_grad[1](i,j),1e-5);
					 EXPECT_NEAR(Mass_manuel(i,j),Mass(i,j),1e-5);					
				}				
		}
	
	
		
	}	



 }
