namespace Test_FEValues
{
	using namespace dealii;

	// we would like to know whether all the cells have the same fe_values and the fe_shapevalues
	TEST(FEValues,HandlesFEValues)
	{
                const unsigned int dim = 2;
		const unsigned int p = 2;
                ASSERT_EQ(dim,2) << "3D not implemented" << std::endl;

		SphericalManifold<dim> boundary;
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

       		triangulation.set_all_manifold_ids_on_boundary(0);
        	triangulation.set_manifold(0,boundary);

		
		FE_DGQ<dim> finite_element(p);
		DoFHandler<dim> dof_handler(triangulation);
		
		dof_handler.distribute_dofs(finite_element);
		MappingQ<dim> mapping(p+1);

		UpdateFlags update_flags = update_quadrature_points |
                	           update_values            |
                        	   update_gradients;
	

		QGauss<dim> quadrature(p+1);
	
		FEValues<dim> fe_values(mapping,finite_element,quadrature,update_flags);

		typename Triangulation<dim>::active_cell_iterator cell=triangulation.begin_active(), endc = triangulation.end();

		const unsigned int num_dof = finite_element.dofs_per_cell;
		const unsigned int num_gauss_points = quadrature.size();

		FullMatrix<double> shape_value(num_dof,num_gauss_points);
		FullMatrix<double> shape_grad_x(num_dof,num_gauss_points);
		FullMatrix<double> shape_grad_y(num_dof,num_gauss_points);
		

		fe_values.reinit(cell);

		for (unsigned int i = 0 ; i < num_dof ; i++)
			for (unsigned int q = 0 ; q < num_gauss_points ; q++)
				shape_value(i,q) = fe_values.shape_value(i,q);
		


		for(; cell != endc ; cell++)
		{
			fe_values.reinit(cell);

			for (unsigned int i = 0 ; i < num_dof ; i++)
				for (unsigned int q = 0 ; q < num_gauss_points ; q++)
					EXPECT_NEAR(shape_value(i,q),fe_values.shape_value(i,q),1e-5) << "Cell Index " << cell->index() << std::endl;
					

		}
		
	}

	TEST(DofIndex,HandlesDofIndex)
	{
		const unsigned int dim = 2;
		const unsigned int p = 2;
		const unsigned int nEqn = 6;


		FESystem<dim> finite_element(FE_DGQ<dim>(p),nEqn);

		// or the indices per cell
		const unsigned int dofs_per_component = finite_element.dofs_per_cell / finite_element.n_components();
		FullMatrix<double> component_to_system(nEqn,dofs_per_component);

		unsigned int count = 0;

		// loop over a particular component
		for(unsigned int i = 0 ; i < nEqn ; i++)

		// loop over all the dofs of that particular component
			for (unsigned int j = 0 ; j < dofs_per_component ; j++)
			{
				// compute the dof to which the it corresponds to 
				component_to_system(i,j) = finite_element.component_to_system_index(i,j);
				EXPECT_EQ(count,component_to_system(i,j));
				count ++;
			}

	}

	TEST(MassMatrix,HandlesMassMatrix)
	{
        const unsigned int dim = 2;
		const unsigned int p = 2;
        ASSERT_EQ(dim,2) << "3D not implemented" << std::endl;

		SphericalManifold<dim> boundary;
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

       		triangulation.set_all_manifold_ids_on_boundary(0);
        	triangulation.set_manifold(0,boundary);

		
		FE_DGQ<dim> finite_element(p);
		DoFHandler<dim> dof_handler(triangulation);
		
		dof_handler.distribute_dofs(finite_element);
		MappingQ<dim> mapping(p+1);

		UpdateFlags update_flags = update_quadrature_points |
                	           		update_values            |
                        	   		update_gradients		  |
                        	   		update_JxW_values;
	

		QGauss<dim> quadrature(p+1);
	
		FEValues<dim> fe_values(mapping,finite_element,quadrature,update_flags);

		typename DoFHandler<dim>::active_cell_iterator cell=dof_handler.begin_active(),
													 	endc = dof_handler.end();

		const unsigned int num_dof = finite_element.dofs_per_cell;
		const unsigned int num_gauss_points = quadrature.size();

		std::vector<double> J(num_gauss_points);
		FullMatrix<double> Mass_x_manuel(num_dof,num_dof);
		Full_matrix Mass_x(num_dof,num_dof);
		
		FullMatrix<double> Mass_y_manuel(num_dof,num_dof);
		Full_matrix Mass_y(num_dof,num_dof);
		
		FullMatrix<double> Mass_z_manuel(num_dof,num_dof);
		Full_matrix Mass_z(num_dof,num_dof);

		FullMatrix<double> Mass_manuel(num_dof,num_dof);
		Full_matrix Mass(num_dof,num_dof);

		NumericalIntegration::Base_NumericalIntegration<dim> numerical_integration;

		// we first need to compute the shape value so that the shape_values in the class numerical_integration 
		// can be initialized 
		numerical_integration.Compute_Shape_Value(mapping,
												  p,
												  p + 1,
												  cell);

		for (; cell != endc ; cell++)
		{
			fe_values.reinit(cell);

			Mass_x_manuel = 0;
			Mass_y_manuel = 0;
			Mass_z_manuel = 0;
			Mass_manuel = 0;

			J = fe_values.get_JxW_values();

			// we now prepare the mass matrix ourselves and compare  it with that returned by the function
			for (unsigned int i = 0 ; i < num_dof ; i++)
				for (unsigned int j = 0 ; j < num_dof ; j++)
					for (unsigned int q = 0 ;q < num_gauss_points ; q++)
						Mass_x_manuel(i,j) += fe_values.shape_value(i,q ) * fe_values.shape_grad(j,q)[0]
											  * J[q];


			Mass_x = numerical_integration.Compute_Mass_shape_grad_x(fe_values, num_dof, J);


			for (unsigned int i = 0 ; i < num_dof ; i++)
				for (unsigned int j = 0 ; j < num_dof ; j++)
					for (unsigned int q = 0 ;q < num_gauss_points ; q++)
						Mass_y_manuel(i,j) += fe_values.shape_value(i,q) * fe_values.shape_grad(j,q)[1]
											  * J[q];


			Mass_y = numerical_integration.Compute_Mass_shape_grad_y(fe_values, num_dof, J);

			for (unsigned int i = 0 ; i < num_dof ; i++)
				for (unsigned int j = 0 ; j < num_dof ; j++)
					for (unsigned int q = 0 ;q < num_gauss_points ; q++)
						Mass_manuel(i,j) += fe_values.shape_value(i,q) * fe_values.shape_value(j,q)
											  * J[q];


			Mass = numerical_integration.Compute_Mass_shape_value(fe_values, num_dof, J);


			for (unsigned int i = 0 ; i < num_dof ; i++)
				for (unsigned int j = 0 ; j < num_dof ; j ++)
				{
					EXPECT_NEAR(Mass_x_manuel(i,j),Mass_x(i,j),1e-5);
					EXPECT_NEAR(Mass_y_manuel(i,j),Mass_y(i,j),1e-5);
					EXPECT_NEAR(Mass_manuel(i,j),Mass(i,j),1e-5);
				}
		}
		
	}

	// we now test our technique to come up with a faster assemblation
	TEST(FasterAssemblation,TestFasterAssemblation)
	{
        const unsigned int dim = 2;
		const unsigned int p = 2;
        ASSERT_EQ(dim,2) << "3D not implemented" << std::endl;

		SphericalManifold<dim> boundary;
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

       		triangulation.set_all_manifold_ids_on_boundary(0);
        	triangulation.set_manifold(0,boundary);

		// consider a model problem with 2 equation
        const unsigned int nEqn = 2;
		FESystem<dim> finite_element(FE_DGQ<dim>(p),nEqn);
		Sparse_matrix A(nEqn,nEqn);

		DoFHandler<dim> dof_handler(triangulation);
		
		dof_handler.distribute_dofs(finite_element);
		MappingQ<dim> mapping(p+1);

		UpdateFlags update_flags = update_quadrature_points |
                	           		update_values            |
                        	   		update_gradients		  |
                        	   		update_JxW_values;
	

		QGauss<dim> quadrature(p+1);
	
		FEValues<dim> fe_values(mapping,finite_element,quadrature,update_flags);

		typename DoFHandler<dim>::active_cell_iterator cell=dof_handler.begin_active(),
													 	endc = dof_handler.end();

		const unsigned int num_dof = finite_element.dofs_per_cell;
		const unsigned int num_components = finite_element.n_components();
		const unsigned int num_dof_per_comp = num_dof/num_components;
		const unsigned int num_gauss_points = quadrature.size();

		std::vector<double> J(num_gauss_points);
		std::vector<std::vector<double>> component_to_system(num_components,std::vector<double> (num_dof_per_comp));

		Full_matrix Mass_x(num_dof_per_comp,num_dof_per_comp);
		Full_matrix Mass_y(num_dof_per_comp,num_dof_per_comp);
		Sparse_matrix cell_matrix(num_dof,num_dof);
		Full_matrix cell_matrix_manuel(num_dof,num_dof);

		NumericalIntegration::Base_NumericalIntegration<dim> numerical_integration;
		MatrixOpt::Base_MatrixOpt matrix_opt;

		// we first need to compute the shape value so that the shape_values in the class numerical_integration 
		// can be initialized 
		numerical_integration.Compute_Shape_Value(mapping,
												  p,
												  p + 1,
												  cell);

		A.coeffRef(0,0) = 1.0;
		A.coeffRef(0,1) = 2.0;
		A.coeffRef(1,1) = 2.5;

		for (unsigned int i = 0 ; i < num_components ; i ++)
    		for (unsigned int j = 0 ; j < num_dof_per_comp ; j ++)
      			component_to_system[i][j] = finite_element.component_to_system_index(i,j); 

		A.makeCompressed();

		for (; cell != endc ; cell++)
		{
			fe_values.reinit(cell);
			cell_matrix_manuel.setZero();

			J = fe_values.get_JxW_values();
			Mass_x = numerical_integration.Compute_Mass_shape_grad_x(fe_values, num_dof_per_comp, J);
			Mass_y = numerical_integration.Compute_Mass_shape_grad_y(fe_values, num_dof_per_comp, J);

			// first the faster computation
			// The reason it is faster is that we not assemble or compute the component of Mass_x for every equation and every 
			// component independently.
			cell_matrix = matrix_opt.compute_A_outer_B(A,Mass_x) +  matrix_opt.compute_A_outer_B(A,Mass_y);

			for (unsigned int q = 0 ; q < num_gauss_points ; q++)
				for (unsigned int index_test = 0 ; index_test < num_dof_per_comp ; index_test++)
					for (unsigned int index_sol = 0 ; index_sol < num_dof_per_comp ; index_sol++ )
						for (unsigned int space = 0 ; space < dim ; space++)
							for (unsigned int m = 0 ; m < A.outerSize(); m++)
							{
								const int dof_test = component_to_system[m][index_test];
          //const double shape_value_test = fe_v.shape_value(dof_test,q);
								const double shape_value_test = numerical_integration.shape_values(index_test,q);

								for (Sparse_matrix::InnerIterator n(A,m); n ; ++n)
								{
									int dof_sol = component_to_system[n.col()][index_sol];
									double grad_value_sol = fe_values.shape_grad(dof_sol,q)[space];

									cell_matrix_manuel(dof_test,dof_sol) += shape_value_test * n.value()
																		* grad_value_sol * J[q];  
								}
							}

				for (unsigned int i = 0 ; i < num_dof ; i++)
					for (unsigned int j = 0 ; j < num_dof ; j++)
						EXPECT_NEAR(cell_matrix.coeffRef(i,j),cell_matrix_manuel.coeffRef(i,j),1e-5);
			}
	}

	
}
