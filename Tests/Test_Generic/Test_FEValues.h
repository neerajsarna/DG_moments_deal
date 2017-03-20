namespace Test_FEValues
{
	using namespace dealii;

	// construct the block strucutre for the allocation of fe nothing
	void construct_block_structure(std::vector<int> &block_structure,std::vector<int> &nEqn)
	{
		Assert(nEqn.size() != 0, ExcNotInitialized());
		Assert(std::is_sorted(std::begin(nEqn),std::end(nEqn)),ExcMessage("number of equations not sorted"));

		// the very first entry should be the number of equations in the first system
		block_structure.push_back(nEqn[0]);

		for (int i = 1 ; i < nEqn.size() ; i++)
		{
			block_structure.push_back(nEqn[i]-nEqn[i-1]);
			std::cout << "Block structure " << block_structure[i] << std::endl;
		}


		AssertDimension(block_structure.size(),nEqn.size());
	}




	TEST(DofIndex,HandlesDofIndex)
	{
		const unsigned int dim = 2;
		const unsigned int p = 2;
		const unsigned int nEqn = 6;


		FESystem<dim> finite_element(FE_DGQ<dim>(p),3,FE_DGQ<dim>(p),3);

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

	TEST(QTrapez,HandlesQTrapez)
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


	

		QTrapez<dim> quadrature;
		
	
		FEValues<dim> fe_values(mapping,finite_element,quadrature,update_flags);
		

		typename DoFHandler<dim>::active_cell_iterator cell=dof_handler.begin_active(),
													 	endc = dof_handler.end();

		fe_values.reinit(cell);

        std::vector<Point<dim>> quadrature_points = fe_values.get_quadrature_points();

        for (unsigned int q = 0 ; q < quadrature.size() ; q++)
        	for (unsigned long int i = 0 ; i < finite_element.dofs_per_cell ; i ++)
        		std::cout << "quad points " << q << " Dof " << i << " shape_value " << fe_values.shape_value(i,q)<< std::endl;		
	}

	TEST(MassMatrix,HandlesMassMatrix)
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




			// for (unsigned int face = 0 ; face < GeometryInfo<dim>::faces_per_cell; face++)
			// {
			// 	fe_values_face.reinit(cell,face);
			// 	J_face = fe_values_face.get_JxW_values();
			// 	Mass_face_manuel = 0;

			// 	for (unsigned int i = 0 ; i < num_dof ; i++)
			// 		for (unsigned int j = 0 ; j < num_dof ; j++)
			// 			for (long unsigned int q = 0 ; q < num_gauss_points_face ; q++)
			// 				Mass_face_manuel(i,j) += fe_values_face.shape_value(i,q) * fe_values_face.shape_value(j,q)
			// 								  		 * J_face[q];



			// 	// Since we are considering a one dimensional equation therefore the number of dofs_per_component
			// 	// will be equal to the total number of degrees of freedom being considered 
			// 	Mass_face = numerical_integration.Compute_Mass_shape_value_face(fe_values_face,num_dof,J_face);			  	

			// 	for (unsigned int i = 0 ; i < num_dof ; i++)
			// 		for (unsigned int j = 0 ; j < num_dof ; j ++)	
			// 		 	EXPECT_NEAR(Mass_face_manuel(i,j),Mass_face(i,j),1e-5);
				
			//   }


			

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

	// in the following test, we understand the working of FE nothing
	TEST(FENothing,HandlesFENothing)
	{
		const int dim = 2;
		Triangulation<dim> triangulation;

		GridGenerator::hyper_cube(triangulation,-0.5,0.5);
		hp::FECollection<dim> fe_collection;
		hp::DoFHandler<dim> dof_handler(triangulation);
		
		FESystem<dim> fe_sys1(FE_DGQ<dim>(1),3,
							  FE_Nothing<dim>(),3,
							  FE_Nothing<dim>(),4);

		FESystem<dim> fe_sys2(FE_DGQ<dim>(1),3,
							  FE_DGQ<dim>(1),3,
							  FE_Nothing<dim>(),4);

		FESystem<dim> fe_sys3(FE_DGQ<dim>(1),3,
							  FE_DGQ<dim>(1),3,
							  FE_DGQ<dim>(1),4);


		fe_collection.push_back(fe_sys1);
		fe_collection.push_back(fe_sys2);
		fe_collection.push_back(fe_sys3);

		typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();

		for (; cell != endc ; cell++)
		{
			if (cell->center()(0) > 0.0 && cell->center()(0) < 0.1)
				cell->set_active_fe_index(0);

			if (cell->center()(0) >= 0.1)
				cell->set_active_fe_index(1);

			if (cell->center()(0) <= 0.0)
				cell->set_active_fe_index(2);
		}


		dof_handler.distribute_dofs(fe_collection);

		std::cout << "dofs sys1 " << fe_sys1.dofs_per_cell << std::endl;
		std::cout << "dofs sys2 " << fe_sys2.dofs_per_cell << std::endl;
		std::cout << "dofs sys3 " << fe_sys3.dofs_per_cell << std::endl;

		std::cout << "components sys1 " << fe_sys1.n_components() << std::endl;
		std::cout << "components sys2 " << fe_sys2.n_components() << std::endl;
		std::cout << "components sys3 " << fe_sys3.n_components() << std::endl;



	}

	TEST(BlockStructure,HandlesBlockStructure)
	{
		const int dim = 2;
		std::vector<int> nEqn;
		std::vector<int> block_structure;
		hp::FECollection<dim> finite_element;

		nEqn.push_back(3);
		nEqn.push_back(6);
		nEqn.push_back(10);
		nEqn.push_back(15);

		construct_block_structure(block_structure,nEqn);

		std::vector<FESystem<dim>> fe_sys;

		if (block_structure.size() == 1)
			fe_sys.push_back(FESystem<dim>(FE_DGQ<dim>(1),block_structure[0]));

		if (block_structure.size() == 2)
		{

			fe_sys.push_back(FESystem<dim>(FE_DGQ<dim>(1),block_structure[0],FE_Nothing<dim>(),block_structure[1]));
			fe_sys.push_back(FESystem<dim>(FE_DGQ<dim>(1),block_structure[0],FE_DGQ<dim>(1),block_structure[1]));
		}

		if (block_structure.size() == 3)
		{
			fe_sys.push_back(FESystem<dim>(FE_DGQ<dim>(1),block_structure[0],FE_Nothing<dim>(),block_structure[1],FE_Nothing<dim>(),block_structure[2]));
			fe_sys.push_back(FESystem<dim>(FE_DGQ<dim>(1),block_structure[0],FE_DGQ<dim>(1),block_structure[1],FE_Nothing<dim>(),block_structure[2]));
			fe_sys.push_back(FESystem<dim>(FE_DGQ<dim>(1),block_structure[0],FE_DGQ<dim>(1),block_structure[1],FE_DGQ<dim>(1),block_structure[2]));
		}

		if (block_structure.size() == 4)
		{
			fe_sys.push_back(FESystem<dim>(FE_DGQ<dim>(1),block_structure[0],FE_Nothing<dim>(),block_structure[1],
							 FE_Nothing<dim>(),block_structure[2],FE_Nothing<dim>(),block_structure[3]));

			fe_sys.push_back(FESystem<dim>(FE_DGQ<dim>(1),block_structure[0],FE_DGQ<dim>(1),block_structure[1],
							 FE_Nothing<dim>(),block_structure[2],FE_Nothing<dim>(),block_structure[3]));

			fe_sys.push_back(FESystem<dim>(FE_DGQ<dim>(1),block_structure[0],FE_DGQ<dim>(1),block_structure[1],
							 FE_DGQ<dim>(1),block_structure[2],FE_Nothing<dim>(),block_structure[3]));

			fe_sys.push_back(FESystem<dim>(FE_DGQ<dim>(1),block_structure[0],FE_DGQ<dim>(1),block_structure[1],
							 FE_DGQ<dim>(1),block_structure[2],FE_DGQ<dim>(1),block_structure[3]));
		}

		for (int i = 0 ; i < fe_sys.size() ; i++)
		{
			finite_element.push_back(fe_sys[i]);
			std::cout << "Dofs: " << finite_element[i].dofs_per_cell << std::endl;
		}
	}
}
