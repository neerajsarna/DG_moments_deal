namespace FEM_Solver
{
	using namespace dealii;
// the following routines asemble for a periodic 2D grid
	template<int dim>
	class
	Assembly_Manager_Periodic:public NumericalIntegration::Base_NumericalIntegration<dim>,
                            public PeriodicityHandler::Base_PeriodicityHandler<dim>
	{
	public:
		typedef MeshWorker::DoFInfo<dim> DoFInfo;
		typedef MeshWorker::IntegrationInfo<dim> CellInfo;



        DeclException2 (ExcCellCenter, double, double,
                        << "Cell Center = " << arg1 << "Neighbor Center = " << arg2 << "Mesh not lexiographical ");

        DeclException2 (ExcNoElementInSparsity, size_t, size_t,
                        << "Dof-1 " << arg1 << " Dof-2 " << arg2 << " Entry does not exist in sparsity pattern");

		Assembly_Manager_Periodic(const std::string &output_file_name,
						            const constant_numerics &constants,
						            std::vector<Develop_System::System<dim>> &equation_info,
                        const std::vector<int> &nEqn,
                        const std::vector<int> &nBC,
                        const unsigned int system_to_solve,
                        Triangulation<dim> &triangulation);
		    
		     // finite element data structure 
     		
     	fe_data<dim> fe_data_structure;

      const unsigned int dofs_per_cell;
      const unsigned int dofs_per_component;

			const int nEqn;
			const int nBC;
			const constant_numerics constants;

			Develop_System::System<dim> system_info;

			// A matrix in the final Ax = b
			TrilinosWrappers::SparseMatrix global_matrix;
            

			// The rhs in final Ax = b.
        	Vector<double> system_rhs;



        	const unsigned int ngp;
        	const unsigned int ngp_face;
            unsigned int integrate_inflow = 0;



			void assemble_system_manuel();


			// routines for manuel assembly
       		void integrate_cell_manuel(FullMatrix<double> &cell_matrix,
                                        Vector<double> &cell_rhs,
                                        const FEValuesBase<dim> &fe_v,
                                        std::vector<double> &J,
                                        std::vector<Vector<double>> &source_term_value,
                                        std::vector<Vector<double>> &component_to_system);

       // // integrate the boundary manuelly using odd boundary implementation
       void integrate_boundary_manuel_odd(FullMatrix<double> &cell_matrix,
                                          Vector<double> &cell_rhs,
                                          const FEValuesBase<dim> &fe_v,
                                          std::vector<double> &J,
                                          const unsigned int b_id);

       void integrate_boundary_manuel_kinetic(FullMatrix<double> &cell_matrix,
                                          Vector<double> &cell_rhs,
                                          const FEValuesBase<dim> &fe_v,
                                          std::vector<double> &J,
                                          const unsigned int b_id);
        	
      void integrate_face_manuel(FullMatrix<double> &u1_v1,
                        FullMatrix<double> &u1_v2,
                        FullMatrix<double> &u2_v1,
                        FullMatrix<double> &u2_v2,
                        const FEValuesBase<dim> &fe_v,
                        const FEValuesBase<dim> &fe_v_neighbor,
                        std::vector<double> &J);


     MatrixOpt::Base_MatrixOpt matrix_opt;


	};

	template<int dim>
	Assembly_Manager_Periodic<dim>::Assembly_Manager_Periodic(const std::string &output_file_name,
											const constant_numerics &constants,
											std::vector<Develop_System::System<dim>> &equation_info,
                        					const std::vector<int> &nEqn,
                        					const std::vector<int> &nBC,          
                                  const unsigned int system_to_solve,
                                  Triangulation<dim> &triangulation)
	:
  PeriodicityHandler::Base_PeriodicityHandler<dim>(constants.xl,constants.xr),
	fe_data_structure(output_file_name,
					         constants,nEqn[system_to_solve],triangulation),
  dofs_per_cell(fe_data_structure.finite_element.dofs_per_cell),
  dofs_per_component(dofs_per_cell/nEqn[system_to_solve]),
  // we only save data corresponding to that system which is needed to be solved
	nEqn(nEqn[system_to_solve]),
	nBC(nBC[system_to_solve]),
	constants(constants),
	system_info(equation_info[system_to_solve]),
	ngp(constants.p + 1),
	ngp_face(constants.p + 1)
	{
  	}


	template<int dim>
	void 
	Assembly_Manager_Periodic<dim>::assemble_system_manuel()
	{
      const QGauss<dim> quadrature(ngp);
      const QGauss<dim-1> face_quadrature(ngp_face);


      const UpdateFlags update_flags               = update_values
      | update_gradients
      | update_q_points
      | update_JxW_values,
      face_update_flags          = update_values
      | update_q_points
      | update_JxW_values
      | update_normal_vectors,
      neighbor_face_update_flags = update_values;

      FEValues<dim>  fe_v(this->fe_data_structure.mapping,this->fe_data_structure.finite_element,
      						quadrature, update_flags);
      FEFaceValues<dim> fe_v_face(this->fe_data_structure.mapping,this->fe_data_structure.finite_element, 
      							  face_quadrature, face_update_flags);
      FEFaceValues<dim> fe_v_face_neighbor(this->fe_data_structure.mapping,this->fe_data_structure.finite_element, 
      									face_quadrature, neighbor_face_update_flags);
      FESubfaceValues<dim> fe_v_subface_neighbor(this->fe_data_structure.mapping,this->fe_data_structure.finite_element, 
      									face_quadrature, neighbor_face_update_flags);  


      // the total number of quadrature points are the same for b
      const unsigned int total_ngp = quadrature.size();
      const unsigned int total_ngp_face = face_quadrature.size();

      typename DoFHandler<dim>::active_cell_iterator cell = this->fe_data_structure.dof_handler.begin_active(),
      												 endc = this->fe_data_structure.dof_handler.end();

      typename DoFHandler<dim>::cell_iterator neighbor;


      // it is much more efficient to declare the memory and then use it again and again
      FullMatrix<double> cell_matrix(dofs_per_cell,dofs_per_cell);                  // contribution from the interior
      FullMatrix<double> boundary_matrix(dofs_per_cell,dofs_per_cell);
      Vector<double>   cell_rhs(dofs_per_cell);                                   // rhs from the current cell

      FullMatrix<double> u1_v1(dofs_per_cell,dofs_per_cell);      // relating current to current
      FullMatrix<double> u1_v2(dofs_per_cell,dofs_per_cell);  // relating current to neighbor
      FullMatrix<double> u2_v1(dofs_per_cell,dofs_per_cell);  // relating neighbor to current
      FullMatrix<double> u2_v2(dofs_per_cell,dofs_per_cell);  // relating neighbor to current

      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      std::vector<types::global_dof_index> local_dof_indices_neighbor(dofs_per_cell);

      // std::vectors to make computations faster
      std::vector<double> Jacobians_interior(total_ngp);
      std::vector<double> Jacobian_face(total_ngp_face);


      // source term
      std::vector<Vector<double>> source_term_value(total_ngp,Vector<double>(nEqn));

      this->Compute_Shape_Value(fe_data_structure.mapping,ngp,cell);
      
      std::vector<Vector<double>> component_to_system(nEqn,
                                                          Vector<double>(dofs_per_component));

      for (int i = 0 ; i < nEqn ; i ++)
        for (unsigned int j = 0 ; j < dofs_per_component ; j ++)
          component_to_system[i](j) = this->fe_data_structure.finite_element.component_to_system_index(i,j); 

  for (; cell != endc ; cell++) 
      {
        cell_rhs = 0;
        cell->get_dof_indices(local_dof_indices);

        fe_v.reinit(cell);
        Jacobians_interior = fe_v.get_JxW_values();
        system_info.source_term(fe_v.get_quadrature_points(),source_term_value);

        
        integrate_cell_manuel(cell_matrix,
                              cell_rhs,
                              fe_v,
                              Jacobians_interior,
                              source_term_value,
                              component_to_system);
     


            for(unsigned int face  = 0; face< GeometryInfo<dim>::faces_per_cell; face++ )
            {
            
              fe_v_face.reinit(cell,face);
              const typename DoFHandler<dim>::face_iterator face_itr = cell->face(face);
              

              Jacobian_face = fe_v_face.get_JxW_values();
              
              boundary_matrix = 0;

           
              if (face_itr->at_boundary())
              { 
             

                // id of the boundary
              	const unsigned int b_id = face_itr->boundary_id();

              	const double xcord_face_center = face_itr->center()(0);
              	const double ycord_cell_center = cell->center()(1);

              	if ( b_id == 100 || b_id == 101) 
              	{
              		neighbor = this->get_periodic_neighbor(xcord_face_center,ycord_cell_center);
              		neighbor->get_dof_indices(local_dof_indices_neighbor);

              		Assert(!neighbor->has_children(),
              			ExcMessage("periodic boundary only to be used with zero level meshes"));
              		Assert(fabs(neighbor->center()(1) - ycord_cell_center) < 1e-8,
              			ExcCellCenter(ycord_cell_center,neighbor->center()(1)));

              		const unsigned int neighbor_face = this->get_periodic_neighbor_face(xcord_face_center
              																			,ycord_cell_center);

              		fe_v_face_neighbor.reinit(neighbor,neighbor_face);


                      // integrate the face between the cell and the periodic neighbor
              		integrate_face_manuel(u1_v1,u1_v2,
              							  u2_v1,u2_v2,
              							  fe_v_face,
              							  fe_v_face_neighbor,
              							 Jacobian_face);  

              		for (unsigned int i = 0 ; i < dofs_per_cell ; i++)
              			for (unsigned int j = 0 ; j < dofs_per_cell ; j++)
              			{

              				global_matrix.add(local_dof_indices[i],local_dof_indices[j],u1_v1(i,j));
              				global_matrix.add(local_dof_indices_neighbor[i],local_dof_indices_neighbor[j],u2_v2(i,j));

              				global_matrix.add(local_dof_indices[i],local_dof_indices_neighbor[j],u2_v1(i,j));
              				global_matrix.add(local_dof_indices_neighbor[i],local_dof_indices[j],u1_v2(i,j)) ;
              			}                  


              		}
              	else
              	{
                integrate_boundary_manuel_kinetic(boundary_matrix, 
                                              cell_rhs,
                                              fe_v_face, 
                                              Jacobian_face,
                                              b_id); 

                      
                // assemble the terms on the boundary
                for (unsigned int i = 0 ; i < dofs_per_cell ; i++)
                  for (unsigned int j = 0 ; j < dofs_per_cell ; j++)
                   global_matrix.add(local_dof_indices[i],local_dof_indices[j],boundary_matrix(i,j));
           		}
                
              }
             
              
                // if the face is not at the boundary then we need to integrate over the face
                else
               {            
                
                 neighbor = cell->neighbor(face);

   
                 if (cell->neighbor_is_coarser(face))
                 {

                   Assert(!cell->has_children(), ExcInternalError());
                   Assert(!neighbor->has_children(), ExcInternalError());

                   neighbor->get_dof_indices(local_dof_indices_neighbor);

                   const std::pair<unsigned int, unsigned int> neighbor_face_no
                   											= cell->neighbor_of_coarser_neighbor(face);

                   fe_v_subface_neighbor.reinit(neighbor,neighbor_face_no.first,neighbor_face_no.second);

                   integrate_face_manuel(u1_v1,u1_v2,
                                         u2_v1,u2_v2,
                                         fe_v_face,
                                         fe_v_subface_neighbor,
                                         Jacobian_face);

                   for (unsigned int i = 0 ; i < dofs_per_cell ; i++)
                    for (unsigned int j = 0 ; j < dofs_per_cell ; j++)
                    {

                     global_matrix.add(local_dof_indices[i],local_dof_indices[j],u1_v1(i,j));
                     global_matrix.add(local_dof_indices_neighbor[i],local_dof_indices_neighbor[j],u2_v2(i,j));

                     global_matrix.add(local_dof_indices[i],local_dof_indices_neighbor[j],u2_v1(i,j));
                     global_matrix.add(local_dof_indices_neighbor[i],local_dof_indices[j],u1_v2(i,j)) ;
                   } 
                  }

                  
                 
                 else
                 {
                  if (neighbor->has_children())
                    continue;


                  Assert(cell->level() == neighbor->level(),ExcMessage("cell and the neighbor not on the same level"));

                  // compute the integral from only one side
                  if (neighbor < cell)
                    continue;

                  neighbor->get_dof_indices(local_dof_indices_neighbor);

                  const unsigned int neighbor_face_no = cell->neighbor_of_neighbor(face);
                  fe_v_face_neighbor.reinit(neighbor,neighbor_face_no);
                  
                  integrate_face_manuel(u1_v1,u1_v2,
                                        u2_v1,u2_v2,
                                        fe_v_face,fe_v_face_neighbor,
                                        Jacobian_face);  


                  for (unsigned int i = 0 ; i < dofs_per_cell ; i++)
                    for (unsigned int j = 0 ; j < dofs_per_cell ; j++)
                    {
                     
                
                     global_matrix.add(local_dof_indices[i],local_dof_indices[j],u1_v1(i,j));
                     global_matrix.add(local_dof_indices_neighbor[i],local_dof_indices_neighbor[j],u2_v2(i,j));

                     global_matrix.add(local_dof_indices[i],local_dof_indices_neighbor[j],u2_v1(i,j));
                     global_matrix.add(local_dof_indices_neighbor[i],local_dof_indices[j],u1_v2(i,j)) ;
                   }                   
  
                 }


               }
            }
              

            
              for (unsigned int i = 0 ; i < dofs_per_cell ; i++)
                  for (unsigned int j = 0 ; j < dofs_per_cell ; j++)          
                    global_matrix.add(local_dof_indices[i],local_dof_indices[j],cell_matrix(i,j));

              
             for(unsigned int i = 0 ; i < dofs_per_cell ; i++)
                  system_rhs(local_dof_indices[i]) += cell_rhs(i);

  
      }
	}


	template<int dim>
	void 
	Assembly_Manager_Periodic<dim>
	::integrate_cell_manuel(FullMatrix<double> &cell_matrix,
		Vector<double> &cell_rhs,
		const FEValuesBase<dim> &fe_v,
		std::vector<double> &J,
		std::vector<Vector<double>> &source_term_value,
    std::vector<Vector<double>> &component_to_system)
	{
    	
    	
    
    const FiniteElement<dim> &fe_in_cell = this->fe_data_structure.finite_element;
    const unsigned int total_ngp = J.size();
   
      std::vector<FullMatrix<double>> Mass_grad(dim);
      FullMatrix<double> Mass(dofs_per_component,dofs_per_component);

      Mass_grad = this->Compute_Mass_shape_grad(fe_v, dofs_per_component, J);
      Mass = this->Compute_Mass_shape_value(fe_v, dofs_per_component, J);

        // this is for initializing the matrix, it is necessary since we have not put cell_matrix to zero
      cell_matrix = this->matrix_opt.compute_A_outer_B(system_info.system_data.P.matrix,Mass);

      for (int space = 0 ;space < dim ; space ++ )
        cell_matrix.add(0,cell_matrix,1,
                    this->matrix_opt.compute_A_outer_B(system_info.system_data.A[space].matrix,Mass_grad[space]));

    for (unsigned int q = 0 ; q < total_ngp ; q++)
    {
    // value of the jacobian at the gauss point
      double jacobian_value = J[q];

    // index is the local id of the degree of freedom
      for (unsigned int index_test = 0 ; index_test < dofs_per_component ; index_test++)
      {
     // for the right hand side we iterate over all the components of the test function
        for (int m = 0 ; m < nEqn ; m++)
        {
          const int dof_test = component_to_system[m](index_test);
          const double shape_value_test = this->shape_values(index_test,q);

          cell_rhs(dof_test) += shape_value_test * source_term_value[q][m] * jacobian_value;
        }
      }

    }

	}



//Integrate the boundary term using manuel implementation
	template<int dim>
	void 
	Assembly_Manager_Periodic<dim>
	::integrate_boundary_manuel_odd(FullMatrix<double> &cell_matrix,
								Vector<double> &cell_rhs,
								const FEValuesBase<dim> &fe_v,
								std::vector<double> &J,
								const unsigned int b_id)
	{

		std::vector<unsigned int> component(dofs_per_cell);

		for (unsigned int i = 0 ; i < dofs_per_cell ; i ++)
			component[i] = this->fe_data_structure.finite_element.system_to_component_index(i).first;

		Vector<double> boundary_rhs_value;
		boundary_rhs_value.reinit(nBC);
		
		Assert(system_info.system_data.B.matrix.rows() == system_info.system_data.Sigma.matrix.cols(),
			ExcMessage("Incorrect dimension"));

// we use a temporary matrix to determine whether inflow or outflow
		Sparse_matrix B_temp;

		Assert(b_id == 0 || b_id == 1,ExcMessage("The poisson problem cannot have any other boundary ids"));

		B_temp = system_info.system_data.B.matrix;

		for (unsigned int q = 0 ; q < fe_v.n_quadrature_points ; q++)
		{
			const double jacobian_value = J[q];

			boundary_rhs_value = 0;                 


			system_info.bcrhs_wall.BCrhs(fe_v.quadrature_point(q),fe_v.normal_vector(q),
					boundary_rhs_value,b_id);


			Sparse_matrix Projector = system_info.build_Projector(fe_v.normal_vector(q));



  // Simga(as given in PDF) * B * Projector
  // Sigma in the PDF = Projector.transpose * Sigma(In the code)
			Full_matrix Sigma_B_P =       Projector.transpose()
			* system_info.system_data.Sigma.matrix 
			* B_temp
			* Projector;

			Full_matrix Sigma = Projector.transpose()
			* system_info.system_data.Sigma.matrix ;

			for (unsigned int i = 0 ; i < dofs_per_cell ; i ++)
			{
				const double shape_value_test = fe_v.shape_value(i,q);
				for (unsigned int j = 0 ; j < dofs_per_cell ; j ++)
					cell_matrix(i,j) -= shape_value_test
				* Sigma_B_P(component[i],component[j])
				* fe_v.shape_value(j,q) 
				* jacobian_value;                                    


				for (unsigned int j = 0 ; j < boundary_rhs_value.size() ; j++)
					cell_rhs(i) -= shape_value_test 
				* Sigma(component[i],j) 
				* boundary_rhs_value[j] 
				* jacobian_value;

			}


      }



	}

  template<int dim>
  void 
  Assembly_Manager_Periodic<dim>
  ::integrate_boundary_manuel_kinetic(FullMatrix<double> &cell_matrix,
                                  Vector<double> &cell_rhs,
                                  const FEValuesBase<dim> &fe_v,
                                  std::vector<double> &J,
                                  const unsigned int b_id)
  {
      std::vector<unsigned int> component(dofs_per_cell);

    for (unsigned int i = 0 ; i < dofs_per_cell ; i ++)
      component[i] = this->fe_data_structure.finite_element.system_to_component_index(i).first;

    Vector<double> boundary_rhs_value;
    boundary_rhs_value.reinit(nEqn);
    
    Assert(system_info.system_data.B.matrix.rows() == system_info.system_data.Sigma.matrix.cols(),
      ExcMessage("Incorrect dimension"));

    // we use a temporary matrix to determine whether inflow or outflow
    Sparse_matrix rho_temp;

    Assert(b_id == 0 || b_id == 1,ExcMessage("The poisson problem cannot have any other boundary ids"));

    rho_temp = system_info.system_data.rhoW.matrix;

    for (unsigned int q = 0 ; q < fe_v.n_quadrature_points ; q++)
    {
      const double jacobian_value = J[q];

      boundary_rhs_value = 0;                 
      Tensor<1,dim> outward_normal = fe_v.normal_vector(q);

      system_info.bcrhs_wall_kinetic.BCrhs(fe_v.quadrature_point(q),outward_normal,
                                            boundary_rhs_value,b_id);

      Sparse_matrix Projector = system_info.build_Projector(outward_normal);

      Eigen::MatrixXd Am = system_info.build_Aminus_kinetic(outward_normal);

      Eigen::MatrixXd Am_rho = Am * rho_temp * Projector;

      for (unsigned int i = 0 ; i < dofs_per_cell ; i ++)
      {
        const double shape_value_test = fe_v.shape_value(i,q);
        for (unsigned int j = 0 ; j < dofs_per_cell ; j ++)
          cell_matrix(i,j) += 0.5 * shape_value_test
                              * (Am(component[i],component[j]) - Am_rho(component[i],component[j]))
                              * fe_v.shape_value(j,q) 
                              * jacobian_value;                                    


        for (unsigned int j = 0 ; j < boundary_rhs_value.size() ; j++)
          cell_rhs(i) += 0.5 * shape_value_test 
                        * Am(component[i],j) 
                        * boundary_rhs_value[j] 
                        * jacobian_value;

      }


      }

  }

			// integrate the face term using manuel computation
	template<int dim>
	void 
	Assembly_Manager_Periodic<dim>
	::integrate_face_manuel(FullMatrix<double> &u1_v1,
							FullMatrix<double> &u1_v2,
							FullMatrix<double> &u2_v1,
							FullMatrix<double> &u2_v2,
							const FEValuesBase<dim> &fe_v,
							const FEValuesBase<dim> &fe_v_neighbor,
							std::vector<double> &J)
	{


		// both the cell and its neighbor must have the same finite element object
		const FiniteElement<dim> &fe_in_cell1 = this->fe_data_structure.finite_element;
		const FiniteElement<dim> &fe_in_cell2 = this->fe_data_structure.finite_element;

		const unsigned int dofs_per_cell1 = fe_in_cell1.dofs_per_cell;
		const unsigned int n_components1 = fe_in_cell1.n_components();
		const unsigned int dofs_per_component1 = dofs_per_cell1/n_components1;

		const unsigned int dofs_per_cell2 = fe_in_cell2.dofs_per_cell;
		const unsigned int n_components2 = fe_in_cell2.n_components();
		const unsigned int dofs_per_component2 = dofs_per_cell2/n_components2;


    //CAUTION: ASSUMPTION OF STRAIGHT EDGES IN THE INTERIOR
		Tensor<1,dim> outward_normal = fe_v.normal_vector(0);
		Eigen::MatrixXd Am = system_info.build_Aminus(outward_normal);
		Eigen::MatrixXd Am_neighbor = system_info.build_Aminus(-outward_normal);


		FullMatrix<double> Mass_u1_v1(dofs_per_component1,dofs_per_component1);
		FullMatrix<double> Mass_u2_v1(dofs_per_component1,dofs_per_component2);
		FullMatrix<double> Mass_u2_v2(dofs_per_component2,dofs_per_component2);
		FullMatrix<double> Mass_u1_v2(dofs_per_component2,dofs_per_component1);

		Mass_u1_v1 = this->Compute_Mass_cell_neighbor(fe_v,
			fe_v,dofs_per_component1,
			dofs_per_component1,J);

		Mass_u2_v1 = this->Compute_Mass_cell_neighbor(fe_v,
			fe_v_neighbor,dofs_per_component1,
			dofs_per_component2,J);

		Mass_u2_v2 =  this->Compute_Mass_cell_neighbor(fe_v_neighbor,
			fe_v_neighbor,dofs_per_component2,
			dofs_per_component2,J);

		Mass_u1_v2 =  this->Compute_Mass_cell_neighbor(fe_v_neighbor,
			fe_v,dofs_per_component2,
			dofs_per_component1,J);

		u1_v1 = matrix_opt.multiply_scalar(0.5,matrix_opt.compute_A_outer_B(Am,Mass_u1_v1));
		u2_v1 = matrix_opt.multiply_scalar(-0.5,matrix_opt.compute_A_outer_B(Am,Mass_u2_v1));
		u1_v2 = matrix_opt.multiply_scalar(-0.5,matrix_opt.compute_A_outer_B(Am_neighbor,Mass_u1_v2));
		u2_v2 = matrix_opt.multiply_scalar(0.5,matrix_opt.compute_A_outer_B(Am_neighbor,Mass_u2_v2));

	}
}