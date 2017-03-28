// the following routine will remain the same untill and unless we choose to refine the grid. 
// If the grid refinement is left as it is then the the following routine will remain the same. 
template<int dim>
void 
Base_Solver<dim>::distribute_dof_residual_computation()
{


    // the maximum fe index should not be greater than the total number of systems given
    AssertIndexRange(max_fe_index_residual,nEqn.size());

	  typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler_max_moments.begin_active(),
	  												     endc = dof_handler_max_moments.end();

	  // fe index corresponding to the moment method which we wish to approximate
	  for (; cell != endc ; cell++)
	  	cell->set_active_fe_index(max_fe_index_residual);

	  dof_handler_max_moments.distribute_dofs(finite_element);
	  VelocitySpace_residual.reinit(dof_handler_max_moments.n_dofs());	
	  solution_max_moments.reinit(dof_handler_max_moments.n_dofs());
    VelocitySpace_error_per_cell.reinit(this->triangulation.n_active_cells());

            // the vector which stores the residual 
        DynamicSparsityPattern dsp(dof_handler_max_moments.n_dofs(),dof_handler_max_moments.n_dofs());


  //std::cout << "making flux sparsity pattern " << std::endl;
        DoFTools::make_flux_sparsity_pattern (dof_handler_max_moments, dsp);
        system_rhs.reinit(dof_handler_max_moments.n_dofs());
  
        global_matrix.reinit(dsp);  

}

template<int dim>
void 
Base_Solver<dim>::develop_solution_max_moments()
{
	
	// now we interpolate between the two different dof objects. The solution object from a mixed moment method
	// is now projected onto a space which corresponds to the maximum number of moments
	FETools::interpolate(dof_handler,
							 solution,
							 dof_handler_max_moments,
							 solution_max_moments);


	 //The initialization of the residual is necessary because we add up the values from different cells.
	 VelocitySpace_residual = 0;


}


template<int dim>
void 
Base_Solver<dim>::compute_residual()
{
	  // we distribute the dofs for residual computation
	  distribute_dof_residual_computation();
	  develop_solution_max_moments();



     // we check whether the moment systems are using the same polynomial degree or not
    for (unsigned long int i = 0 ; i < nEqn.size()-1; i++)
      AssertDimension(finite_element[i].dofs_per_cell/nEqn[i],finite_element[i+1].dofs_per_cell/nEqn[i+1]);

      const QGauss<dim> quadrature_basic(ngp);
      const QGauss<dim-1> face_quadrature_basic(ngp_face);

      // integration on the volume
      hp::QCollection<dim> quadrature;

      // integration on the surface
      hp::QCollection<dim-1> face_quadrature;

      for (unsigned long int i = 0 ; i < nEqn.size() ; i++)
      {
        quadrature.push_back(quadrature_basic);
        face_quadrature.push_back(face_quadrature_basic);
      }


      const UpdateFlags update_flags               =  update_gradients
                                                     | update_q_points
                                                     | update_JxW_values
                                                     | update_values,

      face_update_flags          = update_values
      | update_q_points
      | update_JxW_values
      | update_normal_vectors,
      neighbor_face_update_flags = update_values;

      hp::FEValues<dim>  hp_fe_v(mapping,finite_element,quadrature, update_flags);
      hp::FEFaceValues<dim> hp_fe_v_face(mapping,finite_element, face_quadrature, face_update_flags);
      hp::FEFaceValues<dim> hp_fe_v_face_neighbor(mapping,finite_element, face_quadrature, neighbor_face_update_flags);
      hp::FESubfaceValues<dim> hp_fe_v_subface_neighbor(mapping,finite_element, face_quadrature, neighbor_face_update_flags);  


      // the total number of quadrature points are the same for both the quadrature routines
      const unsigned int total_ngp = quadrature_basic.size();
      const unsigned int total_ngp_face = face_quadrature_basic.size();

      // iterator over the cells 
      typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler_max_moments.begin_active(),
      												     endc = dof_handler_max_moments.end();
      typename hp::DoFHandler<dim>::cell_iterator neighbor;

      // it is much more efficient to declare the memory and then use it again and again
      std::vector<FullMatrix<double>> cell_matrix;                  // contribution from the interior
      std::vector<FullMatrix<double>> boundary_matrix;
      std::vector<Vector<double>>   cell_rhs;                                   // rhs from the current cell


      for (unsigned long int i = 0 ; i < nEqn.size() ; i++)
      {
        cell_matrix.push_back(FullMatrix<double>(dofs_per_cell[i],dofs_per_cell[i]));
        boundary_matrix.push_back(FullMatrix<double>(dofs_per_cell[i],dofs_per_cell[i]));
        cell_rhs.push_back(Vector<double>(dofs_per_cell[i]));
      }

      // std::vectors to make computations faster
      std::vector<double> Jacobians_interior(total_ngp);
      std::vector<double> Jacobian_face(total_ngp_face);


      //component to system for all the finite element objects being considered
      std::vector<std::vector<Vector<double>>> component_to_system(nEqn.size());

      // loop over all the different systems available in the system
      for (unsigned long int i = 0 ; i < nEqn.size(); i++)
        // loop over all the equations of the particular system
        for (int j = 0 ; j < nEqn[i]; j++)
          component_to_system[i].push_back(Vector<double>(dofs_per_component));

      // now we allocate the values for component_to_system
      for (unsigned long int i = 0 ; i < nEqn.size(); i++)
          for (int k = 0 ; k < nEqn[i] ;k ++)
              for (unsigned int j = 0 ; j < dofs_per_component ; j++)
                component_to_system[i][k](j) = finite_element[i].component_to_system_index(k,j);
           

      // source term
      std::vector<std::vector<Vector<double>>> source_term_value(nEqn.size());

      for (unsigned long int i = 0 ; i < nEqn.size() ; i ++)
        for (unsigned int j = 0 ; j < total_ngp ; j++)
          source_term_value[i].push_back(Vector<double>(nEqn[i]));
      
     
      // need to initialize fe_v so as to save the values of the shape functions
      hp_fe_v.reinit(cell);
      this->Compute_Shape_Value(hp_fe_v,dofs_per_component);

      // value of the dofs for this particular cell
      Vector<double> dofs_values_this_cell;

      // value of the dofs for the neighbor
      Vector<double> dofs_values_neighbor_cell;

      // local residual of the cell
      Vector<double> local_residual_cell;

      // local residual of the neighbor
      Vector<double> local_residual_neighbor;

      unsigned int counter = 0;



      
for (; cell != endc ; cell++) 
      {
        const unsigned int fe_index = cell->active_fe_index();

        // total number of degrees of freedom in the present cell
        const unsigned int dofs_this_cell = dofs_per_cell[fe_index];

       

        Assert(fe_index <=4 ,ExcNotImplemented());
        cell_rhs[fe_index] = 0;


        // the mapping from the dofs of the present cell to the global dofs
        std::vector<types::global_dof_index> local_dof_indices(dofs_this_cell);
        cell->get_dof_indices(local_dof_indices);

		dofs_values_this_cell.reinit(dofs_this_cell);
        for (unsigned int i = 0 ; i < dofs_this_cell ; i++)
        	dofs_values_this_cell(i) = solution_max_moments(local_dof_indices[i]);

        local_residual_cell.reinit(dofs_this_cell);
        local_residual_cell = 0;

        
        hp_fe_v.reinit(cell);
        const FEValues<dim> &fe_v = hp_fe_v.get_present_fe_values();        

        Jacobians_interior = fe_v.get_JxW_values();
        system_info[fe_index].source_term(fe_v.get_quadrature_points(),source_term_value[fe_index]);

        
        integrate_cell_manuel(cell_matrix[fe_index],
                              cell_rhs[fe_index],
                              fe_v,
                              Jacobians_interior,
                              source_term_value[fe_index],
                              component_to_system[fe_index],
                              fe_index);
     


            for(unsigned int face  = 0; face< GeometryInfo<dim>::faces_per_cell; face++ )
            {
            
              hp_fe_v_face.reinit(cell,face);
              const FEFaceValues<dim> &fe_v_face = hp_fe_v_face.get_present_fe_values();
              const typename hp::DoFHandler<dim>::face_iterator face_itr = cell->face(face);
              

              Jacobian_face = fe_v_face.get_JxW_values();
              
              boundary_matrix[fe_index] = 0;

           
              if (face_itr->at_boundary())
              { 
             

                // id of the boundary
                const unsigned int b_id = face_itr->boundary_id();

                      
                integrate_boundary_manuel_odd(boundary_matrix[fe_index], 
                                              cell_rhs[fe_index],
                                              fe_v_face, 
                                              Jacobian_face,
                                              component_to_system[fe_index], 
                                              b_id,
                                              fe_index); 


                // the local residual from the cell 
        		local_residual_cell = matrix_opt.Sparse_matrix_dot_Vector(boundary_matrix[fe_index],dofs_values_this_cell);

        	   //The residual in the present cell = cell_matrix * dofs in this cell - cell_rhs
        		// the right hand side will be dealth with later
        		for (unsigned int i = 0 ; i < dofs_this_cell ; i++)
        				VelocitySpace_residual(local_dof_indices[i]) += local_residual_cell(i);
                     
                
              }
             
              
                // if the face is not at the boundary then we need to integrate over the face
                else
               {            
                
                 neighbor = cell->neighbor(face);


                 const unsigned int fe_index_neighbor = neighbor->active_fe_index();
                 const unsigned int dofs_neighbor = dofs_per_cell[fe_index_neighbor];
                 std::vector<types::global_dof_index> local_dof_indices_neighbor(dofs_neighbor);


                 dofs_values_neighbor_cell.reinit(dofs_neighbor);
        		 local_residual_neighbor.reinit(dofs_neighbor);
        		 local_residual_neighbor = 0;


                 // we now allocate the memory for the matrices which will store the integrals

                 FullMatrix<double> u1_v1(dofs_this_cell,dofs_this_cell);
                 FullMatrix<double> u1_v2(dofs_neighbor,dofs_this_cell);
                 FullMatrix<double> u2_v1(dofs_this_cell,dofs_neighbor);
                 FullMatrix<double> u2_v2(dofs_neighbor,dofs_neighbor);

                 if (cell->neighbor_is_coarser(face))
                 {

                   Assert(!cell->has_children(), ExcInternalError());
                   Assert(!neighbor->has_children(), ExcInternalError());

                   neighbor->get_dof_indices(local_dof_indices_neighbor);

                   
        		   for (unsigned int i = 0 ; i < dofs_neighbor ; i++)
        					dofs_values_neighbor_cell(i) = solution_max_moments(local_dof_indices_neighbor[i]);

                   const std::pair<unsigned int, unsigned int> neighbor_face_no
                   											  = cell->neighbor_of_coarser_neighbor(face);


                   hp_fe_v_subface_neighbor.reinit(neighbor,neighbor_face_no.first,neighbor_face_no.second);

                   const FESubfaceValues<dim> &fe_v_subface_neighbor = hp_fe_v_subface_neighbor.get_present_fe_values();
                   
                   integrate_face_manuel(u1_v1,
                                        u1_v2,
                                        u2_v1,
                                        u2_v2,
                                        fe_v_face,
                                        fe_v_subface_neighbor,
                                        Jacobian_face,
                                        fe_index,
                                        fe_index_neighbor);

                   // now we compute the residual 
                   // contribution from this cell into the residual of this cell
                   local_residual_cell = matrix_opt.Sparse_matrix_dot_Vector(u1_v1,dofs_values_this_cell);

                   // contribution from the neighbouring cell into the residual of this cell
                   local_residual_cell += matrix_opt.Sparse_matrix_dot_Vector(u2_v1,dofs_values_neighbor_cell);

                   // contribution from this cell into the residual of the neighboring cell
                   local_residual_neighbor = matrix_opt.Sparse_matrix_dot_Vector(u1_v2,dofs_values_this_cell);

                   // contribution from the neighbor into the residual of the neighboring cell
                   local_residual_neighbor += matrix_opt.Sparse_matrix_dot_Vector(u2_v2,dofs_values_neighbor_cell);

                  
                  // assemble the residual for this cell
                   for (unsigned int i = 0 ; i < dofs_this_cell ; i++)
                   		VelocitySpace_residual(local_dof_indices[i]) += local_residual_cell(i);

                  // assemble the residual for the neighbor
                   	for (unsigned int i = 0 ; i < dofs_neighbor ; i++)
                   		VelocitySpace_residual(local_dof_indices_neighbor[i]) += local_residual_neighbor(i);

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


        		   for (unsigned int i = 0 ; i < dofs_neighbor ; i++)
        					dofs_values_neighbor_cell(i) = solution_max_moments(local_dof_indices_neighbor[i]);


                    const unsigned int neighbor_face_no = cell->neighbor_of_neighbor(face);

                    hp_fe_v_face_neighbor.reinit(neighbor,neighbor_face_no);
                    const FEFaceValues<dim> &fe_v_face_neighbor = hp_fe_v_face_neighbor.get_present_fe_values();


                   integrate_face_manuel(u1_v1,
                                        u1_v2,
                                        u2_v1,
                                        u2_v2,
                                        fe_v_face,
                                        fe_v_face_neighbor,
                                        Jacobian_face,
                                        fe_index,
                                        fe_index_neighbor);


                   // now we compute the residual 
                   // contribution from this cell into the residual of this cell
                   local_residual_cell = matrix_opt.Sparse_matrix_dot_Vector(u1_v1,dofs_values_this_cell);

                   // contribution from the neighbouring cell into the residual of this cell
                   local_residual_cell += matrix_opt.Sparse_matrix_dot_Vector(u2_v1,dofs_values_neighbor_cell);

                   // contribution from this cell into the residual of the neighboring cell
                   local_residual_neighbor = matrix_opt.Sparse_matrix_dot_Vector(u1_v2,dofs_values_this_cell);

                   // contribution from the neighbor into the residual of the neighboring cell
                   local_residual_neighbor += matrix_opt.Sparse_matrix_dot_Vector(u2_v2,dofs_values_neighbor_cell);

                  
                  // assemble the residual for this cell
                   for (unsigned int i = 0 ; i < dofs_this_cell ; i++)
                   		VelocitySpace_residual(local_dof_indices[i]) += local_residual_cell(i);

                  // assemble the residual for the neighbor
                   	for (unsigned int i = 0 ; i < dofs_neighbor ; i++)
                   		VelocitySpace_residual(local_dof_indices_neighbor[i]) += local_residual_neighbor(i);
  
                 }


               }
             }
              

            
            // the local residual from the cell interior
        	local_residual_cell = matrix_opt.Sparse_matrix_dot_Vector(cell_matrix[fe_index],dofs_values_this_cell);

        	//The residual in the present cell = cell_matrix * dofs in this cell - cell_rhs
        	// the right hand side will be dealth with later
        		for (unsigned int i = 0 ; i < dofs_this_cell ; i++)
        				VelocitySpace_residual(local_dof_indices[i]) += local_residual_cell(i)-cell_rhs[fe_index](i);

  
      }

      cell = dof_handler_max_moments.begin_active();

        AssertIndexRange(max_fe_index_residual,nEqn.size());
        const int max_equations = nEqn[max_fe_index_residual];


      for (; cell != endc ; cell++)
      {
        const unsigned int fe_index = cell->active_fe_index();

        // total number of degrees of freedom in the present cell
        const unsigned int dofs_this_cell = dofs_per_cell[fe_index];

       

        Assert(fe_index <=4 ,ExcNotImplemented());
        cell_rhs[fe_index] = 0;


        // the mapping from the dofs of the present cell to the global dofs
        std::vector<types::global_dof_index> local_dof_indices(dofs_this_cell);
        cell->get_dof_indices(local_dof_indices);

        for (int eq = 0 ; eq < max_equations ; eq++)
          if (eq == 2)
            for (unsigned int i = 0 ; i < dofs_per_component ; i++)
              VelocitySpace_error_per_cell(counter) += fabs(VelocitySpace_residual(local_dof_indices[component_to_system[fe_index][eq](i)]));
        
        counter++;

      }

      // In the above computation we have computed the dofs of the residual. Now we will compute the 
      // residual per cell through the following computation
      // error per cell of the domain
        // AssertIndexRange(max_fe_index_residual,nEqn.size());
        // const int max_equations = nEqn[max_fe_index_residual];

        // ComponentSelectFunction<dim> weight(2,max_equations);          // used to compute only the error in theta


        // // compute the l1 norm of the residual
        // VectorTools::integrate_difference (mapping,
        //                                    dof_handler_max_moments,
        //                                    VelocitySpace_residual,
        //                                    ZeroFunction<dim>(max_equations),
        //                                    VelocitySpace_error_per_cell,
        //                                    quadrature,
        //                                    VectorTools::L1_norm,
        //                                    &weight);  


}
