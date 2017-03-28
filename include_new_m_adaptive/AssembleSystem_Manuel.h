
// implementation of the odd boundary conditions
template<int dim> 
void 
Base_Solver<dim>::assemble_system_odd()
  {
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
      typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      typename hp::DoFHandler<dim>::cell_iterator neighbor;

      // std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      // std::vector<types::global_dof_index> local_dof_indices_neighbor(dofs_per_cell);

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
      

  for (; cell != endc ; cell++) 
      {
        const unsigned int fe_index = cell->active_fe_index();
        const unsigned int dofs_this_cell = dofs_per_cell[fe_index];

        Assert(fe_index <= max_fe_index ,ExcNotImplemented());
        cell_rhs[fe_index] = 0;


        // the mapping from the dofs of the present cell to the global dofs
        std::vector<types::global_dof_index> local_dof_indices(dofs_this_cell);
        cell->get_dof_indices(local_dof_indices);

        
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

                      
                // assemble the terms on the boundary
                for (unsigned int i = 0 ; i < dofs_per_cell[fe_index] ; i++)
                  for (unsigned int j = 0 ; j < dofs_per_cell[fe_index] ; j++)
                   global_matrix.add(local_dof_indices[i],local_dof_indices[j],boundary_matrix[fe_index](i,j));
                
              }
             
              
                // if the face is not at the boundary then we need to integrate over the face
                else
               {            
                
                 neighbor = cell->neighbor(face);


                 const unsigned int fe_index_neighbor = neighbor->active_fe_index();
                 const unsigned int dofs_neighbor = dofs_per_cell[fe_index_neighbor];
                 std::vector<types::global_dof_index> local_dof_indices_neighbor(dofs_neighbor);

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

                   // assembly for u1_v1
                   for (unsigned int i = 0 ; i < dofs_this_cell ; i++)
                    for (unsigned int j = 0 ; j < dofs_this_cell ; j++)
                      global_matrix.add(local_dof_indices[i],local_dof_indices[j],u1_v1(i,j));

                  // assembly for u2_v2
                   for (unsigned int i = 0 ; i < dofs_neighbor ; i++)
                    for (unsigned int j = 0 ; j < dofs_neighbor ; j++)
                      global_matrix.add(local_dof_indices_neighbor[i],local_dof_indices_neighbor[j],u2_v2(i,j));

                  // assembly for u2_v1
                   for (unsigned int i = 0 ; i < dofs_this_cell ; i++)
                    for (unsigned int j = 0 ; j < dofs_neighbor ; j++)
                      global_matrix.add(local_dof_indices[i],local_dof_indices_neighbor[j],u2_v1(i,j));

                  // assembly for u1_v2
                  for (unsigned int i = 0 ; i < dofs_neighbor ; i++)
                    for (unsigned int j = 0 ; j < dofs_this_cell ; j++)
                      global_matrix.add(local_dof_indices_neighbor[i],local_dof_indices[j],u1_v2(i,j)) ;
                   
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


                   // assembly for u1_v1
                   for (unsigned int i = 0 ; i < dofs_this_cell ; i++)
                    for (unsigned int j = 0 ; j < dofs_this_cell ; j++)
                      global_matrix.add(local_dof_indices[i],local_dof_indices[j],u1_v1(i,j));

                  // assembly for u2_v2
                   for (unsigned int i = 0 ; i < dofs_neighbor ; i++)
                    for (unsigned int j = 0 ; j < dofs_neighbor ; j++)
                      global_matrix.add(local_dof_indices_neighbor[i],local_dof_indices_neighbor[j],u2_v2(i,j));

                  // assembly for u2_v1
                   for (unsigned int i = 0 ; i < dofs_this_cell ; i++)
                    for (unsigned int j = 0 ; j < dofs_neighbor ; j++)
                      global_matrix.add(local_dof_indices[i],local_dof_indices_neighbor[j],u2_v1(i,j));

                  // assembly for u1_v2
                  for (unsigned int i = 0 ; i < dofs_neighbor ; i++)
                    for (unsigned int j = 0 ; j < dofs_this_cell ; j++)
                      global_matrix.add(local_dof_indices_neighbor[i],local_dof_indices[j],u1_v2(i,j)) ;                
  
                 }


               }
            }
              

            
              for (unsigned int i = 0 ; i < dofs_per_cell[fe_index] ; i++)
                  for (unsigned int j = 0 ; j < dofs_per_cell[fe_index] ; j++)          
                    global_matrix.add(local_dof_indices[i],local_dof_indices[j],cell_matrix[fe_index](i,j));

              
             for(unsigned int i = 0 ; i < dofs_per_cell[fe_index] ; i++)
                  system_rhs(local_dof_indices[i]) += cell_rhs[fe_index](i);

  
      }
}



  
  
