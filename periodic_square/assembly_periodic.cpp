template<int dim> 
void 
Base_Solver<dim>::assemble_system_periodic_odd()
  {
      Assert(constants.problem_type == periodic,ExcMessage("Present routine only for periodic squares"));

      // counters 
      unsigned int boundary_wall = 0;
      unsigned int boundary_periodic = 0;
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

      FEValues<dim>  fe_v(mapping,finite_element,quadrature, update_flags);
      FEFaceValues<dim> fe_v_face(mapping,finite_element, face_quadrature, face_update_flags);
      FEFaceValues<dim> fe_v_face_neighbor(mapping,finite_element, face_quadrature, neighbor_face_update_flags);
      FESubfaceValues<dim> fe_v_subface_neighbor(mapping,finite_element, face_quadrature, neighbor_face_update_flags);  

      const unsigned int dofs_per_cell = finite_element.dofs_per_cell;
      const unsigned int total_ngp = quadrature.size();
      const unsigned int total_ngp_face = face_quadrature.size();

      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      typename DoFHandler<dim>::cell_iterator neighbor;

      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      std::vector<types::global_dof_index> local_dof_indices_neighbor(dofs_per_cell);

      Sparse_matrix cell_matrix(dofs_per_cell,dofs_per_cell);                  // contribution from the interior
      FullMatrix<double> boundary_matrix(dofs_per_cell,dofs_per_cell);
      Vector<double>     cell_rhs(dofs_per_cell);                                   // rhs from the current cell

      // matrices for the boundary terms
      Full_matrix u1_v1(dofs_per_cell,dofs_per_cell);      // relating current to current
      Full_matrix u1_v2(dofs_per_cell,dofs_per_cell);  // relating current to neighbor
      Full_matrix u2_v1(dofs_per_cell,dofs_per_cell);  // relating neighbor to current
      Full_matrix u2_v2(dofs_per_cell,dofs_per_cell);  // relating neighbor to current

      // std::vectors to make computations faster
      std::vector<double> Jacobians_interior(total_ngp);
      std::vector<double> Jacobian_face(total_ngp_face);

      // component corresponding to every dof
      Vector<double> component(dofs_per_cell);
      std::vector<Vector<double>> source_term_value(total_ngp,Vector<double>(constants.nEqn));

      for (unsigned int i = 0 ; i < dofs_per_cell ; i++)
        component(i) = finite_element.system_to_component_index(i).first;
      
      this->Compute_Shape_Value(mapping,constants.p,ngp,cell);

      for (; cell != endc ; cell++) 
      {
        cell_rhs = 0;
        cell->get_dof_indices(local_dof_indices);

        fe_v.reinit(cell);
        Jacobians_interior = fe_v.get_JxW_values();
        system_info->source_term(fe_v.get_quadrature_points(),source_term_value);

        integrate_cell_manuel(cell_matrix,cell_rhs,
                              fe_v,Jacobians_interior,
                              source_term_value,cell);
      

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

                    // check whether the cell is at the periodic boundary or at the wall

                    if ( fabs(xcord_face_center - constants.xl) <= 1e-10
                        || fabs(xcord_face_center - constants.xr) <= 1e-10 ) 
                    {
                      neighbor = this->get_periodic_neighbor(xcord_face_center,ycord_cell_center);
                      neighbor->get_dof_indices(local_dof_indices_neighbor);

                      Assert(!neighbor->has_children(),
                             ExcMessage("periodic boundary only to used with zero level meshes"));
                      Assert(fabs(neighbor->center()(1) - ycord_cell_center) < 1e-8,
                            ExcCellCenter(ycord_cell_center,neighbor->center()(1)));

                      const unsigned int neighbor_face = this->get_periodic_neighbor_face(xcord_face_center
                                                                                    ,ycord_cell_center);

                      fe_v_face_neighbor.reinit(neighbor,neighbor_face);


                      // integrate the face between the cell and the periodic neighbor
                      integrate_face_manuel(u1_v1,u1_v2,
                                        u2_v1,u2_v2,
                                        fe_v_face,fe_v_face_neighbor,
                                        Jacobian_face,component,  
                                        cell);  
                
                for (unsigned int i = 0 ; i < dofs_per_cell ; i++)
                  for (unsigned int j = 0 ; j < dofs_per_cell ; j++)
                  {
               
                   Assert(sparsity_pattern.exists(local_dof_indices[i],local_dof_indices[j]),ExcNoElementInSparsity(local_dof_indices[i],local_dof_indices[j]));
                   Assert(sparsity_pattern.exists(local_dof_indices[i],local_dof_indices_neighbor[j]),ExcNoElementInSparsity(local_dof_indices[i],local_dof_indices_neighbor[j]));
                   Assert(sparsity_pattern.exists(local_dof_indices_neighbor[i],local_dof_indices[j]),ExcNoElementInSparsity(local_dof_indices_neighbor[i],local_dof_indices[j]));
                   Assert(sparsity_pattern.exists(local_dof_indices_neighbor[i],local_dof_indices_neighbor[j]),ExcNoElementInSparsity(local_dof_indices_neighbor[i],local_dof_indices_neighbor[j]));

                   global_matrix.add(local_dof_indices[i],local_dof_indices[j],u1_v1(i,j));
                   global_matrix.add(local_dof_indices_neighbor[i],local_dof_indices_neighbor[j],u2_v2(i,j));

                   global_matrix.add(local_dof_indices[i],local_dof_indices_neighbor[j],u2_v1(i,j));
                   global_matrix.add(local_dof_indices_neighbor[i],local_dof_indices[j],u1_v2(i,j)) ;
                  }                  

                     boundary_periodic++;

                    }

                    // if the face is not on the peridic boundary 
                    else 
                    {
                      integrate_boundary_manuel_odd(boundary_matrix, cell_rhs,
                                                      fe_v_face, Jacobian_face,
                                                      component, cell,b_id); 

                      boundary_wall ++;
                                       // assemble the terms on the boundary
                      for (unsigned int i = 0 ; i < dofs_per_cell ; i++)
                        for (unsigned int j = 0 ; j < dofs_per_cell ; j++)
                        {

                         Assert(sparsity_pattern.exists(local_dof_indices[i],local_dof_indices[j]),ExcNoElementInSparsity(local_dof_indices[i],local_dof_indices[j]));
                         global_matrix.add(local_dof_indices[i],local_dof_indices[j],boundary_matrix(i,j));
                        }

                     }  
                
                
                //print_dealii_matrix(boundary_matrix,"boundary matrix");
              }
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
                                         fe_v_face,fe_v_subface_neighbor,
                                         Jacobian_face,component,  
                                         cell);

                   for (unsigned int i = 0 ; i < dofs_per_cell ; i++)
                    for (unsigned int j = 0 ; j < dofs_per_cell ; j++)
                    {

                     Assert(sparsity_pattern.exists(local_dof_indices[i],local_dof_indices[j]),ExcNoElementInSparsity(local_dof_indices[i],local_dof_indices[j]));
                     Assert(sparsity_pattern.exists(local_dof_indices[i],local_dof_indices_neighbor[j]),ExcNoElementInSparsity(local_dof_indices[i],local_dof_indices_neighbor[j]));
                     Assert(sparsity_pattern.exists(local_dof_indices_neighbor[i],local_dof_indices[j]),ExcNoElementInSparsity(local_dof_indices_neighbor[i],local_dof_indices[j]));
                     Assert(sparsity_pattern.exists(local_dof_indices_neighbor[i],local_dof_indices_neighbor[j]),ExcNoElementInSparsity(local_dof_indices_neighbor[i],local_dof_indices_neighbor[j]));

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
                                        Jacobian_face,component,  
                                        cell);  


                  for (unsigned int i = 0 ; i < dofs_per_cell ; i++)
                    for (unsigned int j = 0 ; j < dofs_per_cell ; j++)
                    {
                     
                     Assert(sparsity_pattern.exists(local_dof_indices[i],local_dof_indices[j]),ExcNoElementInSparsity(local_dof_indices[i],local_dof_indices[j]));
                     Assert(sparsity_pattern.exists(local_dof_indices[i],local_dof_indices_neighbor[j]),ExcNoElementInSparsity(local_dof_indices[i],local_dof_indices_neighbor[j]));
                     Assert(sparsity_pattern.exists(local_dof_indices_neighbor[i],local_dof_indices[j]),ExcNoElementInSparsity(local_dof_indices_neighbor[i],local_dof_indices[j]));
                     Assert(sparsity_pattern.exists(local_dof_indices_neighbor[i],local_dof_indices_neighbor[j]),ExcNoElementInSparsity(local_dof_indices_neighbor[i],local_dof_indices_neighbor[j]));

                     global_matrix.add(local_dof_indices[i],local_dof_indices[j],u1_v1(i,j));
                     global_matrix.add(local_dof_indices_neighbor[i],local_dof_indices_neighbor[j],u2_v2(i,j));

                     global_matrix.add(local_dof_indices[i],local_dof_indices_neighbor[j],u2_v1(i,j));
                     global_matrix.add(local_dof_indices_neighbor[i],local_dof_indices[j],u1_v2(i,j)) ;
                   }                   
  
                 }
                }
  
               
            }
              

            for (unsigned int m = 0 ; m < cell_matrix.outerSize(); m++)
                for (Sparse_matrix::InnerIterator n(cell_matrix,m); n ; ++n)
                  global_matrix.add(local_dof_indices[n.row()],local_dof_indices[n.col()],n.value());

             
              
             for(unsigned int i = 0 ; i < dofs_per_cell ; i++)
              system_rhs(local_dof_indices[i]) += cell_rhs(i);

     
    }

    std::cout << "integrated wall " << boundary_wall << " times " << std::endl;
    std::cout << "integraed periodic wall " << boundary_periodic << " times " << std::endl;

  }

