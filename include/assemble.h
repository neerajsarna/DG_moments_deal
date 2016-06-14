#include "integrate_per_cell_manuel.h"
template<int num_flux,int dim> 
void 
Solver_DG<num_flux,dim>::assemble_system()
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

      FEValues<dim>  fe_v(mapping,finite_element,quadrature, update_flags);
      FEFaceValues<dim> fe_v_face(mapping,finite_element, face_quadrature, face_update_flags);
      FEFaceValues<dim> fe_v_face_neighbor(mapping,finite_element, face_quadrature, neighbor_face_update_flags);
      FESubfaceValues<dim> fe_v_subface_neighbor(mapping,finite_element, face_quadrature, neighbor_face_update_flags);  

      const unsigned int dofs_per_cell = finite_element.dofs_per_cell;
      const unsigned int total_ngp = quadrature.size();
      const unsigned int total_ngp_face = face_quadrature.size();

      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      typename DoFHandler<dim>::active_cell_iterator endc_trial = cell++;
      typename DoFHandler<dim>::cell_iterator neighbor;

      vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      vector<types::global_dof_index> local_dof_indices_neighbor(dofs_per_cell);

      FullMatrix<double> cell_matrix(dofs_per_cell,dofs_per_cell);                  // contribution from the interior
      Vector<double>     cell_rhs(dofs_per_cell);                                   // rhs from the current cell

      // matrices for the boundary terms
      FullMatrix<double> u1_v1(dofs_per_cell,dofs_per_cell);      // relating current to current
      FullMatrix<double> u1_v2(dofs_per_cell,dofs_per_cell);  // relating current to neighbor
      FullMatrix<double> u2_v1(dofs_per_cell,dofs_per_cell);  // relating neighbor to current
      FullMatrix<double> u2_v2(dofs_per_cell,dofs_per_cell);  // relating neighbor to current

      // vectors to make computations faster
      vector<double> Jacobians_interior(total_ngp);
      vector<double> Jacobian_face(total_ngp_face);

      // component corresponding to every dof
      Vector<double> component(dofs_per_cell);
      vector<Vector<double>> source_term_value(total_ngp,Vector<double>(this->nEqn));

      for (unsigned int i = 0 ; i < dofs_per_cell ; i++)
        component(i) = finite_element.system_to_component_index(i).first;
      
      // end of vector to make computation faster

      for (; cell != endc_trial ; cell++) 
      {
        cell_matrix = 0;
        cell_rhs = 0;
        cell->get_dof_indices(local_dof_indices);

        fe_v.reinit(cell);
        Jacobians_interior = fe_v.get_JxW_values();
        equation_system_data->source_term(fe_v.get_quadrature_points(),source_term_value,solve_system);

        integrate_cell_manuel(cell_matrix,cell_rhs,
                              fe_v,Jacobians_interior,
                              source_term_value,cell);
	
        

            for(unsigned int face  = 0; face< GeometryInfo<dim>::faces_per_cell; face++ )
            {
            
              fe_v_face.reinit(cell,face);
              Jacobian_face = fe_v_face.get_JxW_values();
              
              u1_v1 = 0;
              u1_v2 = 0;
              u2_v1 = 0;
              u2_v2 = 0;


              if (cell->face(face)->at_boundary())
              { 
                integrate_boundary_manuel(cell_matrix, cell_rhs,
                                          fe_v_face, Jacobian_face,
                                          component, cell); 
                
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

  
                 }
                }

               for (unsigned int i = 0 ; i < dofs_per_cell ; i++)
              for (unsigned int j = 0 ; j < dofs_per_cell ; j++)
              {
               
                global_matrix.add(local_dof_indices[i],local_dof_indices[j],u1_v1(i,j));
                global_matrix.add(local_dof_indices_neighbor[i],local_dof_indices_neighbor[j],u2_v2(i,j));

                global_matrix.add(local_dof_indices[i],local_dof_indices_neighbor[j],u2_v1(i,j));
                global_matrix.add(local_dof_indices_neighbor[i],local_dof_indices[j],u1_v2(i,j)) ;
              }
            }
              

             for (unsigned int i = 0 ; i < dofs_per_cell ; i++)
              for (unsigned int j = 0 ; j < dofs_per_cell ; j++)
                global_matrix.add(local_dof_indices[i],local_dof_indices[j],cell_matrix(i,j));
              
             for(unsigned int i = 0 ; i < dofs_per_cell ; i++)
              system_rhs(local_dof_indices[i]) += cell_rhs(i);
          }
  }
