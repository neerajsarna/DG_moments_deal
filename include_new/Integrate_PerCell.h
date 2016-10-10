template<int dim>
void 
Base_Solver<dim>
::integrate_cell_manuel(FullMatrix<double> &cell_matrix,
                        Vector<double> &cell_rhs,
                        FEValuesBase<dim> &fe_v,
                        std::vector<double> &J,
                        std::vector<Vector<double>> &source_term_value,
                        const typename DoFHandler<dim>::active_cell_iterator &cell)
{
  Assert(cell_matrix.m() !=0 || cell_matrix.n() !=0 ,ExcNotInitialized());
  Assert(cell_matrix.m() == cell_matrix.n() ,ExcMessage("different dof in same cell"));
  Assert(cell_rhs.size() != 0,ExcNotInitialized());
  Assert(J.size() !=0,ExcNotInitialized());
  Assert(source_term_value.size() != 0,ExcNotInitialized());

  const FiniteElement<dim> &fe_in_cell = cell->get_fe();
  const unsigned int total_ngp = J.size();
  const unsigned int dofs_per_cell = fe_in_cell.dofs_per_cell;

  // the following variable is not 
  const unsigned int components_per_cell = fe_in_cell.n_components();
  const unsigned int indices_per_cell = dofs_per_cell/components_per_cell;

  std::vector<std::vector<double>> component_to_system(components_per_cell,std::vector<double> (indices_per_cell));

  for (unsigned int i = 0 ; i < components_per_cell ; i ++)
    for (unsigned int j = 0 ; j < indices_per_cell ; j ++)
      component_to_system[i][j] = fe_in_cell.component_to_system_index(i,j); 

  // loop over all the gauss points
  for (unsigned int q = 0 ; q < total_ngp ; q++)
  {
    // value of the jacobian at the gauss point
    double jacobian_value = J[q];
    
    // index is the local id of the degree of freedom
    for (unsigned int index_test = 0 ; index_test < indices_per_cell ; index_test++)
    {
      for (unsigned int index_sol = 0 ; index_sol < indices_per_cell ; index_sol++ )
      {
      // instead of looping over all the components for all the indices, we only loop over the 
      // non-zero values of the sparse matrix A
      for (unsigned int space = 0 ; space < dim ; space++)
          for (unsigned int m = 0 ; m < system_info->system_data.A[space].matrix.outerSize(); m++)
        {
          const int dof_test = component_to_system[m][index_test];
          const double shape_value_test = fe_v.shape_value(dof_test,q);

          for (Sparse_matrix::InnerIterator n(system_info->system_data.A[space].matrix,m); n ; ++n)
          {
            int dof_sol = component_to_system[n.col()][index_sol];
            double grad_value_sol = fe_v.shape_grad(dof_sol,q)[space];

            cell_matrix(dof_test,dof_sol) += shape_value_test * n.value()
                                             * grad_value_sol * jacobian_value;  
          }
        }
      for (unsigned int m = 0 ; m < system_info->system_data.P.matrix.outerSize(); m++)
      {
        const int dof_test = component_to_system[m][index_test];
        const double shape_value_test = fe_v.shape_value(dof_test,q);

        for (Sparse_matrix::InnerIterator n(system_info->system_data.P.matrix,m); n ; ++n)
          {
            const int dof_sol = component_to_system[n.col()][index_sol];
            const double shape_value_sol = fe_v.shape_value(dof_sol,q);

            cell_matrix(dof_test,dof_sol) += shape_value_test * n.value()
                                             * shape_value_sol * jacobian_value; 
          }
      }

     }

     // for the right hand side we iterate over all the components of the test function
     for (unsigned int m = 0 ; m < components_per_cell ; m++)
     {
       const int dof_test = component_to_system[m][index_test];
       const double shape_value_test = fe_v.shape_value(dof_test,q);
       cell_rhs(dof_test) += shape_value_test * source_term_value[q][m] * jacobian_value;
     }
   }
    
  }

}


// characteristic implementation of the boundary conditions
template<int dim>
void 
Base_Solver<dim>
::integrate_boundary_manuel_char(FullMatrix<double> &cell_matrix,
                            Vector<double> &cell_rhs,
                            FEValuesBase<dim> &fe_v,
                            std::vector<double> &J,
                            Vector<double> &component,
                            const typename DoFHandler<dim>::active_cell_iterator &cell,
                            const unsigned int b_id)
{
    Assert(cell_matrix.m() !=0 || cell_matrix.n() !=0 ,ExcNotInitialized());
    Assert(cell_matrix.m() == cell_matrix.n() ,ExcMessage("different dof in same cell"));
    Assert(cell_rhs.size() != 0,ExcNotInitialized());
    Assert(J.size() !=0,ExcNotInitialized());
    
    const FiniteElement<dim> &fe_in_cell = cell->get_fe();
   const unsigned int dofs_per_cell = fe_in_cell.dofs_per_cell;

     Vector<double> boundary_rhs_value;

      boundary_rhs_value.reinit(this->constants.nBC);

 for (unsigned int q = 0 ; q < fe_v.n_quadrature_points ; q++)
 {
  const double jacobian_value = J[q];

  boundary_rhs_value = 0;                 

         // the outward normal to the boundary
  Tensor<1,dim> outward_normal = fe_v.normal_vector(q);

  system_info->build_BCrhs(fe_v.quadrature_point(q),outward_normal,
                                    boundary_rhs_value,b_id);

        // build the matrices needed
  Full_matrix Am = system_info->build_Aminus(outward_normal);
  Sparse_matrix Projector = system_info->build_Projector(outward_normal);
  Sparse_matrix Inv_Projector = system_info->build_InvProjector(outward_normal);

  Eigen::MatrixXd Am_invP_B_hat_P = Am * Inv_Projector 
                                    * system_info->B_hat 
                                    * Projector;

  Eigen::MatrixXd Am_invP_X_min_B_tild_inv = Am
                                            * Inv_Projector 
                                            * system_info->X_minus
                                            * system_info->B_tilde_inv;


  for (unsigned int i = 0 ; i < dofs_per_cell ; i ++)
  {
    const double shape_value_test = fe_v.shape_value(i,q);
    for (unsigned int j = 0 ; j < dofs_per_cell ; j ++)
      cell_matrix(i,j) += 0.5 * shape_value_test
                        * Am_invP_B_hat_P(component[i],component[j])
                        * fe_v.shape_value(j,q) 
                        * jacobian_value;                                    


    for (unsigned int j = 0 ; j < Am_invP_X_min_B_tild_inv.cols() ; j++)
     cell_rhs(i) += 0.5 * shape_value_test 
                    * Am_invP_X_min_B_tild_inv(component[i],j)
                     * boundary_rhs_value[j] * jacobian_value;

 }

 }
}


// odd implementation of the boundary conditions

template<int dim>
void 
Base_Solver<dim>
::integrate_boundary_manuel_odd(FullMatrix<double> &cell_matrix,
                            Vector<double> &cell_rhs,
                            FEValuesBase<dim> &fe_v,
                            std::vector<double> &J,
                            Vector<double> &component,
                            const typename DoFHandler<dim>::active_cell_iterator &cell,
                            const unsigned int b_id)
{
    Assert(cell_matrix.m() !=0 || cell_matrix.n() !=0 ,ExcNotInitialized());
    Assert(cell_matrix.m() == cell_matrix.n() ,ExcMessage("different dof in same cell"));
    Assert(cell_rhs.size() != 0,ExcNotInitialized());
    Assert(J.size() !=0,ExcNotInitialized());
    
    const FiniteElement<dim> &fe_in_cell = cell->get_fe();
    const unsigned int dofs_per_cell = fe_in_cell.dofs_per_cell;

    Vector<double> boundary_rhs_value;

    boundary_rhs_value.reinit(this->constants.nEqn);

 for (unsigned int q = 0 ; q < fe_v.n_quadrature_points ; q++)
 {
  const double jacobian_value = J[q];

  boundary_rhs_value = 0;                 
  system_info->build_BCrhs(fe_v.quadrature_point(q),fe_v.normal_vector(q),
                          boundary_rhs_value,b_id);

        // build the matrices needed
  Full_matrix Am = system_info->build_Aminus(fe_v.normal_vector(q));
  Sparse_matrix Projector = system_info->build_Projector(fe_v.normal_vector(q));
  Sparse_matrix Inv_Projector = system_info->build_InvProjector(fe_v.normal_vector(q));

  Full_matrix Am_invP_BC_P = Am * Inv_Projector 
                               * system_info->BC
                               * Projector;

  Full_matrix Am_invP = Am * Inv_Projector;

  for (unsigned int i = 0 ; i < dofs_per_cell ; i ++)
  {
    const double shape_value_test = fe_v.shape_value(i,q);
    for (unsigned int j = 0 ; j < dofs_per_cell ; j ++)
      cell_matrix(i,j) += 0.5 * shape_value_test
                          * (Am(component[i],component[j])-Am_invP_BC_P(component[i],component[j]))
                          * fe_v.shape_value(j,q) * jacobian_value;                                    


    for (unsigned int j = 0 ; j < Am_invP.cols() ; j++)
     cell_rhs(i) += 0.5 * shape_value_test 
                    * Am_invP(component[i],j) * boundary_rhs_value[j] * jacobian_value;

 }

 }
}


template<int dim>
void 
Base_Solver<dim>
::integrate_face_manuel(FullMatrix<double> &u1_v1,
                        FullMatrix<double> &u1_v2,
                        FullMatrix<double> &u2_v1,
                        FullMatrix<double> &u2_v2,
                        FEValuesBase<dim> &fe_v,
                        FEValuesBase<dim> &fe_v_neighbor,
                        std::vector<double> &J,
                        Vector<double> &component,
                        const typename DoFHandler<dim>::active_cell_iterator &cell)
{

    Assert(u1_v1.m() !=0 || u1_v1.n() !=0 ,ExcNotInitialized());
    Assert(u1_v2.m() !=0 || u1_v2.n() !=0 ,ExcNotInitialized());
    Assert(u2_v1.m() !=0 || u2_v1.n() !=0 ,ExcNotInitialized());
    Assert(u2_v2.m() !=0 || u2_v2.n() !=0 ,ExcNotInitialized());
    Assert(J.size() !=0,ExcNotInitialized());
    
    const FiniteElement<dim> &fe_in_cell = cell->get_fe();
  const unsigned int dofs_per_cell = fe_in_cell.dofs_per_cell;

    Assert(component.size() == dofs_per_cell, ExcMessage("value mismatch"));

    for (unsigned int q = 0 ; q < fe_v.n_quadrature_points ; q++)
    {
      Tensor<1,dim> outward_normal = fe_v.normal_vector(q);
      // build the matrices needed
      Eigen::MatrixXd Am = system_info->build_Aminus(outward_normal);
      Eigen::MatrixXd Am_neighbor = system_info->build_Aminus(-outward_normal);
      double jacobian_value = J[q];

      for (unsigned int i = 0 ; i < dofs_per_cell ; i ++)
      {
        const double shape_value_test = fe_v.shape_value(i,q);
        const double shape_value_test_neighbor = fe_v_neighbor.shape_value(i,q);

        for (unsigned int j = 0 ; j < dofs_per_cell ; j ++)
        {
         const double shape_value_sol = fe_v.shape_value(j,q);
         const double shape_value_neighbor = fe_v_neighbor.shape_value(j,q);

         u1_v1(i,j) += 0.5 * shape_value_test * Am(component[i],component[j]) 
                        * shape_value_sol * jacobian_value;

         u2_v1(i,j) -= 0.5 * shape_value_test * Am(component[i],component[j]) 
                      * shape_value_neighbor * jacobian_value;

         u2_v2(i,j) += 0.5 * shape_value_test_neighbor * Am_neighbor(component[i],component[j])
                      * shape_value_neighbor * jacobian_value;

         u1_v2(i,j) -= 0.5 * shape_value_test_neighbor * Am_neighbor(component[i],component[j]) 
                      * shape_value_sol * jacobian_value;

       }
     }
   }
}
