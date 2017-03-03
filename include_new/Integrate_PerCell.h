template<int dim>
void 
Base_Solver<dim>
::integrate_cell_manuel(Sparse_matrix &cell_matrix,
                        Vector<double> &cell_rhs,
                        FEValuesBase<dim> &fe_v,
                        std::vector<double> &J,
                        std::vector<Vector<double>> &source_term_value,
                        const typename DoFHandler<dim>::active_cell_iterator &cell)
{
  Assert(cell_matrix.rows() !=0 || cell_matrix.cols() !=0 ,ExcNotInitialized());
  Assert(cell_matrix.rows() == cell_matrix.cols() ,ExcMessage("different dof in same cell"));
  Assert(cell_rhs.size() != 0,ExcNotInitialized());
  Assert(J.size() !=0,ExcNotInitialized());
  Assert(source_term_value.size() != 0,ExcNotInitialized());

  const FiniteElement<dim> &fe_in_cell = cell->get_fe();
  const unsigned int total_ngp = J.size();
  const unsigned int dofs_per_cell = fe_in_cell.dofs_per_cell;

  // the following variable is not 
  const unsigned int components_per_cell = fe_in_cell.n_components();
  const unsigned int indices_per_cell = dofs_per_cell/components_per_cell;

  AssertDimension(indices_per_cell,this->shape_values.rows());
  AssertDimension(total_ngp,this->shape_values.cols());

  Full_matrix Mass_x(indices_per_cell,indices_per_cell);
  Full_matrix Mass_y(indices_per_cell,indices_per_cell);
  Full_matrix Mass(indices_per_cell,indices_per_cell);

  // matrix which stores Integrate phi * phi_x dx. Mass_y is the same but for the y-derivative and Mass is the normal mass
  // matrix.
  Mass_x = this->Compute_Mass_shape_grad_x(fe_v, indices_per_cell, J);
  Mass_y = this->Compute_Mass_shape_grad_y(fe_v, indices_per_cell, J);
  Mass = this->Compute_Mass_shape_value(fe_v, indices_per_cell, J);

  cell_matrix = matrix_opt.compute_A_outer_B(system_info->system_data.A[0].matrix,Mass_x) 
                    +  matrix_opt.compute_A_outer_B(system_info->system_data.A[1].matrix,Mass_y)
                    +  matrix_opt.compute_A_outer_B(system_info->system_data.P.matrix,Mass);

  std::vector<std::vector<double>> component_to_system(components_per_cell,std::vector<double> (indices_per_cell));

  AssertDimension(indices_per_cell,this->shape_values.rows());
  AssertDimension(total_ngp,this->shape_values.cols());

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
     // for the right hand side we iterate over all the components of the test function
     for (unsigned int m = 0 ; m < components_per_cell ; m++)
     {
       const int dof_test = component_to_system[m][index_test];
       const double shape_value_test = this->shape_values(index_test,q);

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

    boundary_rhs_value.reinit(this->constants.nBC);

for (unsigned int q = 0 ; q < fe_v.n_quadrature_points ; q++)
{
  const double jacobian_value = J[q];

  boundary_rhs_value = 0;                 


    system_info->build_BCrhs(fe_v.quadrature_point(q),fe_v.normal_vector(q),
                          boundary_rhs_value,b_id);


  Sparse_matrix Projector = system_info->build_Projector(fe_v.normal_vector(q));
  


  // Simga(as given in PDF) * B * Projector
  // Sigma in the PDF = Projector.transpose * Sigma(In the code)
  Full_matrix Sigma_B_P =       Projector.transpose()
                               * system_info->system_data.Sigma.matrix 
                               * system_info->system_data.B.matrix
                               * Projector;

  Full_matrix Sigma = Projector.transpose()
                      * system_info->system_data.Sigma.matrix ;

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
Base_Solver<dim>
::integrate_face_manuel(Full_matrix &u1_v1,
                        Full_matrix &u1_v2,
                        Full_matrix &u2_v1,
                        Full_matrix &u2_v2,
                        FEValuesBase<dim> &fe_v,
                        FEValuesBase<dim> &fe_v_neighbor,
                        std::vector<double> &J,
                        Vector<double> &component,
                        const typename DoFHandler<dim>::active_cell_iterator &cell)
{

    Assert(u1_v1.rows() !=0 || u1_v1.cols() !=0 ,ExcNotInitialized());
    Assert(u1_v2.rows() !=0 || u1_v2.cols() !=0 ,ExcNotInitialized());
    Assert(u2_v1.rows() !=0 || u2_v1.cols() !=0 ,ExcNotInitialized());
    Assert(u2_v2.rows() !=0 || u2_v2.cols() !=0 ,ExcNotInitialized());
    Assert(J.size() !=0,ExcNotInitialized());
    
    const FiniteElement<dim> &fe_in_cell1 = fe_v.get_fe();
    const FiniteElement<dim> &fe_in_cell2 = fe_v_neighbor.get_fe();

    const unsigned int dofs_per_cell1 = fe_in_cell1.dofs_per_cell;
    const unsigned int n_components1 = fe_in_cell1.n_components();
    const unsigned int dofs_per_component1 = dofs_per_cell1/n_components1;

    const unsigned int dofs_per_cell2 = fe_in_cell2.dofs_per_cell;
    const unsigned int n_components2 = fe_in_cell2.n_components();
    const unsigned int dofs_per_component2 = dofs_per_cell2/n_components2;


    //CAUTION: ASSUMPTION OF STRAIGHT EDGES IN THE INTERIOR
    Tensor<1,dim> outward_normal = fe_v.normal_vector(0);
    Eigen::MatrixXd Am = system_info->build_Aminus(outward_normal);
    Eigen::MatrixXd Am_neighbor = system_info->build_Aminus(-outward_normal);


    Full_matrix Mass_u1_v1(dofs_per_component1,dofs_per_component1);
    Full_matrix Mass_u2_v1(dofs_per_component1,dofs_per_component2);
    Full_matrix Mass_u2_v2(dofs_per_component2,dofs_per_component2);
    Full_matrix Mass_u1_v2(dofs_per_component2,dofs_per_component1);

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

    u1_v1 = 0.5 * matrix_opt.compute_A_outer_B(Am,Mass_u1_v1);
    u2_v1 = -0.5 * matrix_opt.compute_A_outer_B(Am,Mass_u2_v1);
    u1_v2 = -0.5 * matrix_opt.compute_A_outer_B(Am_neighbor,Mass_u1_v2);
    u2_v2 = 0.5 * matrix_opt.compute_A_outer_B(Am_neighbor,Mass_u2_v2);

}