template<int dim>
void 
Base_Solver<dim>
::integrate_cell_manuel(FullMatrix<double> &cell_matrix,
                        Vector<double> &cell_rhs,
                        const FEValuesBase<dim> &fe_v,
                        std::vector<double> &J,
                        std::vector<Vector<double>> &source_term_value,
                        const typename hp::DoFHandler<dim>::active_cell_iterator &cell,
                        std::vector<Vector<double>> &component_to_system,
                        const unsigned int fe_index)
{

  Assert(cell_matrix.m() !=0 || cell_matrix.n() !=0 ,ExcNotInitialized());
  Assert(cell_matrix.m() == cell_matrix.n() ,ExcMessage("different dof in same cell"));
  Assert(cell_rhs.size() != 0,ExcNotInitialized());
  Assert(J.size() !=0,ExcNotInitialized());
  Assert(source_term_value.size() != 0,ExcNotInitialized());
  Assert(fe_index < nEqn.size(),ExcMessage("Incorrect fe index"));
  AssertDimension(component_to_system.size(),nEqn[fe_index]);
  AssertDimension(cell_matrix.m(),dofs_per_cell[fe_index]);
  AssertDimension(source_term_value.size(),J.size());
  AssertDimension(fe_v.get_fe().dofs_per_cell,dofs_per_cell[fe_index]);


  const unsigned int total_ngp = J.size();


  AssertDimension(dofs_per_component,this->shape_values.rows());
  AssertDimension(total_ngp,this->shape_values.cols());

  std::vector<FullMatrix<double>> Mass_grad(dim);
  FullMatrix<double> Mass(dofs_per_component,dofs_per_component);

  const int components_per_cell = nEqn[fe_index];
  Mass_grad = this->Compute_Mass_shape_grad(fe_v, dofs_per_component, J);
  Mass = this->Compute_Mass_shape_value(fe_v, dofs_per_component, J);

  // this is for initializing the matrix, it is necessary since we have not put cell_matrix to zero
  cell_matrix = matrix_opt.compute_A_outer_B(system_info[fe_index].system_data.P.matrix,Mass);

  for (int space = 0 ;space < dim ; space ++ )
        cell_matrix.add(0,cell_matrix,1,matrix_opt.compute_A_outer_B(system_info[fe_index].system_data.A[space].matrix,Mass_grad[space]));
                    

  AssertDimension(dofs_per_component,this->shape_values.rows());
  AssertDimension(total_ngp,this->shape_values.cols());

  // loop over all the gauss points
  for (unsigned int q = 0 ; q < total_ngp ; q++)
  {
    // value of the jacobian at the gauss point
    double jacobian_value = J[q];
    
    // index is the local id of the degree of freedom
    for (unsigned int index_test = 0 ; index_test < dofs_per_component ; index_test++)
    {
     // for the right hand side we iterate over all the components of the test function
     for (int m = 0 ; m < nEqn[fe_index] ; m++)
     {
       const int dof_test = component_to_system[m](index_test);
       const double shape_value_test = this->shape_values(index_test,q);

       cell_rhs(dof_test) += shape_value_test * source_term_value[q][m] * jacobian_value;
     }
   }
    
  }

  // std::string filename = "cell_matrix_internal" + std::to_string(fe_index);

  // matrix_opt.print_dealii_full(cell_matrix,filename);


}


// characteristic implementation of the boundary conditions
// template<int dim>
// void 
// Base_Solver<dim>
// ::integrate_boundary_manuel_char(FullMatrix<double> &cell_matrix,
//                             Vector<double> &cell_rhs,
//                             FEValuesBase<dim> &fe_v,
//                             std::vector<double> &J,
//                             Vector<double> &component,
//                             const typename DoFHandler<dim>::active_cell_iterator &cell,
//                             const unsigned int b_id)
// {
//     AssertThrow(1 == 0,ExcMessage("The present routine has not been further developed. Use odd boundary implementation."));
//     Assert(cell_matrix.m() !=0 || cell_matrix.n() !=0 ,ExcNotInitialized());
//     Assert(cell_matrix.m() == cell_matrix.n() ,ExcMessage("different dof in same cell"));
//     Assert(cell_rhs.size() != 0,ExcNotInitialized());
//     Assert(J.size() !=0,ExcNotInitialized());
    
//     const FiniteElement<dim> &fe_in_cell = cell->get_fe();
//    const unsigned int dofs_per_cell = fe_in_cell.dofs_per_cell;

//      Vector<double> boundary_rhs_value;

//       boundary_rhs_value.reinit(nBC[0]);

//  for (unsigned int q = 0 ; q < fe_v.n_quadrature_points ; q++)
//  {
//   const double jacobian_value = J[q];

//   boundary_rhs_value = 0;                 

//          // the outward normal to the boundary
//   Tensor<1,dim> outward_normal = fe_v.normal_vector(q);

//   system_info[0].bcrhs_wall.BCrhs(fe_v.quadrature_point(q),outward_normal,
//                                     boundary_rhs_value,b_id);

//         // build the matrices needed
//   Full_matrix Am = system_info[0].build_Aminus(outward_normal);
//   Sparse_matrix Projector = system_info[0].build_Projector(outward_normal);
//   Sparse_matrix Inv_Projector = system_info[0].build_InvProjector(outward_normal);

//   Eigen::MatrixXd Am_invP_B_hat_P = Am * Inv_Projector 
//                                     * system_info[0].B_hat 
//                                     * Projector;

//   Eigen::MatrixXd Am_invP_X_min_B_tild_inv = Am
//                                             * Inv_Projector 
//                                             * system_info[0].X_minus
//                                             * system_info[0].B_tilde_inv;


//   for (unsigned int i = 0 ; i < dofs_per_cell ; i ++)
//   {
//     const double shape_value_test = fe_v.shape_value(i,q);
//     for (unsigned int j = 0 ; j < dofs_per_cell ; j ++)
//       cell_matrix(i,j) += 0.5 * shape_value_test
//                         * Am_invP_B_hat_P(component[i],component[j])
//                         * fe_v.shape_value(j,q) 
//                         * jacobian_value;                                    


//     for (unsigned int j = 0 ; j < Am_invP_X_min_B_tild_inv.cols() ; j++)
//      cell_rhs(i) += 0.5 * shape_value_test 
//                     * Am_invP_X_min_B_tild_inv(component[i],j)
//                      * boundary_rhs_value[j] * jacobian_value;

//  }

//  }
// }


// // odd implementation of the boundary conditions

// template<int dim>
// void 
// Base_Solver<dim>
// ::integrate_boundary_manuel_odd(FullMatrix<double> &cell_matrix,
//                             Vector<double> &cell_rhs,
//                             FEValuesBase<dim> &fe_v,
//                             std::vector<double> &J,
//                             Vector<double> &component,
//                             const typename DoFHandler<dim>::active_cell_iterator &cell,
//                             const unsigned int b_id)
// {


//     Assert(cell_matrix.m() !=0 || cell_matrix.n() !=0 ,ExcNotInitialized());
//     Assert(cell_matrix.m() == cell_matrix.n() ,ExcMessage("different dof in same cell"));
//     Assert(cell_rhs.size() != 0,ExcNotInitialized());
//     Assert(J.size() !=0,ExcNotInitialized());
    
//     const FiniteElement<dim> &fe_in_cell = cell->get_fe();
//     const unsigned int dofs_per_cell = fe_in_cell.dofs_per_cell;

//     Vector<double> boundary_rhs_value;

//     boundary_rhs_value.reinit(nBC[0]);

//     // we use a temporary matrix to determine whether inflow or outflow
//     Sparse_matrix B_temp;

//     if(b_id == 101 || b_id == 102)
//     {
//       integrate_inflow++;
//       B_temp = system_info[0].system_data.Binflow.matrix;
//     }
//     else
//     {
//     // B matrix for specular reflection
//       if (b_id == 50)
//         B_temp = system_info[0].system_data.B_specular.matrix;

//     // B matrix for full accmmodation
//       else
//         B_temp = system_info[0].system_data.B.matrix;
//     }

// for (unsigned int q = 0 ; q < fe_v.n_quadrature_points ; q++)
// {
//   const double jacobian_value = J[q];

//   boundary_rhs_value = 0;                 


//   // check for inflow or outflow
//   // Incase of inflow provide the inflow rhs
//   if(b_id == 101 || b_id == 102)
//     system_info[0].bcrhs_inflow.BCrhs(fe_v.quadrature_point(q),fe_v.normal_vector(q),
//                           boundary_rhs_value,b_id);

//   else
//     system_info[0].bcrhs_wall.BCrhs(fe_v.quadrature_point(q),fe_v.normal_vector(q),
//                           boundary_rhs_value,b_id);


//   Sparse_matrix Projector = system_info[0].build_Projector(fe_v.normal_vector(q));
  


//   // Simga(as given in PDF) * B * Projector
//   // Sigma in the PDF = Projector.transpose * Sigma(In the code)
//   Full_matrix Sigma_B_P =       Projector.transpose()
//                                * system_info[0].system_data.Sigma.matrix 
//                                * system_info[0].system_data.B.matrix
//                                * Projector;

//   Full_matrix Sigma = Projector.transpose()
//                       * system_info[0].system_data.Sigma.matrix ;

//   for (unsigned int i = 0 ; i < dofs_per_cell ; i ++)
//   {
//     const double shape_value_test = fe_v.shape_value(i,q);
//     for (unsigned int j = 0 ; j < dofs_per_cell ; j ++)
//       cell_matrix(i,j) -= shape_value_test
//                           * Sigma_B_P(component[i],component[j])
//                           * fe_v.shape_value(j,q) 
//                           * jacobian_value;                                    


//     for (unsigned int j = 0 ; j < boundary_rhs_value.size() ; j++)
//      cell_rhs(i) -= shape_value_test 
//                     * Sigma(component[i],j) 
//                     * boundary_rhs_value[j] 
//                     * jacobian_value;

//  }


// }
// }


// template<int dim>
// void 
// Base_Solver<dim>
// ::integrate_face_manuel(FullMatrix<double> &u1_v1,
//                         FullMatrix<double> &u1_v2,
//                         FullMatrix<double> &u2_v1,
//                         FullMatrix<double> &u2_v2,
//                         FEValuesBase<dim> &fe_v,
//                         FEValuesBase<dim> &fe_v_neighbor,
//                         std::vector<double> &J,
//                         Vector<double> &component,
//                         const typename DoFHandler<dim>::active_cell_iterator &cell)
// {

//     Assert(u1_v1.m() !=0 || u1_v1.n() !=0 ,ExcNotInitialized());
//     Assert(u1_v2.m() !=0 || u1_v2.n() !=0 ,ExcNotInitialized());
//     Assert(u2_v1.m() !=0 || u2_v1.n() !=0 ,ExcNotInitialized());
//     Assert(u2_v2.m() !=0 || u2_v2.n() !=0 ,ExcNotInitialized());
//     Assert(J.size() !=0,ExcNotInitialized());
    
//     const FiniteElement<dim> &fe_in_cell1 = fe_v.get_fe();
//     const FiniteElement<dim> &fe_in_cell2 = fe_v_neighbor.get_fe();

//     const unsigned int dofs_per_cell1 = fe_in_cell1.dofs_per_cell;
//     const unsigned int n_components1 = fe_in_cell1.n_components();
//     const unsigned int dofs_per_component1 = dofs_per_cell1/n_components1;

//     const unsigned int dofs_per_cell2 = fe_in_cell2.dofs_per_cell;
//     const unsigned int n_components2 = fe_in_cell2.n_components();
//     const unsigned int dofs_per_component2 = dofs_per_cell2/n_components2;


//     //CAUTION: ASSUMPTION OF STRAIGHT EDGES IN THE INTERIOR
//     Tensor<1,dim> outward_normal = fe_v.normal_vector(0);
//     Eigen::MatrixXd Am = system_info[0].build_Aminus(outward_normal);
//     Eigen::MatrixXd Am_neighbor = system_info[0].build_Aminus(-outward_normal);


//     FullMatrix<double> Mass_u1_v1(dofs_per_component1,dofs_per_component1);
//     FullMatrix<double> Mass_u2_v1(dofs_per_component1,dofs_per_component2);
//     FullMatrix<double> Mass_u2_v2(dofs_per_component2,dofs_per_component2);
//     FullMatrix<double> Mass_u1_v2(dofs_per_component2,dofs_per_component1);

//     Mass_u1_v1 = this->Compute_Mass_cell_neighbor(fe_v,
//                                                    fe_v,dofs_per_component1,
//                                                   dofs_per_component1,J);

//     Mass_u2_v1 = this->Compute_Mass_cell_neighbor(fe_v,
//                                                    fe_v_neighbor,dofs_per_component1,
//                                                   dofs_per_component2,J);

//     Mass_u2_v2 =  this->Compute_Mass_cell_neighbor(fe_v_neighbor,
//                                                    fe_v_neighbor,dofs_per_component2,
//                                                   dofs_per_component2,J);

//     Mass_u1_v2 =  this->Compute_Mass_cell_neighbor(fe_v_neighbor,
//                                                    fe_v,dofs_per_component2,
//                                                   dofs_per_component1,J);

//     u1_v1 = matrix_opt.multiply_scalar(0.5,matrix_opt.compute_A_outer_B(Am,Mass_u1_v1));
//     u2_v1 = matrix_opt.multiply_scalar(-0.5,matrix_opt.compute_A_outer_B(Am,Mass_u2_v1));
//     u1_v2 = matrix_opt.multiply_scalar(-0.5,matrix_opt.compute_A_outer_B(Am_neighbor,Mass_u1_v2));
//     u2_v2 = matrix_opt.multiply_scalar(0.5,matrix_opt.compute_A_outer_B(Am_neighbor,Mass_u2_v2));

// }