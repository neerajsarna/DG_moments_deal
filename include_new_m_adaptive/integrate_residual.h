template<int dim>
void 
Base_Solver<dim>
::residual_cell_manuel(Vector<double> &residual,
                         std::vector<types::global_dof_index> &local_dof_indices,
                        const FEValuesBase<dim> &fe_v,
                        std::vector<double> &J,
                        std::vector<Vector<double>> &source_term_value,
                        std::vector<Vector<double>> &component_to_system,
                        const unsigned int fe_index)
{

  Assert(source_term_value.size() != 0,ExcNotInitialized());
  Assert(fe_index < nEqn.size(),ExcMessage("Incorrect fe index"));
  AssertDimension((int)component_to_system.size(),nEqn[fe_index]);
  AssertDimension(source_term_value.size(),J.size());
  AssertDimension(fe_v.get_fe().dofs_per_cell,dofs_per_cell[fe_index]);


  const unsigned int total_ngp = J.size();

  // number of equations in the present cell
  const int components_per_cell = nEqn[fe_index];

  // value of the gradient of the solution
  std::vector<Tensor<1,dim>> value_gradient(components_per_cell);

  // we need to take a transpose of the above array so as to perform the matrix vector dot product
  std::vector<Vector<double>> temp_value_gradient(dim,Vector<double>(components_per_cell));

  // the flux matrix dotted with the value of the gradient of the solution
  Vector<double> A_dot_gradient(components_per_cell);

  for (unsigned int q = 0 ; q < total_ngp ; q++)
  {
    // value of the gradient of the solution for different dimensions
    VectorTools::point_gradient ( dof_handler_max_moments,
                                  solution_max_moments,
                                  fe_v.quadrature_point(q),
                                  value_gradient);

    const double jacobian_value = J[q];

    // now we take the tranpose of the value_gradients matrix. Such a structure would be preferrable 
    // to take the dot product
    for (int eq = 0 ; eq < components_per_cell ; eq++)
      for (unsigned int space = 0 ; space < dim ; space ++)
        temp_value_gradient[space][eq] = value_gradient[eq][space];

      A_dot_gradient = 0;
    // we now compute the flux matrix times the gradient of the solution
    for (unsigned int space = 0 ; space < dim ; space ++)
      A_dot_gradient += matrix_opt.Sparse_matrix_dot_Vector(system_info[fe_index].system_data.A[space].matrix,temp_value_gradient[space]);

   // We now compute the residual generated from the integral in the interior of the cell.
    // We first loop over the degrees of freedom of a particular component.
   for (unsigned int index_test = 0 ; index_test < dofs_per_component ; index_test++)
    // now we loop over all the equations present in the system
     for (int m = 0 ; m < nEqn[fe_index] ; m++)
     {
            // computation of the local degree of freedom corresponding to a particular component and dof per component
             const int dof_test = component_to_system[m](index_test);

             // we bring the right hand side to the left and perform the integration. 
             residual(local_dof_indices[dof_test])  += fe_v.shape_value(dof_test,q) * ( A_dot_gradient(m)
                                                      -source_term_value[q][m] ) *  jacobian_value;
     }
      
  }

}

template<int dim>
void 
Base_Solver<dim>
::residual_boundary_manuel_odd(Vector<double> &residual,
                         std::vector<types::global_dof_index> &local_dof_indices,
                            const FEValuesBase<dim> &fe_v,
                            std::vector<double> &J,
                            std::vector<Vector<double>> &component_to_system,
                            const unsigned int b_id,
                            const unsigned int fe_index)
{

    Vector<double> boundary_rhs_value;

    boundary_rhs_value.reinit(nBC[fe_index]);

    // we use a temporary matrix to determine whether inflow or outflow
    Sparse_matrix B_temp;

    if(b_id == 101 || b_id == 102)
    {
      integrate_inflow++;
      B_temp = system_info[fe_index].system_data.Binflow.matrix;
    }
    else
    {
    // B matrix for specular reflection
      if (b_id == 50)
        B_temp = system_info[fe_index].system_data.B_specular.matrix;

    // B matrix for full accmmodation
      else
        B_temp = system_info[fe_index].system_data.B.matrix;
    }

    Vector<double> Sigma_B_P_dot_sol(nEqn[fe_index]);
    Vector<double> Sigma_dot_bcrhs(nEqn[fe_index]);

    // stores the value of the solution at the quadrature points
    Vector<double> solution_value(nEqn[fe_index]);

for (unsigned int q = 0 ; q < fe_v.n_quadrature_points ; q++)
{
  const double jacobian_value = J[q];

  boundary_rhs_value = 0;                 

  // value of the gradient of the solution for different dimensions
  VectorTools::point_value ( dof_handler_max_moments,
                                  solution_max_moments,
                                  fe_v.quadrature_point(q),
                                  solution_value);


  // check for inflow or outflow
  // Incase of inflow provide the inflow rhs
  if(b_id == 101 || b_id == 102)
    system_info[fe_index].bcrhs_inflow.BCrhs(fe_v.quadrature_point(q),fe_v.normal_vector(q),
                          boundary_rhs_value,b_id);

  else
    system_info[fe_index].bcrhs_wall.BCrhs(fe_v.quadrature_point(q),fe_v.normal_vector(q),
                          boundary_rhs_value,b_id);         


  Sparse_matrix Projector = system_info[fe_index].build_Projector(fe_v.normal_vector(q));
  


  // Simga(as given in PDF) * B * Projector
  // Sigma in the PDF = Projector.transpose * Sigma(In the code)
  Full_matrix Sigma_B_P =       Projector.transpose()
                               * system_info[fe_index].system_data.Sigma.matrix 
                               * B_temp
                               * Projector;

  Full_matrix Sigma = Projector.transpose()
                      * system_info[fe_index].system_data.Sigma.matrix ;


   Sigma_B_P_dot_sol = matrix_opt.Sparse_matrix_dot_Vector(Sigma_B_P,solution_value);
  // Sigma_dot_bcrhs = matrix_opt.Sparse_matrix_dot_Vector(Sigma,boundary_rhs_value);


  // we first loop over all the components of the test function
  for (int component_test = 0; component_test < nEqn[fe_index] ; component_test++)
    // we now loop over all the dof for our component
    for (unsigned long int index_test = 0 ; index_test < dofs_per_component ; index_test ++)
    {
      const unsigned long int dof_test = component_to_system[component_test][index_test];
      const double shape_value_test = fe_v.shape_value(dof_test,q);

      // we bring the right hand side from the boundary conditions to the left and add to the residual.
      residual(local_dof_indices[dof_test]) += shape_value_test * (
                                              -Sigma_B_P_dot_sol(component_test) 
                                              +
                                               Sigma_dot_bcrhs(component_test)) 
                                              * jacobian_value; 

    }





}

}

// compute the residual arising from the face integrals
template<int dim>
void 
Base_Solver<dim>
::residual_face_manuel( Vector<double> &residual,
                         std::vector<types::global_dof_index> &local_dof_indices,
                         std::vector<types::global_dof_index> &local_dof_indices_neighbor,
                        const FEValuesBase<dim> &fe_v,
                        const FEValuesBase<dim> &fe_v_neighbor,
                        std::vector<double> &J,
                        const unsigned int fe_index1,
                        const unsigned int fe_index2,
                        std::vector<Vector<double>> &component_to_system)
{

    //CAUTION: ASSUMPTION OF STRAIGHT EDGES IN THE INTERIOR
    Tensor<1,dim> outward_normal = fe_v.normal_vector(0);
    const unsigned int max_fe_index = std::max(fe_index1,fe_index2);
    Full_matrix Am = system_info[max_fe_index].build_Aminus(outward_normal);
    Full_matrix Am_neighbor = system_info[max_fe_index].build_Aminus(-outward_normal);

    // number of equations in the present cell
    const unsigned int nEqn_this = nEqn[fe_index1];

    // number of equations in the neighbouring cell
    const unsigned int nEqn_neighbor = nEqn[fe_index2];

    // number of dofs in the present cell
    const unsigned int dofs_this = dofs_per_cell[fe_index1];

    // number of dofs in the neighbouring cell
    const unsigned int dofs_neighbor = dofs_per_cell[fe_index2];

    // since we are computing for the maximum possible moment system, therefore 
    // the number of equations and the number of degrees of freedom from both the sides should be the same.
    AssertDimension(nEqn_this,nEqn_neighbor);
    AssertDimension(dofs_this,dofs_neighbor);

    // value of the solution from this cell
    Vector<double> solution_value(nEqn_this);

    // value of the solution from the side of the neighbor. The build in function of the library i.e. the point_value
    // function because for discontinous functions it  takes the average of the value of the solution from neighbouring cells
    // therefore we will loose the contribution of the face terms into the residual. 
    Vector<double> solution_value_neighbor(nEqn_neighbor);

    // value of the dofs for this cell
    Vector<double> dof_values(dofs_per_cell[fe_index1]);

    // value of the dofs for the neighbor of this cell 
    Vector<double> dof_values_neighbor(dofs_per_cell[fe_index2]);


    // Extract the degrees of freedom for the present cell and the neighboring cell.
    for (unsigned long int i = 0 ; i < local_dof_indices.size(); i ++)
      dof_values(i) = solution_max_moments(local_dof_indices[i]);

    for (unsigned long int i = 0 ; i < local_dof_indices_neighbor.size(); i ++)
      dof_values_neighbor(i) = solution_max_moments(local_dof_indices_neighbor[i]);


    const unsigned int total_ngp = J.size();

      // Dot product between the numerical flux matrix and the solution values. The number of equations will remain the 
    // same since we are computing for the maximum number of moments. 
      Vector<double> Am_dot_solution(nEqn_this);
      Vector<double> Am_dot_solution_neighbor(nEqn_this);
      Vector<double> Am_neighbor_dot_solution(nEqn_this);
      Vector<double> Am_neighbor_dot_solution_neighbor(nEqn_this);

    for (unsigned int q = 0 ; q < total_ngp ; q++)
    {
      solution_value = 0;
      solution_value_neighbor = 0;
      const double jacobian_value = J[q];


      // Dot product between the numerical flux matrix and the solution values
      Am_dot_solution = 0;
      Am_dot_solution_neighbor = 0;
      Am_neighbor_dot_solution = 0;
      Am_neighbor_dot_solution_neighbor = 0;



      // the number of equations for this cell and the neighbor are the same therefore we loop only once
      for (unsigned int eq = 0 ; eq < nEqn_this ; eq++)
        for (unsigned int m = 0 ; m < dofs_per_component ; m++)
        {

          // the component to system is the same for this cell and the neighbouring cell. 
          // This is the dof id corresponding to a particular component and its corresponding dof
          const unsigned int dof_id = component_to_system[eq][m];
  

          // Solution = Sum(dofs * basis_functions). Value of the solution at the quadrature point = 
          // value of the dofs * value of the basis functions at the quadrature points.
          solution_value(eq) += dof_values(dof_id) * fe_v.shape_value(dof_id,q);
          solution_value_neighbor(eq) += dof_values_neighbor(dof_id) * fe_v_neighbor.shape_value(dof_id,q); 
        }


      Am_dot_solution = matrix_opt.Sparse_matrix_dot_Vector(Am,solution_value);
      Am_dot_solution_neighbor = matrix_opt.Sparse_matrix_dot_Vector(Am,solution_value_neighbor);
      Am_neighbor_dot_solution = matrix_opt.Sparse_matrix_dot_Vector(Am_neighbor,solution_value);
      Am_neighbor_dot_solution_neighbor = matrix_opt.Sparse_matrix_dot_Vector(Am_neighbor,solution_value_neighbor);

      // now we scale them will half or minus half factors
      Am_dot_solution *= 0.5;
      Am_dot_solution_neighbor *= -0.5;
      Am_neighbor_dot_solution *= -0.5;
      Am_neighbor_dot_solution_neighbor *= 0.5;

      
      for (unsigned int eq = 0 ; eq < nEqn_this ; eq++)
        for (unsigned int m = 0 ; m < dofs_per_component ; m++)
        {
            // the dof index corresponding to a particular component and dof per component
             const unsigned int dof_id = component_to_system[eq][m];
            
             // residual of the present cell. The factors of plus and minus have already been included in the 
             // definition of the vectors describing the dot products
             residual(local_dof_indices[dof_id]) += fe_v.shape_value(dof_id,q) * (Am_dot_solution(eq) 
                                                    + Am_dot_solution_neighbor(eq)) * jacobian_value;

             // residual of the neighbouring cell
             residual(local_dof_indices_neighbor[dof_id]) += fe_v_neighbor.shape_value(dof_id,q) * (Am_neighbor_dot_solution(eq) 
                                                            + Am_neighbor_dot_solution_neighbor(eq)) * jacobian_value;


        }

    }


}
