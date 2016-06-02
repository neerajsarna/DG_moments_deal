template<int num_flux,int dim>
void 
Solver_DG<num_flux,dim>::assemble_rhs()
{

	const QGauss<dim> quadrature(ngp);
	const UpdateFlags update_flags  = update_values | update_JxW_values | update_quadrature_points;

	FEValues<dim>  fe_v(mapping,finite_element,quadrature, update_flags);
	typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();

	const unsigned int total_ngp = quadrature.size();
	const unsigned int dofs_per_cell = finite_element.dofs_per_cell;
    
    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
	 Vector<double> cell_rhs(dofs_per_cell);
	 vector<double> Jacobians_interior(total_ngp);
	 vector<Vector<double>> source_term_value(total_ngp,Vector<double>(nEqn));
    Vector<double> component(dofs_per_cell);

    for (unsigned int i = 0 ; i < dofs_per_cell ; i++)
        component[i] = finite_element.system_to_component_index(i).first;

	for (; cell != endc ; cell++) 
      {
      	cell_rhs = 0;
      	fe_v.reinit(cell);
        Jacobians_interior = fe_v.get_JxW_values();      	
        equation_system_data->source_term(fe_v.get_quadrature_points(),source_term_value,solve_system);

      	for (unsigned int q = 0 ; q < total_ngp ; q++)
      		for (unsigned int i = 0 ; i < dofs_per_cell ; i ++)
            cell_rhs(i) += fe_v.shape_value(i,q) * source_term_value[q][component[i]] * Jacobians_interior[q];

        cell->get_dof_indices(local_dof_indices);

        for(unsigned int i = 0 ; i < dofs_per_cell ; i++)
              system_rhs(local_dof_indices[i]) += cell_rhs(i);

      }
    
}

template<int num_flux,int dim>
void 
Solver_DG<num_flux,dim>
::assemble_system_meshworker()
{
	// first we will assemble the right hand side 
	Threads::Task<> rhs_task = Threads::new_task (&Solver_DG<num_flux,dim>::assemble_rhs,
                                                *this);


	MeshWorker::IntegrationInfoBox<dim> info_box;

	info_box.initialize_gauss_quadrature(ngp,
                                     	 ngp_face,
                                     	 ngp_face);


	info_box.initialize_update_flags();
	UpdateFlags update_flags = update_quadrature_points |
    	                       update_values            |
        	                   update_gradients;

	info_box.add_update_flags(update_flags, true, true, true, true);
	info_box.initialize(finite_element, mapping);
	MeshWorker::DoFInfo<dim> dof_info(dof_handler);

	MeshWorker::Assembler::SystemSimple<TrilinosWrappers::SparseMatrix, Vector<double>> assembler;
	assembler.initialize(global_matrix, system_rhs);


	typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();
	typename DoFHandler<dim>::active_cell_iterator endc = dof_handler.end();


	MeshWorker::loop<dim, dim, MeshWorker::DoFInfo<dim>, MeshWorker::IntegrationInfoBox<dim> >
  												(cell, endc,
   												 dof_info, info_box,
   												 std_cxx11::bind(&Solver_DG<num_flux,dim>::integrate_cell_term,
														this,
														std_cxx11::_1,std_cxx11::_2),
   												 std_cxx11::bind(&Solver_DG<num_flux,dim>::integrate_boundary_term,
   												 	this,
   												 	std_cxx11::_1,
   												 	std_cxx11::_2),
   												 std_cxx11::bind(&Solver_DG<num_flux,dim>::integrate_face_term,
   												 	this,
   												 	std_cxx11::_1,
   												 	std_cxx11::_2,
   												 	std_cxx11::_3,
   												 	std_cxx11::_4),
   												 assembler);
     rhs_task.join();
}


template<int num_flux,int dim>
 void 
 Solver_DG<num_flux,dim>::integrate_cell_term (DoFInfo &dinfo,
                                   			CellInfo &info)
{
  
	const FEValuesBase<dim> &fe_v = info.fe_values();
	FullMatrix<double> &cell_matrix = dinfo.matrix(0).matrix;
	const std::vector<double> &Jacobians_interior = fe_v.get_JxW_values ();
	const FiniteElement<dim> &fe_in_cell = info.finite_element();

	const unsigned int dofs_per_cell = fe_in_cell.dofs_per_cell;
	const unsigned int components_per_cell  = fe_in_cell.n_components();
	const unsigned int indices_per_cell = dofs_per_cell/components_per_cell;
	vector<vector<double>> component_to_system(components_per_cell,vector<double> (indices_per_cell));

	for (unsigned int i = 0 ; i < components_per_cell ; i ++)
		for (unsigned int j = 0 ; j < indices_per_cell ; j ++)
			component_to_system[i][j] = fe_in_cell.component_to_system_index(i,j); 
	

	for (unsigned int q = 0 ; q < fe_v.n_quadrature_points ; q++)
	{
		double jacobian_value = Jacobians_interior[q];
		
		for (unsigned int index_test = 0 ; index_test < indices_per_cell ; index_test++)
		  for (unsigned int index_sol = 0 ; index_sol < indices_per_cell ; index_sol++ )
		  {
			for (unsigned int space = 0 ; space < dim ; space++)
		   		for (unsigned int m = 0 ; m < equation_system_data->system_data[solve_system].A[space].matrix.outerSize(); m++)
				{
					const int dof_test = component_to_system[m][index_test];
					const double shape_value_test = fe_v.shape_value(dof_test,q);

					for (Sparse_matrix::InnerIterator n(equation_system_data->system_data[solve_system].A[space].matrix,m); n ; ++n)
	
					{
						int dof_sol = component_to_system[n.col()][index_sol];
						double grad_value_sol = fe_v.shape_grad(dof_sol,q)[space];

						cell_matrix(dof_test,dof_sol) += shape_value_test * n.value()
								* grad_value_sol * jacobian_value;	
					}
				}
			for (unsigned int m = 0 ; m < equation_system_data->system_data[solve_system].P.matrix.outerSize(); m++)
			{
				const int dof_test = component_to_system[m][index_test];
				const double shape_value_test = fe_v.shape_value(dof_test,q);

				for (Sparse_matrix::InnerIterator n(equation_system_data->system_data[solve_system].P.matrix,m); n ; ++n)
					{
						const int dof_sol = component_to_system[n.col()][index_sol];
						const double shape_value_sol = fe_v.shape_value(dof_sol,q);

						cell_matrix(dof_test,dof_sol) += shape_value_test * n.value()
								* shape_value_sol * jacobian_value;	
					}
			}

		 }
		
	}


}

template<int num_flux,int dim> 
void 
Solver_DG<num_flux,dim>
::integrate_boundary_term(DoFInfo &dinfo,
                          CellInfo &info)
{
 const FEValuesBase<dim> &fe_v = info.fe_values();
 FullMatrix<double> &cell_matrix = dinfo.matrix(0).matrix;
 Vector<double> &cell_rhs = dinfo.vector(0).block(0);
 const std::vector<double> &Jacobian_face = fe_v.get_JxW_values ();

 const FiniteElement<dim> &fe_in_cell = info.finite_element();

 const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
 std::vector<unsigned int> component(dofs_per_cell);

 for (unsigned int i = 0 ; i < dofs_per_cell ; i ++)
  component[i] = fe_in_cell.system_to_component_index(i).first;

  	  for (unsigned int q = 0 ; q < fe_v.n_quadrature_points ; q++)
                {
		              const double jacobian_value = Jacobian_face[q];
                  Vector<double> boundary_rhs_value(this->nEqn);
                  equation_system_data->build_BCrhs(fe_v.quadrature_point(q),fe_v.normal_vector(q),
                                                  boundary_rhs_value,solve_system);
                  Eigen::MatrixXd Am = equation_system_data->build_Aminus(fe_v.normal_vector(q),solve_system);
                  Sparse_matrix Projector = equation_system_data->build_Projector(fe_v.normal_vector(q),solve_system);
                  Sparse_matrix Inv_Projector = equation_system_data->build_InvProjector(fe_v.normal_vector(q),solve_system);

                  Eigen::MatrixXd Am_invP_BC_P = Am * Inv_Projector 
                                                * equation_system_data->system_data[solve_system].BC.matrix * Projector;

<<<<<<< HEAD
                  Eigen::MatrixXd Am_invP = Am * Inv_Projector;
=======
        Eigen::MatrixXd Am_invP_B_hat_P = Am * Inv_Projector 
                                      * equation_system_data->system_data[solve_system].B_hat 
<<<<<<< HEAD
=======
                                      * Projector;

        Eigen::MatrixXd Am_invP_X_min_B_tild_inv = Am
                                                   * Inv_Projector 
                                                   * equation_system_data->system_data[solve_system].X_minus
                                                   * equation_system_data->system_data[solve_system].B_tilde_inv;


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
                          * Am_invP_X_min_B_tild_inv(component[i],j) * boundary_rhs_value[j] * jacobian_value;

        }

      break;
    }

    case odd:
    {
        equation_system_data->build_BCrhs(fe_v.quadrature_point(q),fe_v.normal_vector(q),
                        boundary_rhs_value,solve_system);

        // build the matrices needed
        Eigen::MatrixXd Am = equation_system_data->build_Aminus(fe_v.normal_vector(q),solve_system);
        Sparse_matrix Projector = equation_system_data->build_Projector(fe_v.normal_vector(q),solve_system);
        Sparse_matrix Inv_Projector = equation_system_data->build_InvProjector(fe_v.normal_vector(q),solve_system);

        Eigen::MatrixXd Am_invP_BC_P = Am * Inv_Projector 
                                      * equation_system_data->system_data[solve_system].BC.matrix 
>>>>>>> parent of 566d0df... updated assemblation
                                      * Projector;
>>>>>>> parent of 566d0df... updated assemblation


                  for (unsigned int i = 0 ; i < dofs_per_cell ; i ++)
                  {
		                const double shape_value_test = fe_v.shape_value(i,q);
                    for (unsigned int j = 0 ; j < dofs_per_cell ; j ++)
                        cell_matrix(i,j) += 0.5 * shape_value_test
                                           * (Am(component[i],component[j])-Am_invP_BC_P(component[i],component[j]))
                                            * fe_v.shape_value(j,q) * jacobian_value;                                    
                      

                    for (unsigned int j = 0 ; j < Am_invP.cols() ; j++)
                     cell_rhs(i) -= 0.5 * shape_value_test 
                                    * Am_invP(component[i],j) * boundary_rhs_value[j] * jacobian_value;


                  }

                }

}

template<int num_flux,int dim>
void 
Solver_DG<num_flux,dim>
::integrate_face_term(DoFInfo &dinfo1,
                      DoFInfo &dinfo2,
                      CellInfo &info1,
                      CellInfo &info2)
{
  
	const FEValuesBase<dim> &fe_v = info1.fe_values();
	const FEValuesBase<dim> &fe_v_neighbor = info2.fe_values();

  FullMatrix<double> &u1_v1 = dinfo1.matrix(0,false).matrix;
  FullMatrix<double> &u2_v1 = dinfo1.matrix(0,true).matrix;
  FullMatrix<double> &u1_v2 = dinfo2.matrix(0,true).matrix;
  FullMatrix<double> &u2_v2 = dinfo2.matrix(0,false).matrix;

  const std::vector<double> &Jacobian_face = fe_v.get_JxW_values ();

  const FiniteElement<dim> &fe_in_cell = info1.finite_element();
  const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
  std::vector<unsigned int> component(dofs_per_cell);

  for (unsigned int i = 0 ; i < dofs_per_cell ; i ++)
    component[i] = fe_in_cell.system_to_component_index(i).first;

  for (unsigned int q = 0 ; q < fe_v.n_quadrature_points ; q++)
  {
                  // build the matrices needed
    Eigen::MatrixXd Am = equation_system_data->build_Aminus(fe_v.normal_vector(q),solve_system);
    Eigen::MatrixXd Am_neighbor = equation_system_data->build_Aminus(-fe_v.normal_vector(q),solve_system);
    double jacobian_value = Jacobian_face[q];

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
