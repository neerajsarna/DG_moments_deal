template<int dim> void Solver_DG<dim>::assemble_rhs()
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
	vector<Vector<double>> source_term_value(total_ngp,Vector<double>(this->nEqn));
    Vector<double> component(dofs_per_cell);

    for (unsigned int i = 0 ; i < dofs_per_cell ; i++)
        component[i] = finite_element.system_to_component_index(i).first;

	for (; cell != endc ; cell++) 
      {
      	cell_rhs = 0;
      	fe_v.reinit(cell);
        Jacobians_interior = fe_v.get_JxW_values();      	
        this->source_term(fe_v.get_quadrature_points(),source_term_value);

      	for (unsigned int q = 0 ; q < total_ngp ; q++)
      		for (unsigned int i = 0 ; i < dofs_per_cell ; i ++)
            cell_rhs(i) += fe_v.shape_value(i,q) * source_term_value[q][component[i]] * Jacobians_interior[q];

        cell->get_dof_indices(local_dof_indices);

        for(unsigned int i = 0 ; i < dofs_per_cell ; i++)
              system_rhs(local_dof_indices[i]) += cell_rhs(i);

      }
}

template<int dim> void Solver_DG<dim>::assemble_system_meshworker()
{
	// first we will assemble the right hand side 
	Threads::Task<> rhs_task = Threads::new_task (&Solver_DG<dim>::assemble_rhs,
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


	  MeshWorker::loop<dim, dim, MeshWorker::DoFInfo<dim>, MeshWorker::IntegrationInfoBox<dim> >
  												(dof_handler.begin_active(), dof_handler.end(),
   												 dof_info, info_box,
   												 std_cxx11::bind(&Solver_DG<dim>::integrate_cell_term,
   												 	this,
   												 	std_cxx11::_1,
   												 	std_cxx11::_2),
   												 std_cxx11::bind(&Solver_DG<dim>::integrate_boundary_term,
   												 	this,
   												 	std_cxx11::_1,
   												 	std_cxx11::_2),
   												 std_cxx11::bind(&Solver_DG<dim>::integrate_face_term,
   												 	this,
   												 	std_cxx11::_1,
   												 	std_cxx11::_2,
   												 	std_cxx11::_3,
   												 	std_cxx11::_4),
   												 assembler);

  	rhs_task.join();
}


template<int dim> void Solver_DG<dim>::integrate_cell_term (DoFInfo &dinfo,
                                   							CellInfo &info)
{
	const FEValuesBase<dim> &fe_v = info.fe_values();
	FullMatrix<double> &cell_matrix = dinfo.matrix(0).matrix;
	const std::vector<double> &Jacobians_interior = fe_v.get_JxW_values ();
	const FiniteElement<dim> &fe_in_cell = info.finite_element();

	const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
	std::vector<unsigned int> component(dofs_per_cell);

	for (unsigned int i = 0 ; i < dofs_per_cell ; i ++)
		component[i] = fe_in_cell.system_to_component_index(i).first;

	for(unsigned int q = 0 ;q < fe_v.n_quadrature_points ; q++)
      {
       for (unsigned int i = 0 ; i < fe_v.dofs_per_cell ; i++)
       {
          for (unsigned int j = 0 ; j < fe_v.dofs_per_cell  ; j++)
          {
            // loop over convection
            for (unsigned int k = 0 ;k < dim ; k ++)
                if (this->exists(this->A[k].Row_Col_Value,component[i],component[j]))
                      cell_matrix(i,j) += fe_v.shape_value(i,q) * this->A[k].matrix.coeffRef(component[i],component[j]) 
                                          * fe_v.shape_grad(j,q)[k] * Jacobians_interior[q];
                
                  
            // loop over production
            if (this->exists(this->P.Row_Col_Value,component[i],component[j]))
                cell_matrix(i,j) += fe_v.shape_value(i,q) * this->P.matrix.coeffRef(component[i],component[j]) 
                                    * fe_v.shape_value(j,q) * Jacobians_interior[q];

            
          }     
       
       }
      }      

}

template<int dim> void Solver_DG<dim>::integrate_boundary_term(DoFInfo &dinfo,
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
                  Vector<double> boundary_rhs_value(this->nEqn);
                  this->build_BCrhs(fe_v.quadrature_point(q),fe_v.normal_vector(q),boundary_rhs_value);

                  // build the matrices needed
                  Eigen::MatrixXd Am = this->build_Aminus(fe_v.normal_vector(q));
                  Sparse_matrix Projector = this->build_Projector(fe_v.normal_vector(q));
                  Sparse_matrix Inv_Projector = this->build_InvProjector(fe_v.normal_vector(q));

                  Eigen::MatrixXd Am_invP_BC_P = Am * Inv_Projector 
                                                * this->BC.matrix * Projector;

                  Eigen::MatrixXd Am_invP = Am * Inv_Projector;

                  for (unsigned int i = 0 ; i < dofs_per_cell ; i ++)
                  {
                    for (unsigned int j = 0 ; j < dofs_per_cell ; j ++)
                        cell_matrix(i,j) += 0.5 * fe_v.shape_value(i,q) 
                                            * (Am(component[i],component[j])-Am_invP_BC_P(component[i],component[j]))
                                            * fe_v.shape_value(j,q) * Jacobian_face[q];                                    
                      

                    for (unsigned int j = 0 ; j < Am_invP.cols() ; j++)
                     cell_rhs(i) -= 0.5 * fe_v.shape_value(i,q) 
                                    * Am_invP(component[i],j) * boundary_rhs_value[j] * Jacobian_face[q];


                  }

                }
}

template<int dim> void Solver_DG<dim>::integrate_face_term(DoFInfo &dinfo1,
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
    Eigen::MatrixXd Am = this->build_Aminus(fe_v.normal_vector(q));
    Eigen::MatrixXd Am_neighbor = this->build_Aminus(-fe_v.normal_vector(q));


    for (unsigned int i = 0 ; i < dofs_per_cell ; i ++)
      for (unsigned int j = 0 ; j < dofs_per_cell ; j ++)
      {
        u1_v1(i,j) += 0.5 * fe_v.shape_value(i,q) * Am(component[i],component[j]) 
        * fe_v.shape_value(j,q) * Jacobian_face[q];

        u2_v1(i,j) -= 0.5 * fe_v.shape_value(i,q) * Am(component[i],component[j]) 
        * fe_v_neighbor.shape_value(j,q) * Jacobian_face[q];

        u2_v2(i,j) += 0.5 * fe_v_neighbor.shape_value(i,q) * Am_neighbor(component[i],component[j])
        * fe_v_neighbor.shape_value(j,q) * Jacobian_face[q];

        u1_v2(i,j) -= 0.5 * fe_v_neighbor.shape_value(i,q) * Am_neighbor(component[i],component[j]) 
        * fe_v.shape_value(j,q) * Jacobian_face[q];

      }

    }

}
