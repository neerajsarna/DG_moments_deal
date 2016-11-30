template<int dim>
void 
Base_Solver<dim>::assemble_rhs()
{

	const QGauss<dim> quadrature(ngp);
	const UpdateFlags update_flags  = update_values | update_JxW_values | update_quadrature_points;

	FEValues<dim>  fe_v(mapping,finite_element,quadrature, update_flags);
	typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();

	const unsigned int total_ngp = quadrature.size();
	const unsigned int dofs_per_cell = finite_element.dofs_per_cell;
    
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
	 Vector<double> cell_rhs(dofs_per_cell);
	 std::vector<double> Jacobians_interior(total_ngp);
	 std::vector<Vector<double>> source_term_value(total_ngp,Vector<double>(constants.nEqn));
    Vector<double> component(dofs_per_cell);

    for (unsigned int i = 0 ; i < dofs_per_cell ; i++)
        component[i] = finite_element.system_to_component_index(i).first;

	for (; cell != endc ; cell++) 
      {
      	cell_rhs = 0;
      	fe_v.reinit(cell);
        Jacobians_interior = fe_v.get_JxW_values();      	
        system_info->source_term(fe_v.get_quadrature_points(),source_term_value);

      	for (unsigned int q = 0 ; q < total_ngp ; q++)
      		for (unsigned int i = 0 ; i < dofs_per_cell ; i ++)
            cell_rhs(i) += fe_v.shape_value(i,q) * source_term_value[q][component[i]] * Jacobians_interior[q];

        cell->get_dof_indices(local_dof_indices);

        for(unsigned int i = 0 ; i < dofs_per_cell ; i++)
              system_rhs(local_dof_indices[i]) += cell_rhs(i);

      }    
}

template<int dim>
void 
Base_Solver<dim>
::assemble_system_meshworker()
{
	// first we will assemble the right hand side 
	Threads::Task<> rhs_task = Threads::new_task (&Base_Solver<dim>::assemble_rhs,
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

  // we will now initialize the shape values. We should now be carefull because of the share memory architecture
  this->Compute_Shape_Value(mapping,constants.p,ngp,cell);

  switch(constants.bc_type)
  {
    case odd:
    {
      MeshWorker::loop<dim, dim, MeshWorker::DoFInfo<dim>, MeshWorker::IntegrationInfoBox<dim> >
                          (cell, endc,
                           dof_info, info_box,
                           std_cxx11::bind(&Base_Solver<dim>::integrate_cell_term,
                            this,
                            std_cxx11::_1,std_cxx11::_2),
                           std_cxx11::bind(&Base_Solver<dim>::integrate_boundary_term_odd,
                            this,
                            std_cxx11::_1,
                            std_cxx11::_2),
                           std_cxx11::bind(&Base_Solver<dim>::integrate_face_term,
                            this,
                            std_cxx11::_1,
                            std_cxx11::_2,
                            std_cxx11::_3,
                            std_cxx11::_4),
                           assembler);
      break;
    }

    case characteristic:
    {
      MeshWorker::loop<dim, dim, MeshWorker::DoFInfo<dim>, MeshWorker::IntegrationInfoBox<dim> >
                          (cell, endc,
                           dof_info, info_box,
                           std_cxx11::bind(&Base_Solver<dim>::integrate_cell_term,
                            this,
                            std_cxx11::_1,std_cxx11::_2),
                           std_cxx11::bind(&Base_Solver<dim>::integrate_boundary_term_char,
                            this,
                            std_cxx11::_1,
                            std_cxx11::_2),
                           std_cxx11::bind(&Base_Solver<dim>::integrate_face_term,
                            this,
                            std_cxx11::_1,
                            std_cxx11::_2,
                            std_cxx11::_3,
                            std_cxx11::_4),
                           assembler);
      break;
    }
  }


     rhs_task.join();
}


template<int dim>
 void 
 Base_Solver<dim>::integrate_cell_term (DoFInfo &dinfo,
                                   			CellInfo &info)
{
	//if (count_meshworker > 0)
	//	Assert(1 == 0, ExcMessage("debugging integrate cell term")); 
 
	const FEValuesBase<dim> &fe_v = info.fe_values();
	FullMatrix<double> &cell_matrix = dinfo.matrix(0).matrix;
	const std::vector<double> &J = fe_v.get_JxW_values ();
	const FiniteElement<dim> &fe_in_cell = info.finite_element();

	const unsigned int dofs_per_cell = fe_in_cell.dofs_per_cell;
	const unsigned int components_per_cell  = fe_in_cell.n_components();
	const unsigned int indices_per_cell = dofs_per_cell/components_per_cell;

  Full_matrix Mass_x(indices_per_cell,indices_per_cell);
  Full_matrix Mass_y(indices_per_cell,indices_per_cell);
  Full_matrix Mass(indices_per_cell,indices_per_cell);
  Full_matrix cell_matrix_dummy(dofs_per_cell,dofs_per_cell);

  // matrix which stores Integrate phi * phi_x dx. Mass_y is the same but for the y-derivative and Mass is the normal mass
  // matrix.
  Mass_x = this->Compute_Mass_shape_grad_x(fe_v, indices_per_cell, J);
  Mass_y = this->Compute_Mass_shape_grad_y(fe_v, indices_per_cell, J);
  Mass = this->Compute_Mass_shape_value(fe_v, indices_per_cell, J);

  cell_matrix_dummy = matrix_opt.compute_A_outer_B(system_info->system_data.A[0].matrix,Mass_x) 
                    +  matrix_opt.compute_A_outer_B(system_info->system_data.A[1].matrix,Mass_y)
                    +  matrix_opt.compute_A_outer_B(system_info->system_data.P.matrix,Mass);


  // now we need to copy the memory, this should not be very slow and the overall algorithm should be faster than before
  for (unsigned int i = 0 ; i < dofs_per_cell ; i++)
    for (unsigned int j = 0 ; j < dofs_per_cell ; j++)
      cell_matrix(i,j) = cell_matrix_dummy(i,j);
}

template<int dim> 
void 
Base_Solver<dim>
::integrate_boundary_term_odd(DoFInfo &dinfo,
                              CellInfo &info)
{
 const FEValuesBase<dim> &fe_v = info.fe_values();
 typename Triangulation<dim>::face_iterator face_itr= dinfo.face;
 FullMatrix<double> &cell_matrix = dinfo.matrix(0).matrix;
 Vector<double> &cell_rhs = dinfo.vector(0).block(0);
 const std::vector<double> &Jacobian_face = fe_v.get_JxW_values ();

 const FiniteElement<dim> &fe_in_cell = info.finite_element();

 const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
 std::vector<unsigned int> component(dofs_per_cell);

 for (unsigned int i = 0 ; i < dofs_per_cell ; i ++)
  component[i] = fe_in_cell.system_to_component_index(i).first;

 Vector<double> boundary_rhs_value;

  boundary_rhs_value.reinit(constants.nEqn);

for (unsigned int q = 0 ; q < fe_v.n_quadrature_points ; q++)
{
  const double jacobian_value = Jacobian_face[q];

  boundary_rhs_value = 0;                 


  system_info->build_BCrhs(fe_v.quadrature_point(q),fe_v.normal_vector(q),
                          boundary_rhs_value,face_itr->boundary_id());

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
::integrate_boundary_term_char(DoFInfo &dinfo,
                          CellInfo &info)
{
 const FEValuesBase<dim> &fe_v = info.fe_values();
 typename Triangulation<dim>::face_iterator face_itr = dinfo.face;
 FullMatrix<double> &cell_matrix = dinfo.matrix(0).matrix;
 Vector<double> &cell_rhs = dinfo.vector(0).block(0);
 const std::vector<double> &Jacobian_face = fe_v.get_JxW_values ();

 const FiniteElement<dim> &fe_in_cell = info.finite_element();

 const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
 std::vector<unsigned int> component(dofs_per_cell);

 for (unsigned int i = 0 ; i < dofs_per_cell ; i ++)
  component[i] = fe_in_cell.system_to_component_index(i).first;

 Vector<double> boundary_rhs_value;

  boundary_rhs_value.reinit(constants.nBC);

for (unsigned int q = 0 ; q < fe_v.n_quadrature_points ; q++)
{
  const double jacobian_value = Jacobian_face[q];

  boundary_rhs_value = 0;                 

  system_info->build_BCrhs(fe_v.quadrature_point(q),fe_v.normal_vector(q),boundary_rhs_value,face_itr->boundary_id());

        // build the matrices needed
  Eigen::MatrixXd Am = system_info->build_Aminus(fe_v.normal_vector(q));
  Sparse_matrix Projector = system_info->build_Projector(fe_v.normal_vector(q));
  Sparse_matrix Inv_Projector = system_info->build_InvProjector(fe_v.normal_vector(q));

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




template<int dim>
void 
Base_Solver<dim>
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

  const std::vector<double> &J = fe_v.get_JxW_values ();

  const FiniteElement<dim> &fe_in_cell1 = info1.finite_element();
  const FiniteElement<dim> &fe_in_cell2 = info2.finite_element();
    
    const unsigned int dofs_per_cell1 = fe_in_cell1.dofs_per_cell;
    const unsigned int n_components1 = fe_in_cell1.n_components();
    const unsigned int dofs_per_component1 = dofs_per_cell1/n_components1;

    const unsigned int dofs_per_cell2 = fe_in_cell2.dofs_per_cell;
    const unsigned int n_components2 = fe_in_cell2.n_components();
    const unsigned int dofs_per_component2 = dofs_per_cell2/n_components2;

    AssertDimension(dofs_per_cell1,dofs_per_cell2);
    Full_matrix u1_v1_dummy(dofs_per_cell1,dofs_per_cell1);
    Full_matrix u1_v2_dummy(dofs_per_cell2,dofs_per_cell1);
    Full_matrix u2_v2_dummy(dofs_per_cell2,dofs_per_cell2);
    Full_matrix u2_v1_dummy(dofs_per_cell1,dofs_per_cell2);

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

    u1_v1_dummy = 0.5 * matrix_opt.compute_A_outer_B(Am,Mass_u1_v1);
    u2_v1_dummy = -0.5 * matrix_opt.compute_A_outer_B(Am,Mass_u2_v1);
    u1_v2_dummy = -0.5 * matrix_opt.compute_A_outer_B(Am_neighbor,Mass_u1_v2);
    u2_v2_dummy = 0.5 * matrix_opt.compute_A_outer_B(Am_neighbor,Mass_u2_v2);

    // now we need to copy into the dealii matrix
    for (unsigned int i = 0 ; i < dofs_per_cell1 ; i++)
      for (unsigned int j = 0 ; j < dofs_per_cell1 ; j++)
      {
        u1_v1(i,j)  = u1_v1_dummy(i,j);
        u2_v1(i,j)  = u2_v1_dummy(i,j);
        u1_v2(i,j)  = u1_v2_dummy(i,j);
        u2_v2(i,j)  = u2_v2_dummy(i,j);
      }
}
