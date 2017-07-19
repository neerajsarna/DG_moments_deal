namespace FEM_Solver
{
	using namespace dealii;

		template<int dim>
	class
	Assembly_Manager_FE:public NumericalIntegration::Base_NumericalIntegration<dim>
	{
	public:
		typedef MeshWorker::DoFInfo<dim> DoFInfo;
		typedef MeshWorker::IntegrationInfo<dim> CellInfo;



        DeclException2 (ExcCellCenter, double, double,
                        << "Cell Center = " << arg1 << "Neighbor Center = " << arg2 << "Mesh not lexiographical ");

        DeclException2 (ExcNoElementInSparsity, size_t, size_t,
                        << "Dof-1 " << arg1 << " Dof-2 " << arg2 << " Entry does not exist in sparsity pattern");

		Assembly_Manager_FE(const std::string &output_file_name,
						            const constant_numerics &constants,
						            std::vector<Develop_System::System<dim>> &equation_info,
                        const std::vector<int> &nEqn,
                        const std::vector<int> &nBC,
                        const unsigned int system_to_solve,
                        Triangulation<dim> &triangulation);
		     // finite element data structure 
     fe_data<dim> fe_data_structure;

			const int nEqn;
			const int nBC;
			const constant_numerics constants;

			Develop_System::System<dim> system_info;

			// A matrix in the final Ax = b
			TrilinosWrappers::SparseMatrix global_matrix;
            

			// The rhs in final Ax = b.
        	Vector<double> system_rhs;



        	const unsigned int ngp;
        	const unsigned int ngp_face;
            unsigned int integrate_inflow = 0;


		void assemble_system_meshworker();



      //Routines for meshworker assembly. 
      void integrate_cell_term (DoFInfo &dinfo,
        		CellInfo &info);

      void integrate_boundary_term_odd (DoFInfo &dinfo,
        									CellInfo &info);

     // integrate boundary term using characteristic variables
      void integrate_boundary_term_char(DoFInfo &dinfo,
        								CellInfo &info);

     void integrate_face_term(DoFInfo &dinfo1,
     			        		DoFInfo &dinfo2,
        						CellInfo &info1,
        						CellInfo &info2);

     void integrate_boundary_kinetic_flux(DoFInfo &dinfo,
                                    CellInfo &info);

     // integrate the face term using the kinetic flux
     void integrate_face_term_kinetic_flux(DoFInfo &dinfo1,
                                        DoFInfo &dinfo2,
                    CellInfo &info1,
                    CellInfo &info2);


     void assemble_rhs();

     MatrixOpt::Base_MatrixOpt matrix_opt;



	};

	template<int dim>
	Assembly_Manager_FE<dim>::Assembly_Manager_FE(const std::string &output_file_name,
											const constant_numerics &constants,
											std::vector<Develop_System::System<dim>> &equation_info,
                        					const std::vector<int> &nEqn,
                        					const std::vector<int> &nBC,          
                                  const unsigned int system_to_solve,
                                  Triangulation<dim> &triangulation)
	:
	fe_data_structure(output_file_name,
					         constants,nEqn[system_to_solve],triangulation),
  // we only save data corresponding to that system which is needed to be solved
	nEqn(nEqn[system_to_solve]),
	nBC(nBC[system_to_solve]),
	constants(constants),
	system_info(equation_info[system_to_solve]),
	ngp(constants.p + 1),
	ngp_face(constants.p + 1)
	{
  }


    // assembly the system using a meshworker
  template<int dim>
  void 
  Assembly_Manager_FE<dim>::assemble_system_meshworker()
  {
  // meshworker can only be used for a single system
  //AssertDimension(nEqn.size(),1);

  // first we will assemble the right hand side 
  Threads::Task<> rhs_task = Threads::new_task (&Assembly_Manager_FE<dim>::assemble_rhs,
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
  info_box.initialize(fe_data_structure.finite_element,fe_data_structure.mapping);
  MeshWorker::DoFInfo<dim> dof_info(fe_data_structure.dof_handler);

  MeshWorker::Assembler::SystemSimple<TrilinosWrappers::SparseMatrix, Vector<double>> assembler;
  assembler.initialize(global_matrix,system_rhs);


  typename DoFHandler<dim>::active_cell_iterator cell = fe_data_structure.dof_handler.begin_active();
  typename DoFHandler<dim>::active_cell_iterator endc = fe_data_structure.dof_handler.end();

  // we also initialize the values of the shape functions on the first cell. 
  // the values remain the same even for all the other cells.
  this->Compute_Shape_Value(fe_data_structure.mapping,ngp,cell);

  
  switch(this->constants.bc_type)
  {
    case odd:
    {
      MeshWorker::loop<dim, dim, MeshWorker::DoFInfo<dim>, MeshWorker::IntegrationInfoBox<dim> >
                          (cell, endc,
                           dof_info, info_box,
                           std_cxx11::bind(&Assembly_Manager_FE<dim>::integrate_cell_term,
                            this,
                            std_cxx11::_1,std_cxx11::_2),
                           std_cxx11::bind(&Assembly_Manager_FE<dim>::integrate_boundary_term_odd,
                            this,
                            std_cxx11::_1,
                            std_cxx11::_2),
                           std_cxx11::bind(&Assembly_Manager_FE<dim>::integrate_face_term,
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
                           std_cxx11::bind(&Assembly_Manager_FE<dim>::integrate_cell_term,
                            this,
                            std_cxx11::_1,std_cxx11::_2),
                           std_cxx11::bind(&Assembly_Manager_FE<dim>::integrate_boundary_term_char,
                            this,
                            std_cxx11::_1,
                            std_cxx11::_2),
                           std_cxx11::bind(&Assembly_Manager_FE<dim>::integrate_face_term,
                            this,
                            std_cxx11::_1,
                            std_cxx11::_2,
                            std_cxx11::_3,
                            std_cxx11::_4),
                           assembler);
      break;
    }

    case kinetic:
    {
      std::cout << "************ Using Kinetic Flux For Boundary *********" << std::endl;


      MeshWorker::loop<dim, dim, MeshWorker::DoFInfo<dim>, MeshWorker::IntegrationInfoBox<dim> >
        (cell, endc,
       dof_info, info_box,
       std_cxx11::bind(&Assembly_Manager_FE<dim>::integrate_cell_term,
        this,
        std_cxx11::_1,std_cxx11::_2),
       std_cxx11::bind(&Assembly_Manager_FE<dim>::integrate_boundary_kinetic_flux,
        this,
        std_cxx11::_1,
        std_cxx11::_2),
       std_cxx11::bind(&Assembly_Manager_FE<dim>::integrate_face_term,
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
  Assembly_Manager_FE<dim>::assemble_rhs()
  {


    const QGauss<dim> quadrature(ngp);
    const UpdateFlags update_flags  = update_values | update_JxW_values | update_quadrature_points;

    FEValues<dim>  fe_v(fe_data_structure.mapping,fe_data_structure.finite_element,quadrature, update_flags);
    typename DoFHandler<dim>::active_cell_iterator cell = fe_data_structure.dof_handler.begin_active(), 
                             endc = fe_data_structure.dof_handler.end();

    const unsigned int total_ngp = quadrature.size();
    const unsigned int dofs_per_cell = fe_data_structure.finite_element.dofs_per_cell;

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);
    std::vector<double> Jacobians_interior(total_ngp);
    std::vector<Vector<double>> source_term_value(total_ngp,Vector<double>(nEqn));
    Vector<double> component(dofs_per_cell);

    for (unsigned int i = 0 ; i < dofs_per_cell ; i++)
      component[i] = fe_data_structure.finite_element.system_to_component_index(i).first;

    for (; cell != endc ; cell++) 
    {
      cell_rhs = 0;
      fe_v.reinit(cell);
      Jacobians_interior = fe_v.get_JxW_values();       
      system_info.source_term(fe_v.get_quadrature_points(),source_term_value);

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
  Assembly_Manager_FE<dim>::integrate_cell_term (DoFInfo &dinfo,
                        CellInfo &info)
  {

    const FEValuesBase<dim> &fe_v = info.fe_values();
    FullMatrix<double> &cell_matrix = dinfo.matrix(0).matrix;
    const std::vector<double> &Jacobians_interior = fe_v.get_JxW_values ();
    const FiniteElement<dim> &fe_in_cell = info.finite_element();

    const unsigned int dofs_per_cell = fe_in_cell.dofs_per_cell;
    const unsigned int components_per_cell  = fe_in_cell.n_components();
    const unsigned int indices_per_cell = dofs_per_cell/components_per_cell;
    std::vector<std::vector<double>> component_to_system(components_per_cell,std::vector<double> (indices_per_cell));

    for (unsigned int i = 0 ; i < components_per_cell ; i ++)
      for (unsigned int j = 0 ; j < indices_per_cell ; j ++)
        component_to_system[i][j] = fe_in_cell.component_to_system_index(i,j); 

      std::vector<FullMatrix<double>> Mass_grad(dim);
      FullMatrix<double> Mass(indices_per_cell,indices_per_cell);

      Mass_grad = this->Compute_Mass_shape_grad(fe_v, indices_per_cell, Jacobians_interior);
      Mass = this->Compute_Mass_shape_value(fe_v, indices_per_cell, Jacobians_interior);

        // this is for initializing the matrix, it is necessary since we have not put cell_matrix to zero
      cell_matrix = this->matrix_opt.compute_A_outer_B(system_info.system_data.P.matrix,Mass);

      for (int space = 0 ;space < dim ; space ++ )
        cell_matrix.add(0,cell_matrix,1,this->matrix_opt.compute_A_outer_B(system_info.system_data.A[space].matrix,Mass_grad[space]));

  }

  // boundary implelmentaion using odd splitting 
  template<int dim> 
  void 
  Assembly_Manager_FE<dim>
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
    boundary_rhs_value.reinit(nBC);
    const unsigned int b_id = face_itr->boundary_id();

    Assert(system_info.system_data.B.matrix.rows() == system_info.system_data.Sigma.matrix.cols(),
            ExcMessage("Incorrect dimension"));

// we use a temporary matrix to determine whether inflow or outflow
    Sparse_matrix B_temp;

    if(b_id == 101 || b_id == 102)
    {
      integrate_inflow++;
      B_temp = system_info.system_data.Binflow.matrix;
    }
    else
    {
    // B matrix for specular reflection
      if (b_id == 50)
        B_temp = system_info.system_data.B_specular.matrix;

    // B matrix for full accmmodation
      else
        B_temp = system_info.system_data.B.matrix;
    }

    for (unsigned int q = 0 ; q < fe_v.n_quadrature_points ; q++)
    {
      const double jacobian_value = Jacobian_face[q];

      boundary_rhs_value = 0;                 

  // check for inflow or outflow
  // Incase of inflow provide the inflow rhs
      if(face_itr->boundary_id() == 101 || face_itr->boundary_id() == 102)
        system_info.bcrhs_inflow.BCrhs(fe_v.quadrature_point(q),fe_v.normal_vector(q),
          boundary_rhs_value,face_itr->boundary_id());

      else
        system_info.bcrhs_wall.BCrhs(fe_v.quadrature_point(q),fe_v.normal_vector(q),
          boundary_rhs_value,face_itr->boundary_id());


      Sparse_matrix Projector = system_info.build_Projector(fe_v.normal_vector(q));



  // Simga(as given in PDF) * B * Projector
  // Sigma in the PDF = Projector.transpose * Sigma(In the code)
      Full_matrix Sigma_B_P =       Projector.transpose()
      * system_info.system_data.Sigma.matrix 
      * B_temp
      * Projector;

      Full_matrix Sigma = Projector.transpose()
      * system_info.system_data.Sigma.matrix ;

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


  // boundary implementation using characteristic variables
  template<int dim> 
  void 
  Assembly_Manager_FE<dim>
  ::integrate_boundary_term_char(DoFInfo &dinfo,
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
    boundary_rhs_value.reinit(nBC);
    const unsigned int b_id = face_itr->boundary_id();

    Assert(system_info.system_data.B.matrix.rows() == system_info.system_data.Sigma.matrix.cols() ,
            ExcMessage("Incorrect dimension"));

// we use a temporary matrix to determine whether inflow or outflow
    Sparse_matrix B_temp;
    Full_matrix penalty_temp;

    if(b_id == 101 || b_id == 102)
    {
      integrate_inflow++;
      B_temp = system_info.system_data.Binflow.matrix;
      penalty_temp = system_info.penalty_char_inflow;
    }
    else
    {
    // B matrix for specular reflection
      if (b_id == 50)
        B_temp = system_info.system_data.B_specular.matrix;

    // B matrix for full accmmodation
      else
      {
        B_temp = system_info.system_data.B.matrix;
        penalty_temp = system_info.penalty_char_wall;
      }
      
    }

    for (unsigned int q = 0 ; q < fe_v.n_quadrature_points ; q++)
    {
      const double jacobian_value = Jacobian_face[q];
      Tensor<1,dim> outward_normal = fe_v.normal_vector(q);

      boundary_rhs_value = 0;                 


  // check for inflow or outflow
  // Incase of inflow provide the inflow rhs
      if(face_itr->boundary_id() == 101 || face_itr->boundary_id() == 102)
        system_info.bcrhs_inflow.BCrhs(fe_v.quadrature_point(q),outward_normal,
          boundary_rhs_value,face_itr->boundary_id());

      else
        system_info.base_bcrhs->BCrhs(fe_v.quadrature_point(q),outward_normal,
                                      boundary_rhs_value,face_itr->boundary_id());


      Sparse_matrix Projector = system_info.build_Projector(outward_normal);
      Sparse_matrix Inv_Projector = system_info.build_InvProjector(outward_normal);
      Eigen::MatrixXd Am = system_info.build_Aminus(outward_normal);


      // the projector was not included in the development of the penalty matrix
      Full_matrix Am_Penalty_B = 0.5 * Am * Inv_Projector * 
                                 penalty_temp * B_temp * Projector;

      Full_matrix Am_Penalty = 0.5 * Am * Inv_Projector * 
                                          penalty_temp;

      // boundary condition has been implemented in the form BU-g. Therefore the right hand side 
      // gets a plus sign
      for (unsigned int i = 0 ; i < dofs_per_cell ; i ++)
      {
        const double shape_value_test = fe_v.shape_value(i,q);
        for (unsigned int j = 0 ; j < dofs_per_cell ; j ++)
          cell_matrix(i,j) += shape_value_test
                              * Am_Penalty_B(component[i],component[j])
                              * fe_v.shape_value(j,q) 
                              * jacobian_value;                                    


        for (unsigned int j = 0 ; j < boundary_rhs_value.size() ; j++)
          cell_rhs(i) += shape_value_test 
                          * Am_Penalty(component[i],j) 
                          * boundary_rhs_value[j] 
                          * jacobian_value;

      }


    }

  }


  template<int dim> 
  void 
  Assembly_Manager_FE<dim>
  ::integrate_boundary_kinetic_flux(DoFInfo &dinfo,
                                    CellInfo &info)
  {

    AssertDimension(dim,1);

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
    boundary_rhs_value.reinit(nEqn);
    const unsigned int b_id = face_itr->boundary_id();

    Sparse_matrix rho_temp;

    if(b_id == 101 || b_id == 102)
    {
      integrate_inflow++;
      rho_temp = system_info.system_data.rhoInflow.matrix;
    }
    else
    {


      Assert(b_id != 50,ExcNotImplemented());

    // B matrix for full accmmodation
      if (b_id != 50)
        rho_temp = system_info.system_data.rhoW.matrix;
    }

    for (unsigned int q = 0 ; q < fe_v.n_quadrature_points ; q++)
    {
      const double jacobian_value = Jacobian_face[q];
      Tensor<1,dim> outward_normal = fe_v.normal_vector(q);
      Sparse_matrix Projector = system_info.build_Projector(outward_normal);

      boundary_rhs_value = 0;                 

  // check for inflow or outflow
  // Incase of inflow provide the inflow rhs
      if(face_itr->boundary_id() == 101 || face_itr->boundary_id() == 102)
        system_info.bcrhs_inflow_kinetic.BCrhs(fe_v.quadrature_point(q),fe_v.normal_vector(q),
                                                boundary_rhs_value,face_itr->boundary_id());

      else
        system_info.bcrhs_wall_kinetic.BCrhs(fe_v.quadrature_point(q),fe_v.normal_vector(q),
                                                boundary_rhs_value,face_itr->boundary_id());

  // Simga(as given in PDF) * B * Projector
  // Sigma in the PDF = Projector.transpose * Sigma(In the code)
      Eigen::MatrixXd Am = system_info.build_Aminus(fe_v.normal_vector(q));

      Eigen::MatrixXd Am_rho = Am * rho_temp * Projector;

      for (unsigned int i = 0 ; i < dofs_per_cell ; i ++)
      {
        const double shape_value_test = fe_v.shape_value(i,q);
        for (unsigned int j = 0 ; j < dofs_per_cell ; j ++)
          cell_matrix(i,j) += 0.5 * shape_value_test
                              * (Am(component[i],component[j]) - Am_rho(component[i],component[j]))
                              * fe_v.shape_value(j,q) 
                              * jacobian_value;                                    


        for (unsigned int j = 0 ; j < boundary_rhs_value.size() ; j++)
          cell_rhs(i) += 0.5 * shape_value_test 
                        * Am(component[i],j) 
                        * boundary_rhs_value[j] 
                        * jacobian_value;

      }


    }

  }



  template<int dim>
  void 
  Assembly_Manager_FE<dim>
  ::integrate_face_term(DoFInfo &dinfo1,
            DoFInfo &dinfo2,
            CellInfo &info1,
            CellInfo &info2)
  {

    const FEValuesBase<dim> &fe_v = info1.fe_values();
    const FEValuesBase<dim> &fe_v_neighbor = info2.fe_values();

  // jacobian values for the face of the current cell
    const std::vector<double> &J = fe_v.get_JxW_values ();

    const FiniteElement<dim> &fe_in_cell1 = fe_v.get_fe();
    const FiniteElement<dim> &fe_in_cell2 = fe_v_neighbor.get_fe();

    FullMatrix<double> &u1_v1 = dinfo1.matrix(0,false).matrix;
    FullMatrix<double> &u2_v1 = dinfo1.matrix(0,true).matrix;
    FullMatrix<double> &u1_v2 = dinfo2.matrix(0,true).matrix;
    FullMatrix<double> &u2_v2 = dinfo2.matrix(0,false).matrix;

    const unsigned int dofs_per_cell1 = fe_in_cell1.dofs_per_cell;
    const unsigned int n_components1 = fe_in_cell1.n_components();
    const unsigned int dofs_per_component1 = dofs_per_cell1/n_components1;

    const unsigned int dofs_per_cell2 = fe_in_cell2.dofs_per_cell;
    const unsigned int n_components2 = fe_in_cell2.n_components();
    const unsigned int dofs_per_component2 = dofs_per_cell2/n_components2;


    //CAUTION: ASSUMPTION OF STRAIGHT EDGES IN THE INTERIOR
    Tensor<1,dim> outward_normal = fe_v.normal_vector(0);
    Eigen::MatrixXd Am = system_info.build_Aminus(outward_normal);
    Eigen::MatrixXd Am_neighbor = system_info.build_Aminus(-outward_normal);


    FullMatrix<double> Mass_u1_v1(dofs_per_component1,dofs_per_component1);
    FullMatrix<double> Mass_u2_v1(dofs_per_component1,dofs_per_component2);
    FullMatrix<double> Mass_u2_v2(dofs_per_component2,dofs_per_component2);
    FullMatrix<double> Mass_u1_v2(dofs_per_component2,dofs_per_component1);

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

    u1_v1 = matrix_opt.multiply_scalar(0.5,matrix_opt.compute_A_outer_B(Am,Mass_u1_v1));
    u2_v1 = matrix_opt.multiply_scalar(-0.5,matrix_opt.compute_A_outer_B(Am,Mass_u2_v1));
    u1_v2 = matrix_opt.multiply_scalar(-0.5,matrix_opt.compute_A_outer_B(Am_neighbor,Mass_u1_v2));
    u2_v2 = matrix_opt.multiply_scalar(0.5,matrix_opt.compute_A_outer_B(Am_neighbor,Mass_u2_v2));

  }

  template<int dim>
  void 
  Assembly_Manager_FE<dim>
  ::integrate_face_term_kinetic_flux(DoFInfo &dinfo1,
            DoFInfo &dinfo2,
            CellInfo &info1,
            CellInfo &info2)
  {

    const FEValuesBase<dim> &fe_v = info1.fe_values();
    const FEValuesBase<dim> &fe_v_neighbor = info2.fe_values();

  // jacobian values for the face of the current cell
    const std::vector<double> &J = fe_v.get_JxW_values ();

    const FiniteElement<dim> &fe_in_cell1 = fe_v.get_fe();
    const FiniteElement<dim> &fe_in_cell2 = fe_v_neighbor.get_fe();

    FullMatrix<double> &u1_v1 = dinfo1.matrix(0,false).matrix;
    FullMatrix<double> &u2_v1 = dinfo1.matrix(0,true).matrix;
    FullMatrix<double> &u1_v2 = dinfo2.matrix(0,true).matrix;
    FullMatrix<double> &u2_v2 = dinfo2.matrix(0,false).matrix;

    const unsigned int dofs_per_cell1 = fe_in_cell1.dofs_per_cell;
    const unsigned int n_components1 = fe_in_cell1.n_components();
    const unsigned int dofs_per_component1 = dofs_per_cell1/n_components1;

    const unsigned int dofs_per_cell2 = fe_in_cell2.dofs_per_cell;
    const unsigned int n_components2 = fe_in_cell2.n_components();
    const unsigned int dofs_per_component2 = dofs_per_cell2/n_components2;


    //CAUTION: ASSUMPTION OF STRAIGHT EDGES IN THE INTERIOR
    Tensor<1,dim> outward_normal = fe_v.normal_vector(0);
    Eigen::MatrixXd Am = system_info.build_Aminus_kinetic(outward_normal);
    Eigen::MatrixXd Am_neighbor = system_info.build_Aminus_kinetic(-outward_normal);


    FullMatrix<double> Mass_u1_v1(dofs_per_component1,dofs_per_component1);
    FullMatrix<double> Mass_u2_v1(dofs_per_component1,dofs_per_component2);
    FullMatrix<double> Mass_u2_v2(dofs_per_component2,dofs_per_component2);
    FullMatrix<double> Mass_u1_v2(dofs_per_component2,dofs_per_component1);

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

    u1_v1 = matrix_opt.multiply_scalar(0.5,matrix_opt.compute_A_outer_B(Am,Mass_u1_v1));
    u2_v1 = matrix_opt.multiply_scalar(-0.5,matrix_opt.compute_A_outer_B(Am,Mass_u2_v1));
    u1_v2 = matrix_opt.multiply_scalar(-0.5,matrix_opt.compute_A_outer_B(Am_neighbor,Mass_u1_v2));
    u2_v2 = matrix_opt.multiply_scalar(0.5,matrix_opt.compute_A_outer_B(Am_neighbor,Mass_u2_v2));

  }


	template<int dim>
	class
	Assembly_Manager_hp_FE:public NumericalIntegration::Base_NumericalIntegration<dim>
	{
	public:
		typedef MeshWorker::DoFInfo<dim> DoFInfo;
		typedef MeshWorker::IntegrationInfo<dim> CellInfo;



        DeclException2 (ExcCellCenter, double, double,
                        << "Cell Center = " << arg1 << "Neighbor Center = " << arg2 << "Mesh not lexiographical ");

        DeclException2 (ExcNoElementInSparsity, size_t, size_t,
                        << "Dof-1 " << arg1 << " Dof-2 " << arg2 << " Entry does not exist in sparsity pattern");

		Assembly_Manager_hp_FE(const std::string &output_file_name,
						const constant_numerics &constants,
						std::vector<Develop_System::System<dim>> &equation_info,
                        const std::vector<int> &nEqn,
                        const std::vector<int> &nBC,
                        Triangulation<dim> &triangulation);

     		hp_fe_data<dim> hp_fe_data_structure;
			const std::vector<int> nEqn;
			const std::vector<int> nBC;
			const constant_numerics constants;

			std::vector<Develop_System::System<dim>> system_info;

			// A matrix in the final Ax = b
			TrilinosWrappers::SparseMatrix global_matrix;
            

			// The rhs in final Ax = b.
        	Vector<double> system_rhs;



        	const unsigned int ngp;
        	const unsigned int ngp_face;
            unsigned int integrate_inflow = 0;



		void assemble_system_manuel();


		// routines for manuel assembly
       	void integrate_cell_manuel(FullMatrix<double> &cell_matrix,
                                        Vector<double> &cell_rhs,
                                        const FEValuesBase<dim> &fe_v,
                                        std::vector<double> &J,
                                        std::vector<Vector<double>> &source_term_value,
                                        std::vector<Vector<double>> &component_to_system,
                                        const unsigned int fe_index);

       // // integrate the boundary manuelly using odd boundary implementation
       void integrate_boundary_manuel_odd(FullMatrix<double> &cell_matrix,
                                          Vector<double> &cell_rhs,
                                          const FEValuesBase<dim> &fe_v,
                                          std::vector<double> &J,
                                          std::vector<Vector<double>> &component_to_system,
                                          const unsigned int b_id,
                                          const unsigned int fe_index);
        	
      void integrate_face_manuel(FullMatrix<double> &u1_v1,
                        FullMatrix<double> &u1_v2,
                        FullMatrix<double> &u2_v1,
                        FullMatrix<double> &u2_v2,
                        const FEValuesBase<dim> &fe_v,
                        const FEValuesBase<dim> &fe_v_neighbor,
                        std::vector<double> &J,
                        const unsigned int fe_index1,
                        const unsigned int fe_index2);



     // integrate boundary term using characteristic variables
      void integrate_boundary_term_char(DoFInfo &dinfo,
        								CellInfo &info);

 

     MatrixOpt::Base_MatrixOpt matrix_opt;


	};

	template<int dim>
	Assembly_Manager_hp_FE<dim>::Assembly_Manager_hp_FE(const std::string &output_file_name,
											const constant_numerics &constants,
											std::vector<Develop_System::System<dim>> &equation_info,
                        					const std::vector<int> &nEqn,
                        					const std::vector<int> &nBC,
                                  Triangulation<dim> &triangulation)
	:
	hp_fe_data_structure(output_file_name,
					constants,nEqn,triangulation),
	nEqn(nEqn),
	nBC(nBC),
	constants(constants),
	system_info(equation_info),
	ngp(constants.p + 1),
	ngp_face(constants.p + 1)
	{}


	template<int dim>
	void 
	Assembly_Manager_hp_FE<dim>::assemble_system_manuel()
	{
	// we check whether the moment systems are using the same polynomial degree or not
    for (unsigned long int i = 0 ; i < nEqn.size()-1; i++)
      AssertDimension(hp_fe_data_structure.finite_element[i].dofs_per_cell/nEqn[i],
      				 hp_fe_data_structure.finite_element[i+1].dofs_per_cell/nEqn[i+1]);

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

      hp::FEValues<dim>  hp_fe_v(hp_fe_data_structure.mapping,
      							 hp_fe_data_structure.finite_element,quadrature, update_flags);

      hp::FEFaceValues<dim> hp_fe_v_face(hp_fe_data_structure.mapping,
      									 hp_fe_data_structure.finite_element, face_quadrature, face_update_flags);

      hp::FEFaceValues<dim> hp_fe_v_face_neighbor(hp_fe_data_structure.mapping,
      											hp_fe_data_structure.finite_element,
      											 face_quadrature, neighbor_face_update_flags);

      hp::FESubfaceValues<dim> hp_fe_v_subface_neighbor(hp_fe_data_structure.mapping,
      													hp_fe_data_structure.finite_element, 
      													face_quadrature, neighbor_face_update_flags);  


      // the total number of quadrature points are the same for both the quadrature routines
      const unsigned int total_ngp = quadrature_basic.size();
      const unsigned int total_ngp_face = face_quadrature_basic.size();

      // iterator over the cells 
      typename hp::DoFHandler<dim>::active_cell_iterator cell = hp_fe_data_structure.dof_handler.begin_active(), 
      													 endc = hp_fe_data_structure.dof_handler.end();
      typename hp::DoFHandler<dim>::cell_iterator neighbor;


      // it is much more efficient to declare the memory and then use it again and again
      std::vector<FullMatrix<double>> cell_matrix;                  // contribution from the interior
      std::vector<FullMatrix<double>> boundary_matrix;
      std::vector<Vector<double>>   cell_rhs;                                   // rhs from the current cell


      for (unsigned long int i = 0 ; i < nEqn.size() ; i++)
      {
        cell_matrix.push_back(FullMatrix<double>(hp_fe_data_structure.dofs_per_cell[i],hp_fe_data_structure.dofs_per_cell[i]));
        boundary_matrix.push_back(FullMatrix<double>(hp_fe_data_structure.dofs_per_cell[i],hp_fe_data_structure.dofs_per_cell[i]));
        cell_rhs.push_back(Vector<double>(hp_fe_data_structure.dofs_per_cell[i]));
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
          component_to_system[i].push_back(Vector<double>(hp_fe_data_structure.dofs_per_component));

      // now we allocate the values for component_to_system
      for (unsigned long int i = 0 ; i < nEqn.size(); i++)
          for (int k = 0 ; k < nEqn[i] ;k ++)
              for (unsigned int j = 0 ; j < hp_fe_data_structure.dofs_per_component ; j++)
                component_to_system[i][k](j) = hp_fe_data_structure.finite_element[i].component_to_system_index(k,j);
           

      // source term
      std::vector<std::vector<Vector<double>>> source_term_value(nEqn.size());

      for (unsigned long int i = 0 ; i < nEqn.size() ; i ++)
        for (unsigned int j = 0 ; j < total_ngp ; j++)
          source_term_value[i].push_back(Vector<double>(nEqn[i]));
      
     
      // need to initialize fe_v so as to save the values of the shape functions
      hp_fe_v.reinit(cell);
      this->Compute_Shape_Value(hp_fe_v,hp_fe_data_structure.dofs_per_component);
      

  for (; cell != endc ; cell++) 
      {
        const unsigned int fe_index = cell->active_fe_index();
        const unsigned int dofs_this_cell = hp_fe_data_structure.dofs_per_cell[fe_index];

        Assert(fe_index <= hp_fe_data_structure.max_fe_index ,ExcNotImplemented());
        cell_rhs[fe_index] = 0;


        // the mapping from the dofs of the present cell to the global dofs
        std::vector<types::global_dof_index> local_dof_indices(dofs_this_cell);
        cell->get_dof_indices(local_dof_indices);

        
        hp_fe_v.reinit(cell);
        const FEValues<dim> &fe_v = hp_fe_v.get_present_fe_values();        

        Jacobians_interior = fe_v.get_JxW_values();
        this->system_info[fe_index].source_term(fe_v.get_quadrature_points(),source_term_value[fe_index]);

        
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
                for (unsigned int i = 0 ; i < dofs_this_cell ; i++)
                  for (unsigned int j = 0 ; j < dofs_this_cell ; j++)
                   global_matrix.add(local_dof_indices[i],local_dof_indices[j],boundary_matrix[fe_index](i,j));
                
              }
             
              
                // if the face is not at the boundary then we need to integrate over the face
                else
               {            
                
                 neighbor = cell->neighbor(face);


                 const unsigned int fe_index_neighbor = neighbor->active_fe_index();
                 const unsigned int dofs_neighbor = hp_fe_data_structure.dofs_per_cell[fe_index_neighbor];
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
              

            
              for (unsigned int i = 0 ; i < dofs_this_cell ; i++)
                  for (unsigned int j = 0 ; j < dofs_this_cell ; j++)          
                    global_matrix.add(local_dof_indices[i],local_dof_indices[j],cell_matrix[fe_index](i,j));

              
             for(unsigned int i = 0 ; i < dofs_this_cell ; i++)
                  system_rhs(local_dof_indices[i]) += cell_rhs[fe_index](i);

  
      }
	}


	template<int dim>
	void 
	Assembly_Manager_hp_FE<dim>
	::integrate_cell_manuel(FullMatrix<double> &cell_matrix,
		Vector<double> &cell_rhs,
		const FEValuesBase<dim> &fe_v,
		std::vector<double> &J,
		std::vector<Vector<double>> &source_term_value,
		std::vector<Vector<double>> &component_to_system,
		const unsigned int fe_index)
	{

		Assert(cell_matrix.m() !=0 || cell_matrix.n() !=0 ,ExcNotInitialized());
		Assert(cell_matrix.m() == cell_matrix.n() ,ExcMessage("different dof in same cell"));
		Assert(cell_rhs.size() != 0,ExcNotInitialized());
		Assert(J.size() !=0,ExcNotInitialized());
		Assert(source_term_value.size() != 0,ExcNotInitialized());
		Assert(fe_index < nEqn.size(),ExcMessage("Incorrect fe index"));
		AssertDimension((int)component_to_system.size(),nEqn[fe_index]);
		AssertDimension(cell_matrix.m(),hp_fe_data_structure.dofs_per_cell[fe_index]);
		AssertDimension(source_term_value.size(),J.size());
		AssertDimension(fe_v.get_fe().dofs_per_cell,hp_fe_data_structure.dofs_per_cell[fe_index]);


		const unsigned int total_ngp = J.size();


		AssertDimension(hp_fe_data_structure.dofs_per_component,
						this->shape_values.rows());
		AssertDimension(total_ngp,this->shape_values.cols());

		std::vector<FullMatrix<double>> Mass_grad(dim);

		// mass matrix corresponding to one particular set of basis functions
		FullMatrix<double> Mass(hp_fe_data_structure.dofs_per_component,
								hp_fe_data_structure.dofs_per_component);

		const int components_per_cell =nEqn[fe_index];
		Mass_grad = this->Compute_Mass_shape_grad(fe_v, hp_fe_data_structure.dofs_per_component, J);
		Mass = this->Compute_Mass_shape_value(fe_v, hp_fe_data_structure.dofs_per_component, J);

  // this is for initializing the matrix, it is necessary since we have not put cell_matrix to zero
		cell_matrix = matrix_opt.compute_A_outer_B(this->system_info[fe_index].system_data.P.matrix,Mass);

		for (int space = 0 ;space < dim ; space ++ )
			cell_matrix.add(0,cell_matrix,1,matrix_opt.compute_A_outer_B(this->system_info[fe_index].system_data.A[space].matrix,Mass_grad[space]));


		AssertDimension(hp_fe_data_structure.dofs_per_component,this->shape_values.rows());
		AssertDimension(total_ngp,this->shape_values.cols());

  // loop over all the gauss points
		for (unsigned int q = 0 ; q < total_ngp ; q++)
		{
    // value of the jacobian at the gauss point
			double jacobian_value = J[q];

    // index is the local id of the degree of freedom
			for (unsigned int index_test = 0 ; index_test < this->hp_fe_data_structure.dofs_per_component ; index_test++)
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

	}



//Integrate the boundary term using manuel implementation
	template<int dim>
	void 
	Assembly_Manager_hp_FE<dim>
	::integrate_boundary_manuel_odd(FullMatrix<double> &cell_matrix,
								Vector<double> &cell_rhs,
								const FEValuesBase<dim> &fe_v,
								std::vector<double> &J,
								std::vector<Vector<double>> &component_to_system,
								const unsigned int b_id,
								const unsigned int fe_index)
	{


		Assert(cell_matrix.m() !=0 || cell_matrix.n() !=0 ,ExcNotInitialized());
		Assert(cell_matrix.m() == cell_matrix.n() ,ExcMessage("different dof in same cell"));
		Assert(cell_rhs.size() != 0,ExcNotInitialized());
		Assert(J.size() !=0,ExcNotInitialized());


		Vector<double> boundary_rhs_value;

		boundary_rhs_value.reinit(nBC[fe_index]);

    // we use a temporary matrix to determine whether inflow or outflow
		Sparse_matrix B_temp;

		if(b_id == 101 || b_id == 102)
		{
			integrate_inflow++;
			B_temp = this->system_info[fe_index].system_data.Binflow.matrix;
		}
		else
		{
    // B matrix for specular reflection
			if (b_id == 50)
				B_temp = this->system_info[fe_index].system_data.B_specular.matrix;

    // B matrix for full accmmodation
			else
				B_temp = this->system_info[fe_index].system_data.B.matrix;
		}

		for (unsigned int q = 0 ; q < fe_v.n_quadrature_points ; q++)
		{
			const double jacobian_value = J[q];

			boundary_rhs_value = 0;                 


  // check for inflow or outflow
  // Incase of inflow provide the inflow rhs
			if(b_id == 101 || b_id == 102)
				this->system_info[fe_index].bcrhs_inflow.BCrhs(fe_v.quadrature_point(q),fe_v.normal_vector(q),
					boundary_rhs_value,b_id);

			else
				this->system_info[fe_index].bcrhs_wall.BCrhs(fe_v.quadrature_point(q),fe_v.normal_vector(q),
					boundary_rhs_value,b_id);         


			Sparse_matrix Projector = this->system_info[fe_index].build_Projector(fe_v.normal_vector(q));



  // Simga(as given in PDF) * B * Projector
  // Sigma in the PDF = Projector.transpose * Sigma(In the code)
			Full_matrix Sigma_B_P =       Projector.transpose()
			* this->system_info[fe_index].system_data.Sigma.matrix 
			* B_temp
			* Projector;

			Full_matrix Sigma = Projector.transpose()
			* this->system_info[fe_index].system_data.Sigma.matrix ;


  // we first loop over all the components of the test function
			for (int component_test = 0; component_test < nEqn[fe_index] ; component_test++)
    // we now loop over all the dof for our component
				for (unsigned long int index_test = 0 ; index_test < this->hp_fe_data_structure.dofs_per_component ; index_test ++)
				{
					const unsigned long int dof_test = component_to_system[component_test][index_test];
					const double shape_value_test = fe_v.shape_value(dof_test,q);

      // we now loop over all the components of the solution 
					for (int component_sol = 0 ; component_sol < nEqn[fe_index] ; component_sol++)

        // a loop over all the indices of the solution vector
						for (unsigned long int index_sol = 0 ; index_sol < this->hp_fe_data_structure.dofs_per_component ; index_sol++)
						{
							const unsigned long int dof_sol = component_to_system[component_sol][index_sol];
							const double shape_value_sol = fe_v.shape_value(dof_sol,q);
							cell_matrix(dof_test,dof_sol) -= shape_value_test
							* Sigma_B_P(component_test,component_sol)
							* shape_value_sol 
							* jacobian_value;  

						}

						for (unsigned int j = 0 ; j < boundary_rhs_value.size() ; j++)
							cell_rhs(dof_test) -= shape_value_test 
						* Sigma(component_test,j) 
						* boundary_rhs_value[j] 
						* jacobian_value;


					}



				}


			}

			// integrate the face term using manuel computation
			template<int dim>
			void 
			Assembly_Manager_hp_FE<dim>
			::integrate_face_manuel(FullMatrix<double> &u1_v1,
									FullMatrix<double> &u1_v2,
									FullMatrix<double> &u2_v1,
									FullMatrix<double> &u2_v2,
									const FEValuesBase<dim> &fe_v,
									const FEValuesBase<dim> &fe_v_neighbor,
									std::vector<double> &J,
									const unsigned int fe_index1,
									const unsigned int fe_index2)
			{

    // all of them must be square matrices
				Assert(u1_v1.m()!=0 && u1_v1.n() != 0,ExcNotInitialized());
				Assert(u1_v2.m() != 0 && u1_v2.n() != 0,ExcNotInitialized());
				Assert(u2_v1.m()!=0 && u2_v1.n()!=0,ExcNotInitialized());
				Assert(u2_v2.m()!=0 && u2_v2.n()!=0,ExcNotInitialized());
				Assert(J.size() !=0,ExcNotInitialized());


    //CAUTION: ASSUMPTION OF STRAIGHT EDGES IN THE INTERIOR
				Tensor<1,dim> outward_normal = fe_v.normal_vector(0);
				const unsigned int max_fe_index = std::max(fe_index1,fe_index2);
				Eigen::MatrixXd Am = this->system_info[max_fe_index].build_Aminus(outward_normal);
				Eigen::MatrixXd Am_neighbor = this->system_info[max_fe_index].build_Aminus(-outward_normal);
				const unsigned int nEqn_this = nEqn[fe_index1];
				const unsigned int nEqn_neighbor = nEqn[fe_index2];

    // mass matrices. Assumption is that all the components of the solution vector are being 
    // solved using the same polynomial degree

				FullMatrix<double> Mass_u1_v1(hp_fe_data_structure.dofs_per_component,hp_fe_data_structure.dofs_per_component);
				FullMatrix<double> Mass_u2_v1(hp_fe_data_structure.dofs_per_component,hp_fe_data_structure.dofs_per_component);
				FullMatrix<double> Mass_u2_v2(hp_fe_data_structure.dofs_per_component,hp_fe_data_structure.dofs_per_component);
				FullMatrix<double> Mass_u1_v2(hp_fe_data_structure.dofs_per_component,hp_fe_data_structure.dofs_per_component);

				Mass_u1_v1 = this->Compute_Mass_cell_neighbor(fe_v,
															fe_v,hp_fe_data_structure.dofs_per_component,
															hp_fe_data_structure.dofs_per_component,J);

				Mass_u2_v1 = this->Compute_Mass_cell_neighbor(fe_v,
															fe_v_neighbor,hp_fe_data_structure.dofs_per_component,
															hp_fe_data_structure.dofs_per_component,J);

				Mass_u2_v2 =  this->Compute_Mass_cell_neighbor(fe_v_neighbor,
																fe_v_neighbor,hp_fe_data_structure.dofs_per_component,
																hp_fe_data_structure.dofs_per_component,J);

				Mass_u1_v2 =  this->Compute_Mass_cell_neighbor(fe_v_neighbor,
															  fe_v,hp_fe_data_structure.dofs_per_component,
															  hp_fe_data_structure.dofs_per_component,J);

    // the number of rows always comes from the test function and the number of columns comes from the solution vector
				u1_v1 = matrix_opt.multiply_scalar(0.5,matrix_opt.compute_A_outer_B_limitA(Am,Mass_u1_v1,nEqn_this,nEqn_this));
				u2_v1 = matrix_opt.multiply_scalar(-0.5,matrix_opt.compute_A_outer_B_limitA(Am,Mass_u2_v1,nEqn_this,nEqn_neighbor));
				u1_v2 = matrix_opt.multiply_scalar(-0.5,matrix_opt.compute_A_outer_B_limitA(Am_neighbor,Mass_u1_v2,nEqn_neighbor,nEqn_this));
				u2_v2 = matrix_opt.multiply_scalar(0.5,matrix_opt.compute_A_outer_B_limitA(Am_neighbor,Mass_u2_v2,nEqn_neighbor,nEqn_neighbor));

			}

}
