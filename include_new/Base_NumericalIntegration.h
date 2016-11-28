namespace NumericalIntegration
{
	using namespace dealii;

	template<int dim>
	class 
	Base_NumericalIntegration
	{
			public:
				Base_NumericalIntegration();
				Full_matrix shape_values;

				// gauss points are the number of gauss points in a particular direction
				void Compute_Shape_Value(MappingQ<dim> &mapping,
										 const unsigned int p,
										 const unsigned int gauss_points,
										 typename DoFHandler<dim>::active_cell_iterator &cell);

				// gauss points are the number of gauss points in one particular direction
				// compute the mass matrix M such that M_{ij} = \int \phi_i \pd_x\phi_i dx
				Full_matrix Compute_Mass_shape_grad_x(FEValuesBase<dim> &fe_v,
							 							const unsigned int dofs_per_component,
							  							std::vector<double> &J,
														const typename DoFHandler<dim>::active_cell_iterator &cell);

				Full_matrix Compute_Mass_shape_grad_y(FEValuesBase<dim> &fe_v,
														const unsigned int dofs_per_component,
														std::vector<double> &J,
										 				const typename DoFHandler<dim>::active_cell_iterator &cell);

				Full_matrix Compute_Mass_shape_grad_z(FEValuesBase<dim> &fe_v,
														const unsigned int dofs_per_component,
														std::vector<double> &J,
										 				const typename DoFHandler<dim>::active_cell_iterator &cell);

				Full_matrix Compute_Mass_shape_value(FEValuesBase<dim> &fe_v,
														const unsigned int dofs_per_component,
														std::vector<double> &J,
										 				const typename DoFHandler<dim>::active_cell_iterator &cell);

				// compute the mass matrix for the edge share between the cell and the neighbor
				Full_matrix Compute_Mass_cell_neighbor(FEValuesBase<dim> &fe_v1,
													   FEValuesBase<dim> &fe_v2,
													   const unsigned int dofs_per_component1,
													   const unsigned int dofs_per_component2,
													   std::vector<double> &J);
	};

	template<int dim>
	Base_NumericalIntegration<dim>::
	Base_NumericalIntegration()
	{;}

	template<int dim>
	void
	Base_NumericalIntegration<dim>::
	Compute_Shape_Value(MappingQ<dim> &mapping,
						const unsigned int p,
						const unsigned int gauss_points,
						typename DoFHandler<dim>::active_cell_iterator &cell)
	{
		QGauss<dim> quadrature(gauss_points);
		const unsigned int total_ngp = quadrature.size();

		// construct the finite element corresponding to just one component
		const FiniteElement<dim> &fe_in_cell = cell->get_fe();
		const unsigned int n_components = fe_in_cell.n_components();
		const unsigned int dofs_per_component = fe_in_cell.dofs_per_cell/n_components;

		// we initialize fe_values because we call this from outside of the assemblation loop
		UpdateFlags update_flags = update_values;
		FEValues<dim> fe_v(mapping,fe_in_cell,quadrature, update_flags);
		fe_v.reinit(cell);

		shape_values.resize(dofs_per_component,total_ngp);

		Assert(shape_values.rows() !=0,ExcNotInitialized());

		for (unsigned int i = 0 ; i < dofs_per_component ; i ++)
			for (unsigned int q = 0 ; q < total_ngp ; q++)
			{
				shape_values(i,q) = 0;
				shape_values(i,q) = fe_v.shape_value(i,q);
			}
	}

	template<int dim>
	Full_matrix
	Base_NumericalIntegration<dim>::
	Compute_Mass_shape_grad_x(FEValuesBase<dim> &fe_v,
							 const unsigned int dofs_per_component,
							  std::vector<double> &J,
							const typename DoFHandler<dim>::active_cell_iterator &cell)
	{

		const unsigned int total_ngp = J.size();	
		const unsigned int num_dof_per_comp = fe_v.get_fe().dofs_per_cell/fe_v.get_fe().n_components();		
		Full_matrix M(dofs_per_component,dofs_per_component);
		M.setZero();

		Assert(shape_values.rows() || shape_values.cols() != 0,ExcNotInitialized());
 		AssertDimension(dofs_per_component,shape_values.rows());
 		AssertDimension(num_dof_per_comp,dofs_per_component);
 		
 		for (unsigned int q = 0 ; q < total_ngp ; q++)
 		{
 		const double jacobian_value = J[q];
		for (unsigned int i = 0 ; i < dofs_per_component ; i ++)
			for (unsigned int j = 0 ; j < dofs_per_component ; j++)	
					M.coeffRef(i,j) += shape_values(i,q) * fe_v.shape_grad(j,q)[0] * jacobian_value;
 		}


		return(M);

	}


	template<int dim>
	Full_matrix
	Base_NumericalIntegration<dim>::
	Compute_Mass_shape_grad_y(FEValuesBase<dim> &fe_v,
							 const unsigned int dofs_per_component,
							  std::vector<double> &J,
							const typename DoFHandler<dim>::active_cell_iterator &cell)
	{

		const unsigned int total_ngp = J.size();	
		const unsigned int num_dof_per_comp = fe_v.get_fe().dofs_per_cell/fe_v.get_fe().n_components();				
		Full_matrix M(dofs_per_component,dofs_per_component);
		M.setZero();

		Assert(shape_values.rows() || shape_values.cols() != 0,ExcNotInitialized());
 		AssertDimension(dofs_per_component,shape_values.rows());
 		AssertDimension(num_dof_per_comp,dofs_per_component);
 		
 		for (unsigned int q = 0 ; q < total_ngp ; q++)
 		{
 			const double jacobian_value = J[q];
			for (unsigned int i = 0 ; i < dofs_per_component ; i ++)
				for (unsigned int j = 0 ; j < dofs_per_component ; j++)
						M.coeffRef(i,j) += shape_values(i,q) * fe_v.shape_grad(j,q)[1] * jacobian_value; 			
 		}


		return(M);
	}

	template<int dim>
	Full_matrix
	Base_NumericalIntegration<dim>::
	Compute_Mass_shape_grad_z(FEValuesBase<dim> &fe_v,
							 const unsigned int dofs_per_component,
							  std::vector<double> &J,
							const typename DoFHandler<dim>::active_cell_iterator &cell)
	{

		const unsigned int total_ngp = J.size();
		const unsigned int num_dof_per_comp = fe_v.get_fe().dofs_per_cell/fe_v.get_fe().n_components();				
		Full_matrix M(dofs_per_component,dofs_per_component);
		M.setZero();

		Assert(shape_values.rows() || shape_values.cols() != 0,ExcNotInitialized());
  		AssertDimension(dofs_per_component,shape_values.rows());
  		AssertDimension(num_dof_per_comp,dofs_per_component);

  		for (unsigned int q = 0 ; q < total_ngp ; q++)
  		{
  		const double jacobian_value = J[q];
		for (unsigned int i = 0 ; i < dofs_per_component ; i ++)
			for (unsigned int j = 0 ; j < dofs_per_component ; j++)	
					M(i,j) += shape_values(i,q) * fe_v.shape_grad(j,q)[2] * jacobian_value;
  		}



		return(M);

	}
	
	// now we compute the mass matrix corresponding to Integrate phi * phi dx
	template<int dim>
	Full_matrix
	Base_NumericalIntegration<dim>::
	Compute_Mass_shape_value(FEValuesBase<dim> &fe_v,
							 const unsigned int dofs_per_component,
							  std::vector<double> &J,
							const typename DoFHandler<dim>::active_cell_iterator &cell)
	{

		const unsigned int total_ngp = J.size();
		const unsigned int num_dof_per_comp = fe_v.get_fe().dofs_per_cell/fe_v.get_fe().n_components();				
		Full_matrix M(dofs_per_component,dofs_per_component);
		M.setZero();

		Assert(shape_values.rows() || shape_values.cols() != 0,ExcNotInitialized());
  		AssertDimension(dofs_per_component,shape_values.rows());
  		AssertDimension(num_dof_per_comp,dofs_per_component);

  		for (unsigned int q = 0 ; q < total_ngp ; q++)
  		{
  		const double jacobian_value = J[q];
			for (unsigned int i = 0 ; i < dofs_per_component ; i ++)
				for (unsigned int j = 0 ; j < dofs_per_component ; j++)	
						M(i,j) += shape_values(i,q) * shape_values(j,q) * jacobian_value;
  		}



		return(M);

	}

	// the following function computes the matrix M = basis of fe_v1 * basis of fe_v2 dx
	template<int dim>
	Full_matrix
	Base_NumericalIntegration<dim>::
	Compute_Mass_cell_neighbor(FEValuesBase<dim> &fe_v1,
							   FEValuesBase<dim> &fe_v2,
						       const unsigned int dofs_per_component1,
							   const unsigned int dofs_per_component2,
							   std::vector<double> &J)
	{
		const unsigned int total_ngp = J.size();
		const unsigned int num_dof_per_comp1 = fe_v1.get_fe().dofs_per_cell/fe_v1.get_fe().n_components();				
		const unsigned int num_dof_per_comp2 = fe_v2.get_fe().dofs_per_cell/fe_v2.get_fe().n_components();				

		Full_matrix M(dofs_per_component1,dofs_per_component2);
		M.setZero();

  		AssertDimension(num_dof_per_comp1,dofs_per_component1);
  		AssertDimension(num_dof_per_comp2,dofs_per_component2);

  		for (unsigned int q = 0 ; q < total_ngp ; q++)
  		{
  		const double jacobian_value = J[q];
			for (unsigned int i = 0 ; i < dofs_per_component1 ; i ++)
				for (unsigned int j = 0 ; j < dofs_per_component2 ; j++)	
						M(i,j) += fe_v1.shape_value(i,q) * fe_v2.shape_value(j,q) * jacobian_value;
  		}



		return(M);		
	}

}