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
				void Compute_Shape_Value(const MappingQ<dim> &mapping,
										 const unsigned int gauss_points,
										 typename DoFHandler<dim>::active_cell_iterator &cell);

				void Compute_Shape_Value(hp::FEValues<dim> hp_fe_v,const unsigned int dofs_per_component);


				FullMatrix<double>  Compute_Mass_shape_value(const FEValuesBase<dim> &fe_v,
														const unsigned int dofs_per_component,
														const std::vector<double> &J);

				std::vector<FullMatrix<double>> Compute_Mass_shape_grad(const FEValuesBase<dim> &fe_v,
														const unsigned int dofs_per_component,
														const std::vector<double> &J);

				// compute the mass matrix for the edge share between the cell and the neighbor
				FullMatrix<double> Compute_Mass_cell_neighbor(const FEValuesBase<dim> &fe_v1,
													   const FEValuesBase<dim> &fe_v2,
													   const unsigned int dofs_per_component1,
													   const unsigned int dofs_per_component2,
													   const std::vector<double> &J);

				MatrixOpt::Base_MatrixOpt matrix_opt;
	};

	template<int dim>
	Base_NumericalIntegration<dim>::
	Base_NumericalIntegration()
	{;}

	template<int dim>
	void
	Base_NumericalIntegration<dim>::
	Compute_Shape_Value(const MappingQ<dim> &mapping,
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
	
	// we precompute the values of the shape functions at the various gauss points
	template<int dim>
	void
	Base_NumericalIntegration<dim>::
	Compute_Shape_Value(hp::FEValues<dim> hp_fe_v,const unsigned int dofs_per_component)
	{

		// we get the fe value associated with the present hp_fe_v
		const FEValues<dim> &fe_v = hp_fe_v.get_present_fe_values();

		unsigned int total_ngp = fe_v.get_quadrature().size();
		shape_values.resize(dofs_per_component,total_ngp);


		for (unsigned int i = 0 ; i < dofs_per_component ; i ++)
			for (unsigned int q = 0 ; q < fe_v.get_quadrature().size() ; q++)
			{
				shape_values(i,q) = 0;
				shape_values(i,q) = fe_v.shape_value(i,q);
			}
	}


	template<int dim>
	std::vector<FullMatrix<double>>
	Base_NumericalIntegration<dim>::
	Compute_Mass_shape_grad(const FEValuesBase<dim> &fe_v,
							 const unsigned int dofs_per_component,
							  const std::vector<double> &J)
	{

		const unsigned int total_ngp = J.size();
		std::vector<FullMatrix<double>> M;

		for (int space = 0 ;space < dim ; space++)
		{
			M.push_back(FullMatrix<double>(dofs_per_component,dofs_per_component));
			Assert(M[space].m() != 0,ExcNotInitialized());
			Assert(M[space].n() != 0,ExcNotInitialized());
			M[space] = 0;
		}

	
		AssertDimension(shape_values.rows(),dofs_per_component);
		AssertDimension(shape_values.cols(),total_ngp);

  		// the assumption in the following computation is that all the moments have the same order of 
  		// polynomial degree
  		for (int space = 0 ;space < dim ; space++)
  		{
  			for (unsigned int q = 0 ; q < total_ngp ; q++)
  			{
  				const double jacobian_value = J[q];
  				for (unsigned int i = 0 ; i < dofs_per_component ; i ++)
  					for (unsigned int j = 0 ; j < dofs_per_component ; j++)	
  						M[space](i,j) += shape_values(i,q) * fe_v.shape_grad(j,q)[space] * jacobian_value;
  			}  			
  		}




		return(M);

	}
	// now we compute the mass matrix corresponding to Integrate phi * phi dx
	template<int dim>
	FullMatrix<double>
	Base_NumericalIntegration<dim>::
	Compute_Mass_shape_value(const FEValuesBase<dim> &fe_v,
							 const unsigned int dofs_per_component,
							  const std::vector<double> &J)
	{

		const unsigned int total_ngp = J.size();
		FullMatrix<double> M(dofs_per_component,dofs_per_component);
		M = 0;

		AssertDimension(shape_values.rows(),dofs_per_component);
		AssertDimension(shape_values.cols(),total_ngp);

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
	FullMatrix<double>
	Base_NumericalIntegration<dim>::
	Compute_Mass_cell_neighbor(const FEValuesBase<dim> &fe_v1,
							   const FEValuesBase<dim> &fe_v2,
						       const unsigned int dofs_per_component1,
							   const unsigned int dofs_per_component2,
							   const std::vector<double> &J)
	{
		const unsigned int total_ngp = J.size();

		// we assume that all the cells have the same polynomial degree
		AssertDimension(dofs_per_component1,dofs_per_component2);

		FullMatrix<double> M(dofs_per_component1,dofs_per_component2);
		M = 0;

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