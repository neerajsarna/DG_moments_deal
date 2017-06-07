// in the present file we store the fe data object needed for single system computation and hp computations
namespace FEM_Solver
{
	using namespace dealii;

    // the data structures needed for finite element computation
	template<int dim>
	class
	fe_data:public MeshGenerator::Base_MeshGenerator<dim>
	{
		public:
			fe_data(
							const std::string &output_file_name,
							const constant_numerics &constants,
						  const int nEqn);

			FESystem<dim> finite_element;
            DoFHandler<dim> dof_handler;
            MappingQ<dim> mapping;
	};

	template<int dim>
	fe_data<dim>::fe_data(
							const std::string &output_file_name,
							const constant_numerics &constants,
						  const int nEqn)
	:
	MeshGenerator::Base_MeshGenerator<dim>(output_file_name,constants),
	finite_element(FE_DGQ<dim>(constants.p),nEqn),
	dof_handler(this->triangulation),
	mapping(constants.mapping_order)
	{
	}

	template<int dim>
	class
	hp_fe_data:public MeshGenerator::Base_MeshGenerator<dim>
	{
		public:
			hp_fe_data(
							const std::string &output_file_name,
							const constant_numerics &constants,
						  const std::vector<int> &nEqn);

			    FE_DGQ<dim> fe_dg;

                // block sizes for the finite element objects
                std::vector<int> block_finite_elements;

                // the main finite element object
                hp::FECollection<dim> finite_element;

                // the hp dof handler object
                hp::DoFHandler<dim> dof_handler;

                MappingQ<dim,dim> mapping_basic;
                
                hp::MappingCollection<dim,dim> mapping;

                MatrixOpt::Base_MatrixOpt matrix_opt;

                // fe_index of the maximum moment system which has to be solved in the computation
                const unsigned int max_fe_index;
                const unsigned int max_equations;


                void construct_block_structure(std::vector<int> &block_structure,const std::vector<int> &nEqn);

                // using the block structure constructed in the previous routine, in the following routine we develop 
                // the fe collection object.
                void construct_fe_collection(const std::vector<int> &nEqn);
            
               void compute_distribution_deviation(const int ngp,
                                                const std::vector<int> &nEqn,
                                                Triangulation<dim> &triangulation,
                                                Vector<double> &solution,
                                                const int cycle);

                void compute_error_comparitive(const Vector<double> &solution,
                                               const constant_numerics &constants,
                                               const DoFHandler<dim> &dof_handler_reference,
                                               const Vector<double> &solution_reference);

              // we set the tolerance depending upon which moment theory has to be adapted
                double set_tolerance_distribution_deviation(const unsigned int cycle);

                std::vector<double> set_tolerance_error_comparison();

                void allocate_fe_index_distribution_deviation(const unsigned int present_cycle);

                // allocation of fe index based upon the error obtained by comparision with a higher
                // order moment method
                void allocate_fe_index_error_comparison(const unsigned int present_cycle);

                // compute the norm of the deviation of the distribution function from the local Maxwellian
                Vector<double> VelocitySpace_error_per_cell;


                // we store the number of degrees of freedom in every cell. The assumption is that all the components 
                // of all the moment system are solved with the same polynomial degree.
                std::vector<unsigned int> dofs_per_cell;


            // the total number of degrees of freedom per component will remain the same for every cell since we do not
            // have p-adaptivity presently
                unsigned int dofs_per_component;

	};

	
	template<int dim>
	hp_fe_data<dim>::hp_fe_data(
							const std::string &output_file_name,
							const constant_numerics &constants,
						  const std::vector<int> &nEqn)
	:
	MeshGenerator::Base_MeshGenerator<dim>(output_file_name,constants),
	fe_dg(constants.p),
	dof_handler(this->triangulation),
	mapping_basic(constants.mapping_order),
    max_fe_index(nEqn.size()-1),
    max_equations(nEqn[max_fe_index])
    {
    		// construct the finite element collection for the present class
    	     construct_fe_collection(nEqn);

        	// the total number of finite element objects must be equal to the total number of systems sent to the solver
        	AssertDimension(nEqn.size(),finite_element.size());

        	// now we initialize the hp mapping object which simply contain the same mapping object repeated a several times
        	for (unsigned long int i = 0 ;i < nEqn.size(); i++)
            	mapping.push_back(mapping_basic);

        // we now initialize the dofs per cell and the numer of dofs per component. We can do this since we are looking for 
        // the same polynomial degree for every moment system
              // for every moment system, the number of dofs in a given cell are different
      		for (unsigned long int i = 0 ;i < nEqn.size() ; i++)
        		dofs_per_cell.push_back(finite_element[i].dofs_per_cell);


      	// these are the number of degrees of freedom used for one component. They remain the same for all the equations
        // the assumption is that the polynomial degree remains the same for all the cells. So we do not have p adaptivity
      	dofs_per_component = dofs_per_cell[0]/nEqn[0];


      	// the vector dofs per cell should be sorted since the vector number of equations is sorted
      	Assert(std::is_sorted(std::begin(dofs_per_cell),
      							std::end(dofs_per_cell)),ExcMessage("The number of dofs are not sorted"));

    }




template<int dim>
void 
hp_fe_data<dim>::construct_block_structure(std::vector<int> &block_structure,const std::vector<int> &nEqn)
{
		Assert(nEqn.size() != 0, ExcNotInitialized());
		Assert(std::is_sorted(std::begin(nEqn),std::end(nEqn)),ExcMessage("number of equations not sorted"));

		// the very first entry should be the number of equations in the first system
		block_structure.push_back(nEqn[0]);

		for (unsigned long int i = 1 ; i < nEqn.size() ; i++)
			block_structure.push_back(nEqn[i]-nEqn[i-1]);

		AssertDimension(block_structure.size(),nEqn.size());
}

template<int dim>
void 
hp_fe_data<dim>::construct_fe_collection(const std::vector<int> &nEqn)
{
		std::vector<int> block_structure;

		// first create the block structure for the finite element object which will then be used to construct the fe system
		construct_block_structure(block_structure,nEqn);

		switch (block_structure.size())
		{
			case 1:
			{
				finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0]));
				break;
			}

			case 2:
			{
			
				finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
													FE_Nothing<dim>(),block_structure[1]));

				finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
													fe_dg,block_structure[1]));

				break;
			}

			case 3:
			{
							finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
													FE_Nothing<dim>(),block_structure[1],
													FE_Nothing<dim>(),block_structure[2]));

			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
												   fe_dg,block_structure[1],
												   FE_Nothing<dim>(),block_structure[2]));

			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
												   fe_dg,block_structure[1],
												   fe_dg,block_structure[2]));
				break;
			}

			case 4:
			{
							finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
												   FE_Nothing<dim>(),block_structure[1],
							 						FE_Nothing<dim>(),block_structure[2],
							 						FE_Nothing<dim>(),block_structure[3]));

			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
				                                   fe_dg,block_structure[1],
							                       FE_Nothing<dim>(),block_structure[2],
							                       FE_Nothing<dim>(),block_structure[3]));

			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
												   fe_dg,block_structure[1],
							 					   fe_dg,block_structure[2],
							 					   FE_Nothing<dim>(),block_structure[3]));

			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
				                                   fe_dg,block_structure[1],
							 					   fe_dg,block_structure[2],
							 					   fe_dg,block_structure[3]));
				break;
			}

			case 5:
			{
			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
												   FE_Nothing<dim>(),block_structure[1],
							 						FE_Nothing<dim>(),block_structure[2],
							 						FE_Nothing<dim>(),block_structure[3],
							 						FE_Nothing<dim>(),block_structure[4]));

			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
												   fe_dg,block_structure[1],
							 						FE_Nothing<dim>(),block_structure[2],
							 						FE_Nothing<dim>(),block_structure[3],
							 						FE_Nothing<dim>(),block_structure[4]));

			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
												   fe_dg,block_structure[1],
							 						fe_dg,block_structure[2],
							 						FE_Nothing<dim>(),block_structure[3],
							 						FE_Nothing<dim>(),block_structure[4]));

			finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
												   fe_dg,block_structure[1],
							 						fe_dg,block_structure[2],
							 						fe_dg,block_structure[3],
							 						FE_Nothing<dim>(),block_structure[4]));

		finite_element.push_back(FESystem<dim>(fe_dg,block_structure[0],
												   fe_dg,block_structure[1],
							 						fe_dg,block_structure[2],
							 						fe_dg,block_structure[3],
							 						fe_dg,block_structure[4]));
				break;
			}

			default:
			{
				AssertThrow(1 == 0, ExcMessage("Should not have reached here"));
				break;
			}
		}

			
}

// deviation of the distribution function from the previous one
template<int dim>
void 
hp_fe_data<dim>::compute_distribution_deviation(const int ngp,
                                                const std::vector<int> &nEqn,
                                                Triangulation<dim> &triangulation,
                                                Vector<double> &solution,
                                                const int cycle)
{

      const QGauss<dim> quadrature_basic(ngp);
      // integration on the volume
      hp::QCollection<dim> quadrature;

      for (unsigned long int i = 0 ; i < nEqn.size() ; i++)
        quadrature.push_back(quadrature_basic);
        

    // ID of the moments from where we need to mask the elements
    unsigned int ID_start;

    // end ID till where we need to mask the moments
    unsigned int ID_end;

    std::pair<unsigned int,unsigned int> selected;

    if (cycle == 0)
    {
        ID_start = dim + 2;     // number of conserved quantities
        ID_end = nEqn[0];    // component ID of the present system
    }
    else
    {
        ID_start = nEqn[cycle-1];   // total equations of the previous cycle
        ID_end = nEqn[cycle];     // total equations in the present cycle
    }


    AssertIndexRange(ID_start,max_equations + 1);
    AssertIndexRange(ID_end,max_equations + 1);

    // a pair containing the IDs of the moments which should be computed [ID_start,ID_end)
    selected = std::make_pair(ID_start,ID_end); 

    ComponentSelectFunction<dim> weight(selected,max_equations);          // used to compute only the error in theta

    VectorTools::integrate_difference (mapping,
                                       dof_handler,
                                       solution,
                                      ZeroFunction<dim>(max_equations),
                                      VelocitySpace_error_per_cell,
                                      quadrature,
                                      VectorTools::L2_norm,
                                      &weight);  


    // we now convert from the l2 norm to the l1 norm. We would like to compute the l1 norm of the deviation and therefore
    // we simply square the value obtained above. The reason for this is that since we have already taken the square 
    // therefore the function is anyhow positive. Therefore to compute the l1 norm we just need to square the expression 
    // obtained from the above expression. 
    typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), 
                                                      endc = triangulation.end();

    unsigned int counter = 0;
    for (; cell != endc ; cell++)
    {
        // by taking the square, we are converting back into the l1 norm. Divide by the cell measure to normalize
        // with respect to the area of the cell.
        VelocitySpace_error_per_cell(counter) = pow(VelocitySpace_error_per_cell(counter),2)/cell->measure();
        counter ++;
    }

            
}

template<int dim>
void 
hp_fe_data<dim>::compute_error_comparitive(const Vector<double> &solution,
                                          const constant_numerics &constants,
                                          const DoFHandler<dim> &dof_handler_reference,
                                          const Vector<double> &solution_reference)
{
      const QGauss<dim> quadrature_basic(constants.p + 1);

      // integration on the volume
      hp::QCollection<dim> quadrature;

      for (unsigned long int i = 0 ; i <= max_fe_index ; i++)
        quadrature.push_back(quadrature_basic);


      const UpdateFlags update_flags               =  update_gradients
                                                     | update_q_points
                                                     | update_JxW_values
                                                     | update_values,

      face_update_flags          = update_values
      | update_q_points
      | update_JxW_values
      | update_normal_vectors,
      neighbor_face_update_flags = update_values;

      hp::FEValues<dim>  hp_fe_v(this->mapping,
                                 this->finite_element,
                                 quadrature,
                                 update_flags);


      typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(),
                                                         endc = this->dof_handler.end();


      const int total_ngp = quadrature_basic.size();
      std::vector<double> Jacobians_interior(total_ngp);
      std::vector<Point<dim>> Quad_points(total_ngp);
      unsigned int counter = 0;

      unsigned int component;
        if (dim == 1)
            component = constants.variable_map_1D.find(constants.error_variable)->second;

        if (dim == 2)
            component = constants.variable_map.find(constants.error_variable)->second;


      for (; cell != endc ; cell++)
      {
                hp_fe_v.reinit(cell);
                const FEValues<dim> &fe_v = hp_fe_v.get_present_fe_values();

                // we now extract the quadrature points and the quadrature weights
                Jacobians_interior = fe_v.get_JxW_values();
                Quad_points = fe_v.get_quadrature_points();

                // we now loop over all the quadrature points
                for (unsigned int q = 0 ; q < total_ngp ; q++)
                {
                    Vector<double> solution_value(max_equations);
                    Vector<double> reference_value(max_equations);

                    VectorTools::point_value(this->dof_handler, solution, Quad_points[q],solution_value); 
                    VectorTools::point_value(dof_handler_reference,
                                            solution_reference, Quad_points[q],reference_value); 


                    const double error_value = fabs(solution_value(component)-reference_value(component));

                    // value of the function multiplied by the jacobian and the weights
                    VelocitySpace_error_per_cell(counter) += pow(error_value,2) * Jacobians_interior[q];
                }

                // convert back into the L2 norm
                VelocitySpace_error_per_cell(counter) = sqrt(VelocitySpace_error_per_cell(counter));

                counter ++;

      }
}



template<int dim>
double
hp_fe_data<dim>::set_tolerance_distribution_deviation(const unsigned int cycle)
{
        typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), 
                                                           endc = dof_handler.end();

        // the error in the velocity space corresponding to the highest present fe index
        std::vector<double> error_VelocitySpace_max_theory;
        int counter = 0;
        int count_lower = 0;

        for (; cell != endc ; cell++)
        {
            if (cell->active_fe_index() == cycle - 1)  // if the fe index is equal to the highest fe index present, then only we store
            {
                error_VelocitySpace_max_theory.push_back(VelocitySpace_error_per_cell(counter));
                count_lower ++ ;
            }

            counter++;            
        }

        // all the cells which have a deviation > min + frac(max - min) will be refined

        const double frac = this->constants.fraction_refine;
        const double max_error = matrix_opt.max_Vector(error_VelocitySpace_max_theory);
        const double min_error = matrix_opt.min_Vector(error_VelocitySpace_max_theory);
        const double tolerance = min_error + frac * (max_error - min_error) ;


        return(tolerance);
}

template<int dim>
std::vector<double>
hp_fe_data<dim>::set_tolerance_error_comparison()
{

        const double max_error = matrix_opt.max_Vector(VelocitySpace_error_per_cell);
        const double min_error = matrix_opt.min_Vector(VelocitySpace_error_per_cell);

        // size is total number of systems + 1
        std::vector<double> tolerance(max_fe_index+2);

        tolerance[0] = min_error;
        tolerance[max_fe_index + 1] = max_error;

        // equipartition the error while creating tolerance bands
        for (unsigned int i = 1 ; i < max_fe_index + 1 ; i++)
            tolerance[i] = tolerance[i - 1] + (max_error-min_error)/(max_fe_index+1);

        return(tolerance);
}


template<int dim>
void 
hp_fe_data<dim>::allocate_fe_index_distribution_deviation(const unsigned int present_cycle)
{
        typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), 
                                                           endc = dof_handler.end();   


        int counter = 0;
        std::vector<int> fe_index_count(max_fe_index+1);

        if (present_cycle == 0)
            for(; cell != endc ; cell++)
            {
                cell->set_active_fe_index(0);
                fe_index_count[cell->active_fe_index()]++;
            }

        else
        {
            const double tolerance = set_tolerance_distribution_deviation(present_cycle);

            for(; cell != endc ; cell++)
            {

            if (cell->active_fe_index()  == present_cycle - 1) // only increase the fe index if we are on a lower order moment theory
                if ((VelocitySpace_error_per_cell(counter) - tolerance) >= 10e-100 )
                {
                    cell->set_active_fe_index(present_cycle);   // increase the fe index if the residual is too high
                }

                fe_index_count[cell->active_fe_index()]++;      // keep a count of the number of systems in the domain
                counter ++;
            }

        }



    for (unsigned long int i = 0 ; i < fe_index_count.size() ; i++)
        std::cout << "Fe Index: " << i << " Times: " << fe_index_count[i] << std::endl;

}

template<int dim>
void 
hp_fe_data<dim>::allocate_fe_index_error_comparison(const unsigned int present_cycle)
{
        typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), 
                                                           endc = dof_handler.end();   


        int counter = 0;
        std::vector<int> fe_index_count(max_fe_index+1);

        if (present_cycle == 0)
            for(; cell != endc ; cell++)
            {
                cell->set_active_fe_index(0);
                fe_index_count[cell->active_fe_index()]++;
            }

        else
        {
            const std::vector<double> tolerance = set_tolerance_error_comparison();

            AssertDimension(tolerance.size(),max_equations);
            for(; cell != endc ; cell++)
             for (unsigned int i = 0 ; i < tolerance.size() - 1 ; i++)
                if ( fabs(VelocitySpace_error_per_cell(counter) - tolerance[i]) >= 1e-16 &&
                     fabs(VelocitySpace_error_per_cell(counter) - tolerance[i + 1]) <= 1e-16)
                {
                    AssertIndexRange(i,max_fe_index);
                    cell->set_active_fe_index(i);
                    fe_index_count[cell->active_fe_index()]++;      // keep a count of the number of systems in the domain
                    counter ++;
                }


            

        }



    for (unsigned long int i = 0 ; i < fe_index_count.size() ; i++)
        std::cout << "Fe Index: " << i << " Times: " << fe_index_count[i] << std::endl;

}

}
