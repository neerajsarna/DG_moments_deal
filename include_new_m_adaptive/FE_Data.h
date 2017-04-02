// in the present file we store the fe data object needed for single system computation and hp computations
namespace FEM_Solver
{
	using namespace dealii;

	template<int dim>
	class
	fe_data:public MeshGenerator::Base_MeshGenerator<dim>
	{
		public:
			fe_data(
							const std::string &output_file_name,
							const constant_numerics &constants,
						  const std::vector<int> &nEqn);

			FESystem<dim> finite_element;
            DoFHandler<dim> dof_handler;
            MappingQ<dim> mapping;
	};

	template<int dim>
	fe_data<dim>::fe_data(
							const std::string &output_file_name,
							const constant_numerics &constants,
						  const std::vector<int> &nEqn)
	:
	MeshGenerator::Base_MeshGenerator<dim>(output_file_name,constants),
	finite_element(FE_DGQ<dim>(constants.p),nEqn[0]),
	dof_handler(this->triangulation),
	mapping(constants.mapping_order)
	{
		AssertDimension(nEqn.size(),1);
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

                // fe_index of the maximum moment system which has to be solved in the computation
                const unsigned int max_fe_index;


                void construct_block_structure(std::vector<int> &block_structure,const std::vector<int> &nEqn);

                // using the block structure constructed in the previous routine, in the following routine we develop 
                // the fe collection object.
                void construct_fe_collection(const std::vector<int> &nEqn);
            
                // compute the deviation from the equilibrium value
               void compute_equilibrium_deviation(const int ngp,
                                                const std::vector<int> &nEqn,
                                                Triangulation<dim> &triangulation,
                                                Vector<double> &solution,
                                                const int cycle);

                // sets the tolerance band for difference moment systems
               // setting of the tolerance based upon the deviation of the distribution function from the local
               // maxwellian.
                void set_tolerance_bands_equilibrium_deviation(const int cycle);

               // setting of the tolerance based upon the distance of the cell from the center.
                void set_tolerance_bands_distance_center(const int cycle);

                void allocate_fe_index_equilibrium_deviation(const unsigned int present_cycle,
                											 const unsigned int total_cycles);

                void allocate_fe_index_distance_center(const unsigned int present_cycle,
                									   const unsigned int total_cycles);

                // compute the norm of the deviation of the distribution function from the local Maxwellian
                Vector<double> VelocitySpace_error_per_cell;


                // we store the number of degrees of freedom in every cell. The assumption is that all the components 
                // of all the moment system are solved with the same polynomial degree.
                std::vector<unsigned int> dofs_per_cell;



            // tolerance for the residual, if we find that the residual in a particular cell is greater than this value
            // then we refine it. 
                std::vector<double> VelocitySpace_error_tolerance;

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
    max_fe_index(nEqn.size()-1)
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


//setting of the tolerance bands based upon the deviation of the distribution function from the equilibrium
template<int dim>
void 
hp_fe_data<dim>::set_tolerance_bands_equilibrium_deviation(const int cycle)
{
	MatrixOpt::Base_MatrixOpt matrix_opt;
    // we only compute the tolerance bands for the 
    
    // size = max_nEqn + 1
    VelocitySpace_error_tolerance.resize(max_fe_index+2);

    //The maximum deviation obtained 
    const double max_error = matrix_opt.max_Vector(VelocitySpace_error_per_cell);

    // the minum deviation obtained
    const double min_error = matrix_opt.min_Vector(VelocitySpace_error_per_cell);
 
    // the first entry will just be the minimum value  and the last entry will we the maximum value
    VelocitySpace_error_tolerance[0] = min_error;

    // the fractions of the total difference for every band
    // size = nEqn
    std::vector<double> frac(max_fe_index + 1);

    // the last element will always be one
    frac[frac.size()-1] = 1.0;


    // equipartition for all the other elements
    // values for all the interior bands
    // for (unsigned long int i = 0 ; i < frac.size()-1; i ++)
    // {
    //     if (i == 0)
    //         frac[i] = 1.0/(max_fe_index+1);

    //     else
    //         frac[i] = 1.0/(max_fe_index+1) + frac[i-1];

    // }

    // adaptive fractions
    switch(max_fe_index)
    {
    	// 2 systems
    	// use a negative value if all of the cells should get the higher order moment
    	// use 2.0 if all the cells should get the lowest order moment
    	case 1:
    	{
    		frac[0] = 1 - (cycle * 0.5 + 1.0)/10;
    		std::cout << "********* fractions of 0 fe index " << frac[0] << std::endl;
    		break;
    	}
    	
    	// 3 systems
    	case 2:
    	{
    		frac[0] = 1 - (cycle * 0.5 + 1.0)/5;
    		frac[1] = 1 - (cycle * 0.5 + 1.0)/10;
    		break;
    	}

    	// 4 systems
    	case 3:
    	{
    		frac[0] = 1 - (cycle * 0.5 + 1.0)/2.5;
    		frac[1] = 1 - (cycle * 0.5 + 1.0)/5;
    		frac[2] = 1 - (cycle * 0.5 + 1.0)/10;
    		break;
    	}



    	default:
    	{
    		AssertThrow(1 == 0 ,ExcMessage("Should not have reached here"));
    		break;
    	}
    }
    // the value for the first one has alreayd been set
    for (unsigned long int i = 0 ; i < frac.size() ; i++)
        VelocitySpace_error_tolerance[i + 1] = min_error + frac[i] * (max_error-min_error);

}

// tolerance band based upon the distance of the cell from the center of the domain
template<>
void 
hp_fe_data<1>::set_tolerance_bands_distance_center(const int cycle)
{ 
	const double min_distance = 0;
	const double max_distance = 0.5;

    // size = max_nEqn + 1
    VelocitySpace_error_tolerance.resize(max_fe_index+2);

    // the first entry will just be the minimum value  and the last entry will we the maximum value
    VelocitySpace_error_tolerance[0] = min_distance;

    // the fractions of the total difference for every band
    // size = nEqn
    std::vector<double> frac(max_fe_index + 1);

    // the last element will always be one
    frac[frac.size()-1] = 1.0;


    // adaptive fractions
    switch(max_fe_index)
    {
    	// 2 systems
    	// use a negative value if all of the cells should get the higher order moment
    	// use 2.0 if all the cells should get the lowest order moment
    	case 1:
    	{
    		frac[0] = 1 - (cycle * 0.5 + 1.0)/10;
    		std::cout << "********* fractions of 0 fe index " << frac[0] << std::endl;
    		break;
    	}
    	
    	// 3 systems
    	case 2:
    	{
    		frac[0] = 1 - (cycle * 0.5 + 1.0)/5;
    		frac[1] = 1 - (cycle * 0.5 + 1.0)/10;
    		break;
    	}

    	// 4 systems
    	case 3:
    	{
    		frac[0] = 1 - (cycle * 0.5 + 1.0)/2.5;
    		frac[1] = 1 - (cycle * 0.5 + 1.0)/5;
    		frac[2] = 1 - (cycle * 0.5 + 1.0)/10;
    		break;
    	}



    	default:
    	{
    		AssertThrow(1 == 0 ,ExcMessage("Should not have reached here"));
    		break;
    	}
    }
    // the value for the first one has alreayd been set
    for (unsigned long int i = 0 ; i < frac.size() ; i++)
        VelocitySpace_error_tolerance[i + 1] = min_distance + frac[i] * (max_distance-min_distance);

}



// in the following function we compute the deviation of the distribution function from the equilibrium
template<int dim>
void 
hp_fe_data<dim>::compute_equilibrium_deviation(const int ngp,
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
        

	// total number of conserved variables = dim (velocity) + rho + theta.
	// -1 to accommodate for the c indices
	const unsigned int ID_num_conserved = dim + 1;
	std::pair<unsigned int,unsigned int> selected;

    // the maximum number of equations present in the moment system
	const int max_equations = nEqn[max_fe_index];

	// zero weightage to all the quantities which are conserved
    // The deviation of the distribution function from the local Maxwellian is only modelled by the higher
    // order moments. Therefore while computing the residuals we do not consider any contribution from 
    // the higher order moments. We mask all the conserved quantities
	selected = std::make_pair(ID_num_conserved+1,max_equations);

    ComponentSelectFunction<dim> weight(selected,max_equations);          // used to compute only the error in theta


    // compute the l1 norm of the residual
    // we first compute the L2 norm per cell. Lets assume that the higher order contribution is being denoted by 
    // only the second order moment R_{ij}. Then the following function will compute \sqrt(\int R_{ij}dx) on the domain
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
        VelocitySpace_error_per_cell(counter) = pow(VelocitySpace_error_per_cell(counter),2)/cell->measure();
        counter ++;
    }


    // now we set the limits for the tolerance
    set_tolerance_bands_distance_center(cycle);
        
    
}


// allocation of the fe index based upon the deviation of the distribution function from the equilibrium.
template<int dim>
void 
hp_fe_data<dim>::allocate_fe_index_equilibrium_deviation(const unsigned int present_cycle,const unsigned int total_cycles)
{
	typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
	unsigned int counter = 0;
	const unsigned int initial_fe_index = 0;
	Vector<double> refined(max_fe_index+1);
	refined = 0;

	// the very first cycle of m refinement
	 if (present_cycle == 0)
	 {
		for (; cell != endc ; cell++)
		{
			refined(initial_fe_index)++;
			cell->set_active_fe_index(initial_fe_index);
		}

		cell = dof_handler.begin_active();

	}


	else
	{

		for (; cell != endc ; cell++)
		{
			const double error_value = VelocitySpace_error_per_cell(counter);

			for (unsigned int i = 0 ; i < VelocitySpace_error_tolerance.size()-1 ; i++)
				if (error_value >= VelocitySpace_error_tolerance[i] && error_value <= VelocitySpace_error_tolerance[i+1])
				{
					refined(i)++;
					cell->set_active_fe_index(i);
				}
 
			counter++;

		}		
	}


	for (unsigned long int i = 0 ; i < refined.size() ; i++)
		std::cout << "Fe Index: " << i << " Times: " << refined(i) << std::endl;



}


template<int dim>
void 
hp_fe_data<dim>::allocate_fe_index_distance_center(const unsigned int present_cycle,const unsigned int total_cycles)
{
	typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
	unsigned int counter = 0;
	const unsigned int initial_fe_index = 0;
	Vector<double> refined(max_fe_index+1);
	refined = 0;

	// the very first cycle of m refinement
	 if (present_cycle == 0)
	 {
		for (; cell != endc ; cell++)
		{
			refined(initial_fe_index)++;
			cell->set_active_fe_index(initial_fe_index);
		}

		cell = dof_handler.begin_active();

	}


	else
	{

		for (; cell != endc ; cell++)
		{
			const double distance_center = cell->center().norm();

			for (unsigned int i = 0 ; i < VelocitySpace_error_tolerance.size()-1 ; i++)
				if (distance_center >= VelocitySpace_error_tolerance[i] && 
					distance_center <= VelocitySpace_error_tolerance[i+1])
				{
					refined(i)++;
					cell->set_active_fe_index(i);
				}
 
			counter++;

		}		
	}


	for (unsigned long int i = 0 ; i < refined.size() ; i++)
		std::cout << "Fe Index: " << i << " Times: " << refined(i) << std::endl;



}



}