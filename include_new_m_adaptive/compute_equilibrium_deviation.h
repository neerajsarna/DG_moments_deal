// in the following function we try to compute the tolerance bands for the differen theoreis
template<int dim>
void 
Base_Solver<dim>::set_tolerance_bands()
{
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
    for (unsigned long int i = 0 ; i < frac.size()-1; i ++)
    {
        if (i == 0)
            frac[i] = 1.0/(max_fe_index+1);

        else
            frac[i] = 1.0/(max_fe_index+1) + frac[i-1];

    }

    // the value for the first one has alreayd been set
    for (unsigned long int i = 0 ; i < frac.size() ; i++)
        VelocitySpace_error_tolerance[i + 1] = min_error + frac[i] * (max_error-min_error);

}

// in the following function we compute the deviation of the distribution function from the equilibrium
template<int dim>
void 
Base_Solver<dim>::compute_equilibrium_deviation()
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
    typename Triangulation<dim>::active_cell_iterator cell = this->triangulation.begin_active(), 
                                                      endc = this->triangulation.end();

    unsigned int counter = 0;
    for (; cell != endc ; cell++)
    {
        VelocitySpace_error_per_cell(counter) = pow(VelocitySpace_error_per_cell(counter),2)/cell->measure();
        counter ++;
    }
        
    
}