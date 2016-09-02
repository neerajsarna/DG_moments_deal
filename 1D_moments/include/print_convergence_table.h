template<int num_flux,int dim> 
void 
Solver_DG<num_flux,dim>::print_convergence_table(string filename)
        {
		   convergence_table.evaluate_convergence_rates("#L2 in u0","#number of cells", ConvergenceTable::reduction_rate_log2);
           convergence_table.evaluate_convergence_rates("#Linfty in u0","#number of cells", ConvergenceTable::reduction_rate_log2);

           ofstream output_convergence(filename);
           convergence_table.write_text(output_convergence);

         }
