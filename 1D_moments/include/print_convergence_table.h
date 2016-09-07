template<int num_flux,int dim> 
void 
Solver_DG<num_flux,dim>::print_convergence_table(string filename)
        {
       	   unsigned int component = this->variable_map.find(error_variable)->second;  

       	   string column_name_L2;
       	   string column_name_Linfty;

       	   column_name_L2 = "#L2 in u" + to_string(component);
       	   column_name_Linfty = "#Linfty in u" + to_string(component);

       	
		   convergence_table.evaluate_convergence_rates(column_name_L2,"#number of cells", ConvergenceTable::reduction_rate_log2);
           convergence_table.evaluate_convergence_rates(column_name_Linfty,"#number of cells", ConvergenceTable::reduction_rate_log2);

           ofstream output_convergence(filename);
           convergence_table.write_text(output_convergence);

         }
