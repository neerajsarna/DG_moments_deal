template<int num_flux,int dim> 
void 
Solver_DG<num_flux,dim>::error_evaluation(const Vector<double> solution)
      {

        unsigned int component = 0;                                                // the component for which the error has to be computed
        Vector<double> error_per_cell(triangulation.n_active_cells());      

        ComponentSelectFunction<dim> weight(component,this->nEqn);                              // used to compute only the error in theta

        // computation of L2 error
        VectorTools::integrate_difference (mapping,dof_handler,solution,
          *exact_solution,
          error_per_cell,
          QGauss<dim>(ngp),
          VectorTools::L2_norm,&weight);  

        switch(equation_system_data->system_type)
        {
          case symmetric:
          {
            Assert(nEqn != 10,ExcNotImplemented());
            error_per_cell /= sqrt(2);
            break;
          }
          case un_symmetric:
          {
            break;
          }
        }

        const double L2_error = error_per_cell.l2_norm();

        // computation of L_inifinity error
        error_per_cell = 0;
        VectorTools::integrate_difference (mapping,dof_handler,solution,
          *exact_solution,
          error_per_cell,
          QGauss<dim>(ngp),
          VectorTools::Linfty_norm,&weight);  
        
        switch(equation_system_data->system_type)
        {
          case symmetric:
          {
            Assert(nEqn != 10,ExcNotImplemented());
            error_per_cell /= sqrt(2);
            break;
          }
          case un_symmetric:
          {
            break;
          }
        }

        const double Linfty_error = error_per_cell.linfty_norm();
                
        string column_name_L2;
        string column_name_Linfty;

        column_name_L2 = "#L2 in u" + to_string(component);
        column_name_Linfty = "#Linfty in u" + to_string(component);

        string error_details;
        error_details = " L2_error: " + to_string(L2_error) +" Linfty_error: " + to_string(Linfty_error) + " #DOF: " + to_string(dof_handler.n_dofs()) +
                        " #Cells " + to_string(triangulation.n_active_cells()); 

        cout << "*******************Error details******************"<< endl ;
        cout << error_details << endl;
        cout << "**************************************************" << endl ;

        convergence_table.add_value(column_name_L2,L2_error);
        convergence_table.add_value(column_name_Linfty,Linfty_error);
        convergence_table.add_value("#degree of freedom",dof_handler.n_dofs());
        convergence_table.add_value("#number of cells",triangulation.n_active_cells());

        convergence_table.set_scientific(column_name_L2,true);
        convergence_table.set_scientific(column_name_Linfty,true);

      }
