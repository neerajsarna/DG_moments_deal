
template<int dim> void Solver_DG<dim>::solve()
         {
          SolverControl           solver_control (10000, 1e-10);
          TrilinosWrappers::SolverGMRES::AdditionalData additional_data;
          TrilinosWrappers::SolverGMRES  solver (solver_control,additional_data);


          TrilinosWrappers::PreconditionILU preconditioner;
          TrilinosWrappers::PreconditionILU::AdditionalData additional_data_PC(0,1e-5,1.01,0);
          preconditioner.initialize(global_matrix,additional_data_PC);
        

          /*PreconditionBlockSSOR<SparseMatrix<double> > preconditioner;
          preconditioner.initialize(system_data.system_matrix, system_data.finite_element.dofs_per_cell);*/

          solver.solve (global_matrix, solution,system_rhs,preconditioner);

        }
