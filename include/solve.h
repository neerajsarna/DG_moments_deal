
template<int system_type,int num_flux,int dim> 
void 
Solver_DG<system_type,num_flux,dim>
::solve()
         {
          SolverControl           solver_control (10000, 1e-10);
          
          /*TrilinosWrappers::SolverDirect::AdditionalData additional_data;
          TrilinosWrappers::SolverDirect solver(solver_control,additional_data);
          solver.solve (global_matrix, solution,system_rhs);*/

          TrilinosWrappers::SolverGMRES::AdditionalData additional_data;
          TrilinosWrappers::SolverGMRES  solver (solver_control,additional_data);

          TrilinosWrappers::PreconditionILU preconditioner;
          TrilinosWrappers::PreconditionILU::AdditionalData additional_data_PC(0,1e-5,1.01,0);
          preconditioner.initialize(global_matrix,additional_data_PC);
        
          solver.solve (global_matrix, solution,system_rhs,preconditioner);
          

        }
