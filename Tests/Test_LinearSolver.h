namespace Test_LinearSolver
{
	using namespace dealii;
	TEST(LinearSolver,HandlesLinearSolver)
	{
		TrilinosWrappers::SparseMatrix global_matrix;
	        SparsityPattern sparsity_pattern;
		DynamicSparsityPattern dsp(5,5);
		Vector<double> system_rhs(5);
		Vector<double> solution(5);

		for(unsigned int i = 0 ; i < 5 ; i ++)
			dsp.add(i,i);

		sparsity_pattern.copy_from(dsp);
		global_matrix.reinit(sparsity_pattern);
	
		for(unsigned int i = 0 ; i < 5 ; i ++)
		{
			global_matrix.add(i,i,1);
			system_rhs(i) = 1;
		}


              
 	      std::cout << "using GMRES" << std::endl;
              SolverControl           solver_control (10000, 1e-10);
              TrilinosWrappers::SolverGMRES::AdditionalData additional_data;
              TrilinosWrappers::SolverGMRES  solver (solver_control,additional_data);

              std::cout << "Preparing Preconditionar" << std::endl;
              TrilinosWrappers::PreconditionILU preconditioner;
              TrilinosWrappers::PreconditionILU::AdditionalData additional_data_PC(0,1e-5,1.01,0);
              preconditioner.initialize(global_matrix,additional_data_PC);

              std::cout << "solving...." << std::endl;
              solver.solve (global_matrix, solution,system_rhs,preconditioner);

           }


}
