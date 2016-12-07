namespace Test_LinearSolver
{
	using namespace dealii;
	TEST(LinearSolver,HandlesLinearSolver)
	{
		TrilinosWrappers::SparseMatrix global_matrix;
	    TrilinosWrappers::SparsityPattern sparsity_pattern;
		DynamicSparsityPattern dsp(5,5);
		Vector<double> system_rhs(5);
		Vector<double> solution(5);

		for(unsigned int i = 0 ; i < 5 ; i ++)
			dsp.add(i,i);

		sparsity_pattern.copy_from(dsp);
		sparsity_pattern.compress();
		global_matrix.reinit(sparsity_pattern);
	
		for(unsigned int i = 0 ; i < 5 ; i ++)
		{
			global_matrix.add(i,i,1);
			system_rhs(i) = 1;
		}


		LinearSolver::LinearSolver linear_solver;
		linear_solver.develop_pardiso_data(global_matrix,sparsity_pattern);
		double res = linear_solver.solve_with_pardiso(system_rhs,solution);

 		std::cout << "Solution " << std::endl;
		std::cout << solution << std::endl;             
		std::cout << "Residual " << res << std::endl;
    }


}
