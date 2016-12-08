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


    TEST(ResidualComputation,HandlesResidualComputation)
    {
    	MKL_INT *ia;
    	MKL_INT *ja;
    	double *values;

    	// consider a diagonal matrix first
    	const unsigned int n_rows = 3;
    	const unsigned int nnz = 3;

    	ia =(MKL_INT*)calloc(n_rows + 1,sizeof(MKL_INT));
    	ja = (MKL_INT*)calloc(nnz,sizeof(MKL_INT));
    	values = (double*)calloc(nnz,sizeof(double));

    	ia[0] = 0;
    	ia[1] = 1;
    	ia[2] = 2;
    	ia[3] = nnz;

    	ja[0] = 0;
    	ja[1] = 1;
    	ja[2] = 2;

    	values[0] = 1;
    	values[1] = 1;
    	values[2] = 1;

    	Vector<double> vec(3);
    	vec(0) = 1;
    	vec(1) = 2;
    	vec(2) = 3;

    	MatrixOpt::Base_MatrixOpt matrix_opt;
    	Vector<double> result = matrix_opt.Sparse_matrix_dot_Vector(ia, ja,
                                            values,vec,n_rows);



    	std::cout << "Results Sparse Matrix Product " << result << std::endl;


    }
}
