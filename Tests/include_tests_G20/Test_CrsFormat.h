namespace TestCRS
{
	using namespace dealii;

	// TEST(TestCrs,HandlesCrs)
	// {
	// 	const unsigned int dim = 2;
	// 	ASSERT_EQ(dim,2) << "3D not implemented" << std::endl;

	// 	std::string folder_name = "../system_matrices/";
	// 	Constants::Base_Constants constants(input_file);
	// 	G20::G20<dim> G20(constants.constants,folder_name);

	// 	ExactSolution::G20_PoissonHeat<dim>  G20_PoissonHeat(constants.constants,G20.base_tensorinfo.S_half);

	// 	FEM_Solver::Base_Solver<dim> base_solver("mesh_file_name",
	// 										 "grid",
	// 										 constants.constants,
	// 										 &G20,
	// 										 &G20_PoissonHeat);


	// 	base_solver.distribute_dof_allocate_matrix();
	// 	base_solver.assemble_system_meshworker();

	// 	typename TrilinosWrappers::SparseMatrix::const_iterator it = base_solver.global_matrix.begin(),
	// 													it_end = base_solver.global_matrix.end();


		
	// 	std::vector<size_t> row_index;

	// 	for (; it != it_end ; it++)
	// 		row_index.push_back(it->row());

	// 	for (unsigned int i = 0 ; i < row_index.size() -1 ; i++)
	// 		Assert(row_index[i] <= row_index[i + 1],ExcMessage("Matrix not in Row Major form"));
	// }

	TEST(TestCOOtoCRS,HandlesCOOtoCRS)
	{
		MatrixOpt::Base_MatrixOpt matrix_opt;
		TrilinosWrappers::SparseMatrix global_matrix;
		DynamicSparsityPattern dsp(3,3);
		SparsityPattern sparsity_pattern;

		dsp.add(0,0);
		dsp.add(0,1);
		dsp.add(1,1);
		dsp.add(2,2);

		sparsity_pattern.copy_from(dsp);
		global_matrix.reinit(sparsity_pattern);

		global_matrix.add(0,0,1);
		global_matrix.add(0,1,1);
		global_matrix.add(1,1,1);
		global_matrix.add(2,2,1);

		MKL_INT *IA;
		MKL_INT *JA;
		double *V;

		const unsigned int n_rows = global_matrix.m();
		const unsigned int nnz = global_matrix.n_nonzero_elements();


		IA = (MKL_INT*)calloc(n_rows+1,sizeof(MKL_INT));
        JA = (MKL_INT*)calloc(nnz,sizeof(MKL_INT));
        V = (double*)calloc(nnz,sizeof(double));

 
 	 	std::cout << "Non Zeros "<< global_matrix.n_nonzero_elements() << std::endl;
        matrix_opt.COO_to_CSR(global_matrix,IA,JA,V);

        for (unsigned int i = 0 ; i < n_rows + 1 ; i++)
        	std::cout << "IA " << IA[i] <<std::endl;

        for (unsigned int i = 0 ; i < nnz ; i++)
        	std::cout << "V " << V[i] << std::endl;

	}

	TEST(TestCOOtoCRSBigMat,HandlesCOOtoCRSBigMat)
	{
		const unsigned int dim = 2;
		ASSERT_EQ(dim,2) << "3D not implemented" << std::endl;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		G20::G20<dim> G20(constants.constants,folder_name);

		ExactSolution::G20_PoissonHeat<dim>  G20_PoissonHeat(constants.constants,G20.base_tensorinfo.S_half);

		FEM_Solver::Base_Solver<dim> base_solver("mesh_file_name",
											 "grid",
											 constants.constants,
											 &G20,
											 &G20_PoissonHeat);


		base_solver.distribute_dof_allocate_matrix();
		base_solver.assemble_system_meshworker();

		const unsigned int n_rows = base_solver.global_matrix.m();
		const unsigned int nnz = base_solver.global_matrix.n_nonzero_elements();

		MKL_INT *IA;
		MKL_INT *JA;
		double *V;


		IA = (MKL_INT*)calloc(n_rows+1,sizeof(MKL_INT));
        JA = (MKL_INT*)calloc(nnz,sizeof(MKL_INT));
        V = (double*)calloc(nnz,sizeof(double));

        MatrixOpt::Base_MatrixOpt matrix_opt;
        matrix_opt.COO_to_CSR(base_solver.global_matrix,IA,JA,V);

		const Epetra_CrsMatrix temp_system_matrix = base_solver.global_matrix.trilinos_matrix();
        Epetra_CrsMatrix temp_system_matrix2 = temp_system_matrix; 
        Epetra_IntSerialDenseVector row_ptr = temp_system_matrix2.ExpertExtractIndexOffset();
        Epetra_IntSerialDenseVector col_ind = temp_system_matrix2.ExpertExtractIndices();
        double *values = temp_system_matrix2.ExpertExtractValues();

    	for (unsigned int i = 0 ; i < n_rows + 1 ;i++)
    		EXPECT_EQ(IA[i],row_ptr[i]);

    	for (unsigned int i = 0 ; i < nnz ; i++)
    	{
    		EXPECT_EQ(JA[i],col_ind[i]);
    		EXPECT_NEAR(V[i],values[i],1e-5);
    	}


	}


	TEST(TestPardiso,HandlesPardiso)
	{
		const unsigned int dim = 2;
		ASSERT_EQ(dim,2) << "3D not implemented" << std::endl;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		G20::G20<dim> G20(constants.constants,folder_name);

		ExactSolution::G20_PoissonHeat<dim>  G20_PoissonHeat(constants.constants,G20.base_tensorinfo.S_half);

		FEM_Solver::Base_Solver<dim> base_solver("mesh_file_name",
											 "grid",
											 constants.constants,
											 &G20,
											 &G20_PoissonHeat);


		std::cout << "Dof Distribution " << std::endl;
		base_solver.distribute_dof_allocate_matrix();

		std::cout << "Assembling matrix" << std::endl;
		base_solver.assemble_system_meshworker();

		std::cout << "Linear Solver " << std::endl;
		LinearSolver::LinearSolver linear_solver;
		linear_solver.develop_pardiso_data(base_solver.global_matrix,base_solver.sparsity_pattern);
		double res = linear_solver.solve_with_pardiso(base_solver.system_rhs,base_solver.solution);

		std::cout << "Residual " << res << std::endl;
	}
}