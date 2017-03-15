// following is simply a solver of linear systems
// It is independent of the assemblation and the system properties

namespace LinearSolver
{
  using namespace dealii;
  
	class
	LinearSolver
	{

	public:
    enum Solver_Type
    {Trilinos_Direct,Trilinos_GMRES};
    const Solver_Type solver_type = Trilinos_Direct;

		LinearSolver(){};


    void solve_trilinos(TrilinosWrappers::SparseMatrix &global_matrix,Vector<double> &system_rhs,
                 Vector<double> &solution);

    void solve_eigen(Sparse_matrix &global_matrix,Vector<double> &system_rhs,
                    Vector<double> &solution);

    Vector<double> solve_with_pardiso(Vector<double> &system_rhs,
                          Vector<double> &solution);
		void PardisoSolve(MKL_INT *ia,MKL_INT *ja,double *a,double *b,double *x,MKL_INT n);

    MatrixOpt::Base_MatrixOpt matrix_opt;

    // data for pardiso
    MKL_INT *IA;
    MKL_INT *JA;
    double *V;
    unsigned int n_rows;
    unsigned int nnz;

    void develop_pardiso_data(TrilinosWrappers::SparseMatrix &matrix);

	};

  void LinearSolver::develop_pardiso_data(TrilinosWrappers::SparseMatrix &matrix)
  {
   n_rows = matrix.m();
   nnz = matrix.n_nonzero_elements();

   Epetra_CrsMatrix mat1 = matrix.trilinos_matrix();

   matrix.clear();
  

   //Epetra_CrsMatrix mat2 = mat1;
   Epetra_IntSerialDenseVector row_ptr = mat1.ExpertExtractIndexOffset();
   Epetra_IntSerialDenseVector col_ind = mat1.ExpertExtractIndices();
   double *values = mat1.ExpertExtractValues();	 

   IA = (MKL_INT*)calloc(n_rows+1,sizeof(MKL_INT));
   JA = (MKL_INT*)calloc(nnz,sizeof(MKL_INT));
   V = (double*)calloc(nnz,sizeof(double));   

   for (unsigned int i = 0 ; i < n_rows + 1; i++)
		IA[i] = row_ptr[i];

   for (unsigned int i = 0 ; i < nnz ; i++)
   {
	JA[i] = col_ind[i];
	V[i] = values[i];
   }
	
 //  matrix_opt.COO_to_CSR(matrix,IA,JA,V);

    // now we release the old memory
    matrix.clear();

  }


  Vector<double> LinearSolver::solve_with_pardiso(Vector<double> &system_rhs,
                                        Vector<double> &solution)
  {
    Assert(system_rhs.size() == n_rows,ExcNotInitialized());
    Assert(solution.size() == n_rows,ExcNotInitialized());

      PardisoSolve(IA,JA,                    
                  V,&system_rhs(0),&solution(0),n_rows);

      Vector<double> res = matrix_opt.Sparse_matrix_dot_Vector(IA,JA,
                                          V,solution,n_rows);

      res -=system_rhs;

      free(IA);
      free(JA);
      free(V);

      return(res);
  }

  // solving routines for a trilinos sparse matrix
	void LinearSolver::solve_trilinos(TrilinosWrappers::SparseMatrix &global_matrix,Vector<double> &system_rhs,
		                        Vector<double> &solution)
	{

    Assert(1 == 0, ExcMessage("Please use pardiso"));
          
          switch(solver_type)
          {
            case Trilinos_Direct:
            {
				        std::cout << "using Direct" << std::endl;
                SolverControl           solver_control (10000, 1e-10);
          
                TrilinosWrappers::SolverDirect::AdditionalData additional_data;
                TrilinosWrappers::SolverDirect solver(solver_control,additional_data);
                solver.solve (global_matrix, solution,system_rhs);

                break;
            }

            case Trilinos_GMRES:
            {
	      	  std::cout << "using GMRES" << std::endl;
              SolverControl           solver_control (10000, 1e-10);
              TrilinosWrappers::SolverGMRES::AdditionalData additional_data;
              TrilinosWrappers::SolverGMRES  solver (solver_control,additional_data);

	      std::cout << "Preparing Preconditionar" << std::endl;
              TrilinosWrappers::PreconditionILU preconditioner;
              TrilinosWrappers::PreconditionILU::AdditionalData additional_data_PC(0,1e-5,1.01,0);
              preconditioner.initialize(global_matrix,additional_data_PC);
        
              MatrixOpt::Base_MatrixOpt matrix_opt;
              matrix_opt.print_dealii_sparse(global_matrix,"global_matrix_inside_lin_solve");
              matrix_opt.print_dealii_vector(system_rhs,"system_rhs_inside_lin_solve");
             
	      fflush(stdout);
 	      std::cout << "solving...." << std::endl;
	      solver.solve (global_matrix, solution,system_rhs,preconditioner);
              break;
            }

            default:
            {
              Assert(1 == 0, ExcMessage("Should not have reached here"));
              break;

            }
          }

	}

  void LinearSolver::solve_eigen(Sparse_matrix &global_matrix,Vector<double> &system_rhs,
                            Vector<double> &solution)
  {
    std::cout << "Solving Eigen Matrix " << std::endl;
    //Assert(solver_type == Pardiso,ExcMessage("No other solver works with eigen matrices presently"));
    Assert(global_matrix.isCompressed(),ExcMessage("Compressed Matrix is a must"));

    const int *row_ptr = global_matrix.outerIndexPtr();
    const int *col_ind = global_matrix.innerIndexPtr();

    const int n_rows = global_matrix.rows();
    const int nnz = global_matrix.nonZeros();

    MKL_INT *ia;
    MKL_INT *ja;


    ia = (MKL_INT*)calloc(n_rows+1,sizeof(MKL_INT));
    ja = (MKL_INT*)calloc(nnz,sizeof(MKL_INT));


    for (long int i = 0 ; i < n_rows + 1 ; i++)
      ia[i]  = row_ptr[i];


    for (long int j = 0 ; j < nnz ; j ++)
      ja[j] = col_ind[j];

    PardisoSolve(ia , ja , global_matrix.valuePtr(),
                 &system_rhs(0),&solution(0),global_matrix.rows());
  }


	

	void 
	LinearSolver
	::PardisoSolve(MKL_INT *ia,MKL_INT *ja,double *a,double *b,double *x,MKL_INT n)
	{
		MKL_INT mtype = 11;       

		MKL_INT nrhs = 1;     
		void *pt[64];


		MKL_INT iparm[64];
		MKL_INT maxfct, mnum, phase, error, msglvl;

		MKL_INT i;
		double ddum;          
		MKL_INT idum;         
		std::string uplo;

		for ( i = 0; i < 64; i++ )
			iparm[i] = 0;

      iparm[0] = 1; // No solver default  
      iparm[1] = 3; // Fill-in reordering from METIS  
      iparm[2] = 3; // Numbers of processors, value of OMP_NUM_THREADS  
      iparm[3] = 0; // No iterative-direct algorithm  
      iparm[4] = 0; // No user fill-in reducing permutation  
      iparm[5] = 0; // 0 = Write solution into x  
      iparm[6] = 0; // Not in use  
      iparm[7] = 200; // Max numbers of iterative refinement steps  
      iparm[8] = 0; // Not in use  
      iparm[9] = 13; // Perturb the pivot elements with 1E-13  
      iparm[10] = 1; // 1 = Use nonsymmetric permutation and scaling MPS  
      iparm[11] = 0; // Not in use  
      iparm[12] = 0; // 1 = Maximum weighted matching algorithm is switched-on (default for non-symmetric)  
      iparm[13] = 0; // Output: Number of perturbed pivots  
      iparm[14] = 0; // Not in use  
      iparm[15] = 0; // Not in use  
      iparm[16] = 0; // Not in use  
      iparm[17] = -1; // Output: Number of nonzeros in the factor LU  
      iparm[18] = -1; // Output: Mflops for LU factorization  
      iparm[19] = 0; // Output: Numbers of CG Iterations  
      iparm[34] = 1; // 0= one-based indices, 1= zero-besed indices  

    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;         /* Which factorization to use. */
    msglvl = 0;           /* Print statistical information in file */
    error = 0;            /* Initialize error flag */

      for ( i = 0; i < 64; i++ )
      	pt[i] = 0;

    // Symbolic Factorization
      phase = 11;
      std::cout << "Symbolic Factorization" << std::endl;
      PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
      	&n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
      if ( error != 0 )
      {
      	printf ("\nERROR during symbolic factorization: %lld", error);
      	exit (1);
      }
      printf ("\nReordering completed ... ");

    // Numerical Factorization
      phase = 22;
      std::cout << "Numerical Factorization" << std::endl;
      PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
      	&n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
      if ( error != 0 )
      {
      	printf ("\nERROR during numerical factorization: %lld", error);
      	exit (2);
      }
      printf ("\nFactorization completed ... \n");

    // Solving the system
      phase = 33;

    iparm[11] = 0;        /* Conjugate transposed/transpose solve */
      uplo = "non-transposed";


     std::cout << "Solving the system " << std::endl;

      pardiso (pt, &maxfct, &mnum, &mtype, &phase,
      	&n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);

      if ( error != 0 )
      {
      	printf ("\nERROR during solution: %lld", error);
      	exit (3);
      }

      


    phase = -1;           /* Release internal memory. */
      PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
      	&n, &ddum, ia, ja, &idum, &nrhs,
      	iparm, &msglvl, &ddum, &ddum, &error);

  }


}
