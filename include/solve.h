
template<int force_type,int system_type,int num_flux,int dim> 
void 
Solver_DG<force_type,system_type,num_flux,dim>
::solve(const enum Solver_Type solver_type)
         {

          switch(solver_type)
          {
            case Trilinos_Direct:
            {
		cout << "using Direct" << endl;
                SolverControl           solver_control (10000, 1e-10);
          
                TrilinosWrappers::SolverDirect::AdditionalData additional_data;
                TrilinosWrappers::SolverDirect solver(solver_control,additional_data);
                solver.solve (global_matrix, solution,system_rhs);

                break;
            }

            case Trilinos_GMRES:
            {
	      cout << "using GMRES" << endl;
              SolverControl           solver_control (10000, 1e-10);
              TrilinosWrappers::SolverGMRES::AdditionalData additional_data;
              TrilinosWrappers::SolverGMRES  solver (solver_control,additional_data);

              TrilinosWrappers::PreconditionILU preconditioner;
              TrilinosWrappers::PreconditionILU::AdditionalData additional_data_PC(0,1e-5,1.01,0);
              preconditioner.initialize(global_matrix,additional_data_PC);
        
              solver.solve (global_matrix, solution,system_rhs,preconditioner);
              break;
            }

            case Pardiso:
            {
	      cout << "Using Pardiso" << endl;
              system_matrix temp_global_matrix;

              const unsigned int size_temp_matrix = global_matrix.m();
              const unsigned int non_zero_in_system_matrix = global_matrix.n_nonzero_elements();

              temp_global_matrix.Row_Col_Value.reserve(non_zero_in_system_matrix);
              temp_global_matrix.matrix.resize(size_temp_matrix,size_temp_matrix);

	      typename TrilinosWrappers::SparseMatrix::const_iterator it = global_matrix.begin(), it_end = global_matrix.end();
	   
	      cout << "trial start " << endl;
   	      for(; it != it_end ; it++)
	      {}
	      cout << "trial end " << endl;

	/*      cout<< "copying system matrix..." << endl;
              for (long long  int row = 0 ; row < global_matrix.m(); row++)
              {
                 it = global_matrix.begin(row);
                                                                        it_end = global_matrix.end(row);

                for (; it != it_end ; it++)
		{}
        //          temp_global_matrix.Row_Col_Value.push_back(triplet(it->row(),it->column(),it->value()));


              }*/

	     cout << "done copying " << endl;
	     fflush(stdout);
              global_matrix.clear();
              temp_global_matrix.matrix.setFromTriplets(temp_global_matrix.Row_Col_Value.begin(),
                                                        temp_global_matrix.Row_Col_Value.end());
              temp_global_matrix.Row_Col_Value.clear();
	      cout << "done clearing data " << endl;

	      cout << "copying to MKL_INT" << endl;
              temp_global_matrix.matrix.makeCompressed();

              const unsigned int n_rows = temp_global_matrix.matrix.rows();
              const unsigned int nnz = temp_global_matrix.matrix.nonZeros();

              MKL_INT *ia;
              MKL_INT *ja;
              

              ia = (MKL_INT*)calloc(n_rows+1,sizeof(MKL_INT));
              ja = (MKL_INT*)calloc(nnz,sizeof(MKL_INT));


              for (long long int i = 0 ; i < n_rows + 1 ; i++)
                ia[i]  = temp_global_matrix.matrix.outerIndexPtr()[i] + 1;
            

              for (long long int j = 0 ; j < nnz ; j ++)
                ja[j] = temp_global_matrix.matrix.innerIndexPtr()[j] + 1;

	      cout << "Done copying to MKL_INT " << endl;
	      cout << "entering pardiso " << endl;
              PardisoSolve(ia,ja,
                            temp_global_matrix.matrix.valuePtr(),&system_rhs(0),&solution(0),n_rows);

              free(ia);
              free(ja);

              break;
            }
          }

          



        }


/*
template<int force_type,int system_type,int num_flux,int dim> 
void 
Solver_DG<force_type,system_type,num_flux,dim>
::PardisoSolve(MKL_INT mtype, MKL_INT n, MKL_INT *ia, MKL_INT *ja, void *A, int nz,
                     void *b, void *x, void *pt, MKL_INT phase, MKL_INT *iparm)
{
  {  
  MKL_INT i, nrhs;
  MKL_INT maxfct, mnum, msglvl, error;       
  MKL_INT idum;    // Integer dummy 
  double ddum;  // Double dummy 
  nrhs = 1;    // number of rhs
  maxfct = 1;  // Maximum number of numerical factorizations
  mnum = 1;    // Which factorization to use.
  msglvl = 0;  // 1 = Print statistical information in file 
  error = 0;   // Initialize error flag 
  

  printf( " PardisoSolve :" );
  switch( phase ) 
  {
    case( -11 ):
      printf( " Setup Phase \n" );

      
      // Pardiso control parameters.  
      for (i = 0; i < 64; i++) iparm[i] = 0;             

      iparm[0] = 1; // No solver default  
      iparm[1] = 3; // Fill-in reordering from METIS  
      iparm[2] = 1; // Numbers of processors, value of OMP_NUM_THREADS  
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

    break;
    case( 11 ):    // Symbolic Factorization
      printf( " Symbolic Factorization \n" );
      pardiso(pt, &maxfct, &mnum, &mtype, &phase, 
        &n, A, (MKL_INT *)ia, (MKL_INT *)ja, &idum, &nrhs,
        iparm, &msglvl, &ddum, &ddum, &error);
      if (error != 0) {
        printf("\nERROR during symbolic factorization: %d\n", error);
        exit(1);
      }
    break;
    case( 22 ):    // Numerical Factorization
      printf( " Numerical Factorization \n" );
      pardiso(pt, &maxfct, &mnum, &mtype, &phase, 
        &n, A, (MKL_INT *)ia, (MKL_INT *)ja, &idum, &nrhs,
        iparm, &msglvl, &ddum, &ddum, &error);
      if (error != 0) {
        printf("\nERROR during numerical factorization: %d\n", error);
        exit(2);
      }
    break;
    case( 33 ):   //Back substitution and iterative refinement
      printf( " Back Substitution and Iterative Refinement \n" );
      pardiso(pt, &maxfct, &mnum, &mtype, &phase, 
        &n, A, (MKL_INT *)ia, (MKL_INT *)ja, &idum, &nrhs,
        iparm, &msglvl, b, x, &error);
      if (error != 0) {
        printf("\nERROR during solution: %d\n", error);
        exit(3);
      }
    break;
    case( -1 ): // Release internal memory. 
      printf( " Release Internal Memory \n" );
      pardiso(pt, &maxfct, &mnum, &mtype, &phase, 
        &n, A, (MKL_INT *)ia, (MKL_INT *)ja, &idum, &nrhs,
        iparm, &msglvl, &ddum, &ddum, &error);
    break;
  };
};

}*/

template<int force_type,int system_type,int num_flux,int dim> 
void 
Solver_DG<force_type,system_type,num_flux,dim>
::PardisoSolve(MKL_INT *ia,MKL_INT *ja,double *a,double *b,double *x,MKL_INT n)
{
      MKL_INT mtype = 11;       

    double *bs;
    double  res, res0;
    MKL_INT nrhs = 1;     
    void *pt[64];

    bs = (double*)calloc(n,sizeof(double));    

    MKL_INT iparm[64];
    MKL_INT maxfct, mnum, phase, error, msglvl;
    
    MKL_INT i, j;
    double ddum;          
    MKL_INT idum;         
    char *uplo;

    for ( i = 0; i < 64; i++ )
        iparm[i] = 0;

    iparm[0] = 1;         /* No solver default */
    iparm[1] = 2;         /* Fill-in reordering from METIS */
    iparm[3] = 0;         /* No iterative-direct algorithm */
    iparm[4] = 0;         /* No user fill-in reducing permutation */
    iparm[5] = 0;         /* Write solution into x */
    iparm[6] = 0;         /* Not in use */
    iparm[7] = 2;         /* Max numbers of iterative refinement steps */
    iparm[8] = 0;         /* Not in use */
    iparm[9] = 15;        /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0;        /* Conjugate transposed/transpose solve */
    iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
    iparm[13] = 0;        /* Output: Number of perturbed pivots */
    iparm[14] = 0;        /* Not in use */
    iparm[15] = 0;        /* Not in use */
    iparm[16] = 0;        /* Not in use */
    iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1;       /* Output: Mflops for LU factorization */
    iparm[19] = 0;        /* Output: Numbers of CG Iterations */
  
    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;         /* Which factorization to use. */
    msglvl = 0;           /* Print statistical information in file */
    error = 0;            /* Initialize error flag */

    for ( i = 0; i < 64; i++ )
        pt[i] = 0;

    // Symbolic Factorization
    phase = 11;
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
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during numerical factorization: %lld", error);
        exit (2);
    }
    printf ("\nFactorization completed ... ");

    // Solving the system
    phase = 33;

    iparm[11] = 0;        /* Conjugate transposed/transpose solve */
    uplo = "non-transposed";
    
    printf ("\n\nSolving %s system...\n", uplo);
    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
       &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);

    if ( error != 0 )
    {
        printf ("\nERROR during solution: %lld", error);
        exit (3);
    }

    printf ("\n");
// Compute residual
    mkl_dcsrgemv (uplo, &n, a, ia, ja, x, bs);
    res = 0.0;
    res0 = 0.0;
    for ( j = 1; j <= n; j++ )
    {
        res += (bs[j - 1] - b[j - 1]) * (bs[j - 1] - b[j - 1]);
        res0 += b[j - 1] * b[j - 1];
    }
    res = sqrt (res) / sqrt (res0);
    printf ("\nRelative residual = %e\n", res);
// Check residual
    if ( res > 1e-10 )
    {
        printf ("Error: residual is too high!\n");
        exit (10);
    }
    
    phase = -1;           /* Release internal memory. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error);

   free(bs);
}
