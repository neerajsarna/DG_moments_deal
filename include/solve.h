
template<int force_type,int system_type,int num_flux,int dim> 
void 
Solver_DG<force_type,system_type,num_flux,dim>
::solve(const enum Solver_Type solver_type)
         {

          switch(solver_type)
          {
            case Trilinos_Direct:
            {
                SolverControl           solver_control (10000, 1e-10);
          
                TrilinosWrappers::SolverDirect::AdditionalData additional_data;
                TrilinosWrappers::SolverDirect solver(solver_control,additional_data);
                solver.solve (global_matrix, solution,system_rhs);

                break;
            }

            case Trilinos_GMRES:
            {
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
              system_matrix temp_global_matrix;

              temp_global_matrix.Row_Col_Value.reserve(global_matrix.n_nonzero_elements());
              temp_global_matrix.matrix.resize(global_matrix.m(),global_matrix.n());

              for (unsigned int i = 0 ; i < global_matrix.m(); i++)
                temp_global_matrix.matrix.coeffRef(i,i) = 1;
             /* for (unsigned int row = 0 ; row < global_matrix.m(); row++)
              {
                typename TrilinosWrappers::SparseMatrix::const_iterator it = global_matrix.begin(row),
                                                                        it_end = global_matrix.end(row);

                for (; it != it_end ; it++)
                  temp_global_matrix.Row_Col_Value.push_back(triplet(it->row(),it->column(),it->value()));


              }

              global_matrix.clear();
              temp_global_matrix.matrix.setFromTriplets(temp_global_matrix.Row_Col_Value.begin(),
                                                        temp_global_matrix.Row_Col_Value.end());
              temp_global_matrix.Row_Col_Value.clear();*/

              MKL_INT iparm[64];
              MKL_INT mtype = 11;       // unsymmetric real matrix
              MKL_INT phase;

              void *pt[64];

              for(int i=0; i<64; i++ ) pt[i]=0;

              phase = -11;;               // setup phase
              PardisoSolve(mtype, temp_global_matrix.matrix.rows(), (MKL_INT *)temp_global_matrix.matrix.outerIndexPtr(),
                           (MKL_INT *)temp_global_matrix.matrix.innerIndexPtr(), temp_global_matrix.matrix.valuePtr(), 
                          temp_global_matrix.matrix.nonZeros(), &system_rhs(0), &solution(0), pt, phase, iparm);

              phase = 11;                 // symbolic factorization
              PardisoSolve(mtype, temp_global_matrix.matrix.rows(), (MKL_INT *)temp_global_matrix.matrix.outerIndexPtr(),
                           (MKL_INT *)temp_global_matrix.matrix.innerIndexPtr(), temp_global_matrix.matrix.valuePtr(), 
                          temp_global_matrix.matrix.nonZeros(), &system_rhs(0), &solution(0), pt, phase, iparm);

              phase = 22;               // numerical factorization
              PardisoSolve(mtype, temp_global_matrix.matrix.rows(), (MKL_INT *)temp_global_matrix.matrix.outerIndexPtr(),
                           (MKL_INT *)temp_global_matrix.matrix.innerIndexPtr(), temp_global_matrix.matrix.valuePtr(), 
                          temp_global_matrix.matrix.nonZeros(), &system_rhs(0), &solution(0), pt, phase, iparm);
              
              phase = 33;             //Solve the system
              PardisoSolve(mtype, temp_global_matrix.matrix.rows(), (MKL_INT *)temp_global_matrix.matrix.outerIndexPtr(),
                           (MKL_INT *)temp_global_matrix.matrix.innerIndexPtr(), temp_global_matrix.matrix.valuePtr(), 
                          temp_global_matrix.matrix.nonZeros(), &system_rhs(0), &solution(0), pt, phase, iparm);
              
              break;
            }
          }



        }


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
  
/*  MKL_INT *ia, *ja;
  ia = (MKL_INT *)malloc(sizeof(MKL_INT)*nz);
  ja = (MKL_INT *)malloc(sizeof(MKL_INT)*nz);
  for( int i=0; i<nz; i++ ) {
    ia[i] = (MKL_INT)ia0[i];
    ja[i] = (MKL_INT)ja0[i];
  };*/
  
  printf( " PardisoSolve :" );
  switch( phase ) 
  {
    case( -11 ):
      printf( " Setup Phase \n" );
      /*
      printf("\nINTERNAL: The matrix is given by: ");fflush(stdout);
      for (i = 0; i < nz; i+=1) {
        printf("\n i=%d, j=%d, A(i,j)=%f", ia[i], ja[i], ((double*)A)[i] );fflush(stdout);
      };  printf ("\n");
      */
      
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

}
