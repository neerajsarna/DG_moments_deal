#include "include_deal.h"
#include "gtest/gtest.h"
#include "EigenSetup.h"
#include "basic_data_structures.h"
#include "basic_routines.h"
#include "mkl.h"

// core routines
#include "Base_Constants.h"
#include "Base_TensorInfo.h"
 #include "Base_ForceType.h"
 #include "Base_MatrixOpt.h"
 #include "Base_BoundaryHandler.h"
 #include "BCrhs_systemA.h"
 #include "Base_ExactSolution.h"
 #include "ExactSolution_SystemA.h"
 #include "Base_EquationGenerator.h"
 #include "LinearSolver.h"
 #include "SystemA.h"
 #include "Base_MeshGenerator.h"
 #include "Base_PeriodicityHandler.h"
 #include "Base_PostProc.h"
 #include "Base_Solver.h"

//tests
//#include "Test_BaseConstants.h"
 #include "Test_ForceType.h" 
 #include "Test_EquationGenerator.h"
 #include "Test_MatrixOpt.h"
#include "Test_BoundaryHandler_Char.h"
#include "Test_TensorInfo.h"
#include "Test_BoundaryHandler_Odd.h"
#include "Test_ExactSolution.h"
#include "Test_Solver.h"

int main(int argc,char **argv)
{
	dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc,argv, 1);

	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}