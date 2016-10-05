#include "include_deal.h"
#include "gtest/gtest.h"
#include "EigenSetup.h"
#include "basic_data_structures.h"
#include "basic_routines.h"
#include "parse_command_line.h"
#include "mkl.h"

std::string input_file;

// core routines
#include "Base_Constants.h"
#include "Base_TensorInfo.h"
 #include "Base_ForceType.h"
 #include "Base_MatrixOpt.h"
 #include "Base_BoundaryHandler.h"
#include "Base_BCrhs.h"
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

/**********TESTS FOR SYTEM A *************/
#include "Test_BaseConstants.h"
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
	const unsigned int num_threads = atoi(argv[1]);
	dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc,argv, num_threads);

	parse_command_line commmand_line_parser(argc,argv);
	/*store the filename read from the command line*/
	commmand_line_parser.read_input_filename(input_file);

	::testing::InitGoogleTest(&argc, argv);
	const unsigned int result = RUN_ALL_TESTS();

	return 0;
}