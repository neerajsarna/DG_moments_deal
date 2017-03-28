#include "include_deal.h"
#include "gtest/gtest.h"
#include "EigenSetup.h"
#include "basic_data_structures.h"
#include "basic_routines.h"
#include "parse_command_line.h"
#include "mkl.h"

std::vector<std::string> input_file;

// core routines
#include "Base_Constants.h"
#include "Base_TensorInfo.h"
#include "Base_ForceType.h"
#include "Base_MatrixOpt.h"

#include "Base_BCrhs.h"
#include "BCrhs_systemA.h"

#include "Base_BoundaryHandler.h"
#include "Base_EquationGenerator.h"
#include "Develop_System.h"


#include "Base_ExactSolution.h"
#include "ExactSolution_Dummy.h"
#include "ExactSolution_SystemA.h"

#include "LinearSolver.h"
#include "Base_MeshGenerator.h"
#include "Base_PeriodicityHandler.h"
#include "Base_PostProc.h"
#include "Base_NumericalIntegration.h"
#include "Base_Solver.h"

#include "ExactSolution_PoissonHeat.h"

// certain generic tests which have to be performed independent of the system
//#include "Test_LinearSolver.h"
//#include "Test_FEValues.h"
//#include "Test_MatrixOpt.h"
//#include "Test_MeshGenerator.h"
//#include "Test_VertexDof.h"
//#include "Test_Solver_Generic.h"
//#include "Test_Solver_m_adaptive.h"
//#include "Test_Assembly.h"
//#include "Test_Assembly_MAdaptive.h"
#include "Test_Run.h"
//#include "Test_ResidualCompution.h"
//#include "Test_Projection.h"

#include "SystemA.h"

//tests

/**********Standard Tests for any system. Change The include directory from the make file*************/
//#include "Test_BaseConstants.h"
 //#include "Test_ForceType.h" 
 //#include "Test_EquationGenerator.h"
 //#include "Test_TensorInfo.h"
//#include "Test_BoundaryHandler_Odd.h"
 //#include "Test_ExactSolution.h"
//#include "Test_SquareGrid.h"
//  #include "Test_CrsFormat.h"
 //#include "Test_Production.h"
// #include "Test_ProjectorData.h"
// #include "Test_SymmetrizerData.h"

int main(int argc,char **argv)
{
	const unsigned int num_threads = atoi(argv[1]);
	input_file.resize(2);

	dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc,argv, num_threads);

	parse_command_line commmand_line_parser(argc,argv);
	/*store the filename read from the command line*/
	commmand_line_parser.read_input_filename(input_file);

	::testing::InitGoogleTest(&argc, argv);
	const unsigned int result = RUN_ALL_TESTS();

    std::cout <<"Results " << result<< std::endl;
	return 0;
}
