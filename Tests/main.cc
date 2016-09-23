#include "include_deal.h"
#include "gtest/gtest.h"
#include "EigenSetup.h"
#include "basic_data_structures.h"

// core routines
#include "Base_Constants.h"
 #include "Base_ForceType.h"
 #include "Base_MatrixOpt.h"
 #include "Base_EquationGenerator.h"
#include "Base_BoundaryHandler.h"
 #include "SystemA.h"


//tests
#include "Test_BaseConstants.h"
 #include "Test_ForceType.h" 
 #include "Test_EquationGenerator.h"
 #include "Test_MatrixOpt.h"

int main(int argc,char **argv)
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}