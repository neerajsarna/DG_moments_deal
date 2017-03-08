#include "Base_Constants.h"
#include "Base_TensorInfo.h"
#include "Base_ForceType.h"
#include "Base_MatrixOpt.h"
#include "Base_BoundaryHandler.h"
#include "Base_BCrhs.h"
#include "Base_ExactSolution.h"
#include "Base_EquationGenerator.h"
#include "LinearSolver.h"
#include "Base_MeshGenerator.h"
#include "Base_PeriodicityHandler.h"
#include "Base_PostProc.h"
#include "Base_NumericalIntegration.h"
#include "Base_Solver.h"
#include "Develop_System.h"


#include "BCrhs_systemA.h"

#include "ExactSolution_Dummy.h"
#include "ExactSolution_PoissonHeat.h"
// #include "ExactSolution_PoissonHeatConduction_G26.h"
// #include "ExactSolution_SystemA.h"
// #include "ExactSolution_PoissonHeatConduction_G20.h"
// #include "ExactSolution_PoissonHeatConduction_G35.h"
// #include "ExactSolution_PoissonHeatConduction_G45.h"
// #include "ExactSolution_PoissonHeatConduction_G56.h"
// #include "ExactSolution_PoissonHeatConduction_G71.h"
// #include "ExactSolution_PoissonHeatConduction_G84.h"
// #include "ExactSolution_PoissonHeatConduction_G105.h"
// #include "ExactSolution_PoissonHeatConduction_G120.h"

#include "SystemA.h"
#include "Run.h"



