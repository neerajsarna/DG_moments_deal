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


#include "BCrhs_systemA.h"

#include "ExactSolution_Dummy.h"
#include "ExactSolution_PoissonHeatConduction_G26.h"
#include "ExactSolution_SystemA.h"
#include "ExactSolution_PoissonHeatConduction_G20.h"
#include "ExactSolution_PoissonHeatConduction_G35.h"
#include "ExactSolution_PoissonHeatConduction_G45.h"
#include "ExactSolution_PoissonHeatConduction_G56.h"
#include "ExactSolution_PoissonHeatConduction_G71.h"
#include "ExactSolution_PoissonHeatConduction_G84.h"
#include "ExactSolution_PoissonHeatConduction_G105.h"
#include "ExactSolution_PoissonHeatConduction_G120.h"

#include "SystemA.h"
#include "G26.h"
#include "G20.h"
#include "G35.h"
#include "G45.h"
#include "G56.h"
#include "G71.h"
#include "G84.h"
#include "G105.h"
#include "G120.h"
#include "Run.h"
#include "MeshGenerator_Square_CircularCavity.h"
#include "MeshGenerator_Square_CircularCavity_Channel.h"
#include "MeshGenerator_NACA.h"
#include "MeshGenerator_cylinder.h"


