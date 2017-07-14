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

#include "Base_PostProc.h"
#include "Base_NumericalIntegration.h"
#include "BCrhs_systemA.h"
#include "Develop_System.h"

#include "ExactSolution_Dummy.h"
#include "ExactSolution_PoissonHeat.h"
#include "ExactSolution_SystemA.h"
#include "SystemA.h"
#include "FE_Data.h"

#include "Base_PeriodicityHandler.h"
#include "Assembler.h"
#include "Assembler_Periodic.h"
#include "Run_Problem.h"



