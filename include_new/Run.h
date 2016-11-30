using namespace dealii;
void Run()
{
	const unsigned int dim = 2;
	Constants::Base_Constants constants(input_file);
	std::string folder_name = "../system_matrices/";

	switch(constants.constants.nEqn)
	{
		case 6:
		{
			// first develop the class which loads the systems
			SystemA::SystemA<dim> systemA(constants.constants,folder_name);			

			// routine to develop the exact solution
			ExactSolution::ExactSolution_SystemA_ring<dim>  exact_solution_systemA(constants.constants,systemA.base_tensorinfo.S_half);

			// the solver which combines the above two classes
			FEM_Solver::Base_Solver<dim> base_solver("mesh_file_name",
													 "grid",
													 constants.constants,
													 &systemA,
													 &exact_solution_systemA);

			// now run the simulations 
			switch(systemA.constants.mesh_type)
			{
				case ring:
				{
					base_solver.run_ring();
					break;
				}

				default:
				{
					Assert(1 == 0,ExcMessage("Should not have reached here"));
					break;
				}
				
			}
			break;
		}

		case 17:
		{
			G26::G26<dim> G26(constants.constants,folder_name);

			ExactSolution::G26_PoissonHeat<dim>  G26_PoissonHeat(constants.constants,G26.base_tensorinfo.S_half);

			FEM_Solver::Base_Solver<dim> base_solver("mesh_file_name",
													 "grid",
													 constants.constants,
													 &G26,
													 &G26_PoissonHeat);

			switch(G26.constants.mesh_type)
			{
				case periodic_square:
				{
					base_solver.run_periodic();
					break;
				}
				default:
				{
					Assert(1 == 0,ExcMessage("Should not have reached here"));
					break;
				}
			}
			break;
		}

		case 22:
		{
			G35::G35<dim> G35(constants.constants,folder_name);

			ExactSolution::G35_PoissonHeat<dim>  G35_PoissonHeat(constants.constants,G35.base_tensorinfo.S_half);

			FEM_Solver::Base_Solver<dim> base_solver("mesh_file_name",
													 "grid",
													 constants.constants,
													 &G35,
													 &G35_PoissonHeat);

			switch(G35.constants.mesh_type)
			{
				case periodic_square:
				{
					base_solver.run_periodic();
					break;
				}
				default:
				{
					Assert(1 == 0,ExcMessage("Should not have reached here"));
					break;
				}
			}
			break;
		}

		case 28:
		{
			G45::G45<dim> G45(constants.constants,folder_name);

			ExactSolution::G45_PoissonHeat<dim>  G45_PoissonHeat(constants.constants,G45.base_tensorinfo.S_half);

			FEM_Solver::Base_Solver<dim> base_solver("mesh_file_name",
													 "grid",
													 constants.constants,
													 &G45,
													 &G45_PoissonHeat);

			switch(G45.constants.mesh_type)
			{
				case periodic_square:
				{
					base_solver.run_periodic();
					break;
				}
				default:
				{
					Assert(1 == 0,ExcMessage("Should not have reached here"));
					break;
				}
			}
			break;
		}

		case 34:
		{
			G56::G56<dim> G56(constants.constants,folder_name);

			ExactSolution::G56_PoissonHeat<dim>  G56_PoissonHeat(constants.constants,G56.base_tensorinfo.S_half);

			FEM_Solver::Base_Solver<dim> base_solver("mesh_file_name",
													 "grid",
													 constants.constants,
													 &G56,
													 &G56_PoissonHeat);

			switch(G56.constants.mesh_type)
			{
				case periodic_square:
				{
					base_solver.run_periodic();
					break;
				}
				default:
				{
					Assert(1 == 0,ExcMessage("Should not have reached here"));
					break;
				}
			}
			break;
		}

		case 43:
		{
			G71::G71<dim> G71(constants.constants,folder_name);

			ExactSolution::G71_PoissonHeat<dim>  G71_PoissonHeat(constants.constants,G71.base_tensorinfo.S_half);

			FEM_Solver::Base_Solver<dim> base_solver("mesh_file_name",
													 "grid",
													 constants.constants,
													 &G71,
													 &G71_PoissonHeat);

			switch(G71.constants.mesh_type)
			{
				case periodic_square:
				{
					base_solver.run_periodic();
					break;
				}
				default:
				{
					Assert(1 == 0,ExcMessage("Should not have reached here"));
					break;
				}
			}
			break;
		}

		case 50:
		{
			G84::G84<dim> G84(constants.constants,folder_name);

			ExactSolution::G84_PoissonHeat<dim>  G84_PoissonHeat(constants.constants,G84.base_tensorinfo.S_half);

			FEM_Solver::Base_Solver<dim> base_solver("mesh_file_name",
													 "grid",
													 constants.constants,
													 &G84,
													 &G84_PoissonHeat);

			switch(G84.constants.mesh_type)
			{
				case periodic_square:
				{
					base_solver.run_periodic();
					break;
				}
				default:
				{
					Assert(1 == 0,ExcMessage("Should not have reached here"));
					break;
				}
			}
			break;
		}

		case 62:
		{
			G105::G105<dim> G105(constants.constants,folder_name);

			ExactSolution::G105_PoissonHeat<dim>  G105_PoissonHeat(constants.constants,G105.base_tensorinfo.S_half);

			FEM_Solver::Base_Solver<dim> base_solver("mesh_file_name",
													 "grid",
													 constants.constants,
													 &G105,
													 &G105_PoissonHeat);

			switch(G105.constants.mesh_type)
			{
				case periodic_square:
				{
					base_solver.run_periodic();
					break;
				}
				default:
				{
					Assert(1 == 0,ExcMessage("Should not have reached here"));
					break;
				}
			}
			break;
		}

		case 70:
		{
			G120::G120<dim> G120(constants.constants,folder_name);

			ExactSolution::G120_PoissonHeat<dim>  G120_PoissonHeat(constants.constants,G120.base_tensorinfo.S_half);

			FEM_Solver::Base_Solver<dim> base_solver("mesh_file_name",
													 "grid",
													 constants.constants,
													 &G120,
													 &G120_PoissonHeat);

			switch(G120.constants.mesh_type)
			{
				case periodic_square:
				{
					base_solver.run_periodic();
					break;
				}
				default:
				{
					Assert(1 == 0,ExcMessage("Should not have reached here"));
					break;
				}
			}
			break;
		}


		default:
		{
			Assert(1 == 0,ExcMessage("Should not have reached here"));
			break;
		}
	}


}