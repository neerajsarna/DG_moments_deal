using namespace dealii;
#include "global_matrix_G20_square_inflow_outflow.h"

// we consider a square domain with boundaries and check the assembly 
TEST(AsseblySingleSystem,HandlesSingleSystem)
{
		const unsigned int dim = 2;

		// not implemented for any other dimension
		AssertDimension(dim,2);

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);

	// 	// we construct a vector of all the system data being considered 
		std::vector<Develop_System::System<dim>> System;

		// initialize the vector containing all the systems
		// initialize the vector containing all the systems
		for(int i = 0 ; i < constants.constants_sys.total_systems ; i++)
			 System.push_back(Develop_System::System<dim>(constants.constants_num,constants.constants_sys.nEqn[i],
			 											constants.constants_sys.nBC[i],constants.constants_sys.Ntensors[i],folder_name));


		for (int i = 0 ; i < constants.constants_sys.total_systems ; i++)
			System[i].initialize_system();


		// the exact solution can only be created for one of the systems
		ExactSolution::ExactSolution_Dummy<dim>  exactsolution_dummy(constants.constants_num,System[0].base_tensorinfo.S_half,
																			 constants.constants_sys.nEqn[0],constants.constants_sys.Ntensors[0]);

		FEM_Solver::Base_Solver<dim> base_solver("grid",
											 	 constants.constants_num,
											 	 System,
											 	 &exactsolution_dummy,
											 	 constants.constants_sys.nEqn,
											 	 constants.constants_sys.nBC);

		// following are the specifications which should be met when running this particular test
		Assert(constants.constants_num.problem_type == heat_conduction,ExcNotImplemented());
		Assert(constants.constants_num.mesh_type == square_domain,ExcNotImplemented());
		AssertDimension(constants.constants_num.refine_cycles,1);
		AssertDimension(constants.constants_num.initial_refinement,1);
		AssertDimension(constants.constants_sys.Ntensors[0],6);
		AssertDimension(constants.constants_sys.total_systems,1);

		// we consider a single cell for simplicity
		AssertDimension(constants.constants_num.part_x,1);
		AssertDimension(constants.constants_num.part_y,1);

		// we now change the boundary ids of the cell and replace them by inflow and outflow boundaries
		typename Triangulation<dim>::cell_iterator cell = base_solver.triangulation.begin(),
                                                    endc = base_solver.triangulation.end();

        for (; cell != endc ; cell++)
        {
                for (unsigned int face = 0 ; face < GeometryInfo<dim>::faces_per_cell ; face++)
              
                  if (cell->face(face)->at_boundary())
                  { 
                    double x_cord = cell->face(face)->center()(0);
                    double y_cord = cell->face(face)->center()(1);

                    // Following is the boundary ids:
                    // Left Wall = 0
                    // Bottom Wall = 1
                    // Right Wall = 2
                    // Top Wall = 3
                    // the periodic faces get the id 100 and 101
                    // right edge
                    
                    if (x_cord == constants.constants_num.xr)
                      cell->face(face)->set_boundary_id(1);

                    // left edge
                    if (x_cord == constants.constants_num.xl)
                      cell->face(face)->set_boundary_id(101);

                    // this is the bottom wall
                    if (y_cord == constants.constants_num.yb)
                      cell->face(face)->set_boundary_id(0);

                    // top edge, This is the second wall
                    if (y_cord == constants.constants_num.yt)
                      cell->face(face)->set_boundary_id(102);
                   }
        }



        base_solver.distribute_dof_allocate_matrix();
		base_solver.allocate_vectors();
		base_solver.print_mesh_info();
		

		// the following routine assembles
		std::cout << "Assembling" << std::endl;
		switch (constants.constants_num.assembly_type)
		{
			case meshworker:
			{
				AssertThrow(1 == 0,ExcMessage("Meshworker not available for hp dof_handler. "));
				std::cout << "Using meshwoker " << std::endl;
				//assemble_system_meshworker();
				break;		
			}
			case manuel:
			{
				//AssertThrow(1 == 0,ExcMessage("Should not have reached here"));
				switch(constants.constants_num.bc_type)
				{
					case odd:
					{
						std::cout << "Using odd boundary " << std::endl;
						base_solver.assemble_system_odd();
						break;
					}

					case characteristic:
					default:
					{
						AssertThrow(1 == 0, ExcMessage("Should not have reached here"));
						break;
					}
		
				}
				break;
			}

			default:
			{
				AssertThrow(1==0,ExcMessage("Should not have reached here"));
				break;
			}
		}

		Full_matrix global_matrix;
		global_matrix.resize(base_solver.global_matrix.m(),base_solver.global_matrix.n());
		develop_global_matrix(global_matrix);
		Compare_Float_Mat(base_solver.global_matrix,global_matrix);

}


