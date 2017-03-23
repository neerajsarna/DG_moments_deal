using namespace dealii;
#include "global_matrix_1D_poisson_madaptive.h"

// we compare the global matrix for the 1D case for an m_adaptive system
TEST(AsseblyMAdaptivity,HandlesMAdaptivity)
{
		const unsigned int dim = 1;

		// not implemented for any other dimension
		AssertDimension(dim,1);

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
		Assert(constants.constants_num.mesh_type == line,ExcNotImplemented());
		AssertDimension(constants.constants_num.refine_cycles,1);
		AssertDimension(constants.constants_sys.Ntensors[0],6);
		AssertDimension(constants.constants_sys.Ntensors[1],8);
		AssertDimension(constants.constants_sys.total_systems,2);

		// for simplicity we consider only two cells
		AssertDimension(constants.constants_num.initial_refinement,1);


		// we now change the boundary ids of the cell and replace them by inflow and outflow boundaries
		typename Triangulation<dim>::cell_iterator cell = base_solver.triangulation.begin(),
                                                    endc = base_solver.triangulation.end();

        for (; cell != endc ; cell++)
        {
                for (unsigned int face = 0 ; face < GeometryInfo<dim>::faces_per_cell ; face++)
              
                  if (cell->face(face)->at_boundary())
                  { 
                    double x_cord = cell->face(face)->center()(0);
                   
                    // Following is the boundary ids:
                    // Left Wall = 0
                    // Bottom Wall = 1
                    // Right Wall = 2
                    // Top Wall = 3
                    // the periodic faces get the id 100 and 101
                    // right edge
                    
                    if (x_cord == constants.constants_num.xr)
                      cell->face(face)->set_boundary_id(0);

                    // left edge
                    if (x_cord == constants.constants_num.xl)
                      cell->face(face)->set_boundary_id(101);

                   }
        }



        base_solver.distribute_dof_allocate_matrix();
		base_solver.allocate_vectors();
		base_solver.print_mesh_info();

		// now we check whether the allocation of fe index is as required or not
		typename hp::DoFHandler<dim>::active_cell_iterator cell_check = base_solver.dof_handler.begin_active(),
														  endc_check = base_solver.dof_handler.end();


		for (; cell_check != endc_check ; cell_check++)
		{
			if (cell_check->index() == 0)
				AssertDimension(cell_check->active_fe_index(),0);

			if (cell_check->index() == 1)
				AssertDimension(cell_check->active_fe_index(),1);
		}
		

		// the following routine assembles
		base_solver.assemble_system_odd();

		Full_matrix global_matrix;
		global_matrix.resize(base_solver.global_matrix.m(),base_solver.global_matrix.n());
		develop_global_matrix(global_matrix);
		Compare_Float_Mat(base_solver.global_matrix,global_matrix);

}