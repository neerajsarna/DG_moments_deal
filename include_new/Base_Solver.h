namespace FEM_Solver
{
	using namespace dealii;

	template<int dim>
	class
	Base_Solver:public MeshGenerator::Base_MeshGenerator<dim>,
				public PeriodicityHandler::Base_PeriodicityHandler<dim>
	{
		typedef MeshWorker::DoFInfo<dim> DoFInfo;
		typedef MeshWorker::IntegrationInfo<dim> CellInfo;
		
		public:
			Base_Solver(const std::string &mesh_file_name,
						const std::string &output_file_name,
						const constant_data &constants,
						EquationGenerator::Base_EquationGenerator<dim> *equation_info);

			const constant_data constants;
			EquationGenerator::Base_EquationGenerator<dim> *system_info;

			// the following order is important or else we run into trouble
			FESystem<dim> finite_element;
			DoFHandler<dim> dof_handler;
			SparsityPattern sparsity_pattern;
			
			
			// A matrix in the final Ax = b
			TrilinosWrappers::SparseMatrix global_matrix;

			// x in final Ax = b
			Vector<double> solution;

			// The rhs in final Ax = b.
        	Vector<double> system_rhs;

        	// running routines for a ring
        	void distribute_dof_allocate_matrix_ring();
        	void run_ring(const unsigned int refine_cycles);

        	// running routines for a periodic box
        	// void distribute_dof_allocate_matrix_periodic_box();
        	// void run_periodic();

        	// assembling routines for meshworker
    		// variables and functions for meshworker
        	void assemble_rhs();

        	void integrate_cell_term (DoFInfo &dinfo,
        		CellInfo &info);
        	void integrate_boundary_term_odd (DoFInfo &dinfo,
        									CellInfo &info);

        	// integrate boundary term using characteristic variables
        	void integrate_boundary_term_char(DoFInfo &dinfo,
        									CellInfo &info);

        	void integrate_face_term (DoFInfo &dinfo1,
        		DoFInfo &dinfo2,
        		CellInfo &info1,
        		CellInfo &info2);


        	void assemble_system_meshworker();

        	// routines for manuel assemblation
        	//     // functions for manuel assemblation
        	// void assemble_system();

        	// void integrate_cell_manuel(FullMatrix<double> &cell_matrix, Vector<double> &cell_rhs,
        	// 							FEValuesBase<dim> &fe_v,  vector<double> &J,
        	// 							vector<Vector<double>> &source_term_value, const typename DoFHandler<dim>::active_cell_iterator &cell);

        	// // integrate the boundary manuelly using odd boundary implementation
        	// void integrate_boundary_manuel_odd(FullMatrix<double> &cell_matrix,
        	// 								Vector<double> &cell_rhs,
        	// 								FEValuesBase<dim> &fe_v,
        	// 								vector<double> &J,
        	// 								Vector<double> &component,
        	// 								const typename DoFHandler<dim>::active_cell_iterator &cell);
        	
        	// // integrate the boundary using characteristic variables
        	// void integrate_boundary_manuel_char(FullMatrix<double> &cell_matrix,
        	// 								Vector<double> &cell_rhs,
        	// 								FEValuesBase<dim> &fe_v,
        	// 								vector<double> &J,
        	// 								Vector<double> &component,
        	// 								const typename DoFHandler<dim>::active_cell_iterator &cell);

        	// void integrate_face_manuel(FullMatrix<double> &u1_v1,
        	// 							FullMatrix<double> &u1_v2,
        	// 							FullMatrix<double> &u2_v1,
        	// 							FullMatrix<double> &u2_v2,
        	// 							FEValuesBase<dim> &fe_v,
        	// 							FEValuesBase<dim> &fe_v_neighbor,
        	// 							vector<double> &J,
        	// 							Vector<double> &component,
        	// 							const typename DoFHandler<dim>::active_cell_iterator &cell);


	};

	template<int dim>
	Base_Solver<dim>::Base_Solver(const std::string &mesh_file_name,
								  const std::string &output_file_name,
								  const constant_data &constants,
								  EquationGenerator::Base_EquationGenerator<dim> *equation_info)
	:
	MeshGenerator::Base_MeshGenerator<dim>(mesh_file_name,output_file_name,constants),
	PeriodicityHandler::Base_PeriodicityHandler<dim>(constants.xl,constants.xr)
	constants(constants),
	system_info(equation_info),
	finite_element(FE_DGQ<dim>(constants.p),constants.nEqn),
	dof_handler(this->triangulation)
	{;}

  	#include "AssembleSystem_Meshworker.h"
  	#include "Run_Ring.h"

}
