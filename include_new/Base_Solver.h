namespace FEM_Solver
{
	using namespace dealii;

    
	template<int dim>
	class
	Base_Solver:public MeshGenerator::Base_MeshGenerator<dim>,
				public PeriodicityHandler::Base_PeriodicityHandler<dim>,
                public NumericalIntegration::Base_NumericalIntegration<dim> 
	{
		typedef MeshWorker::DoFInfo<dim> DoFInfo;
		typedef MeshWorker::IntegrationInfo<dim> CellInfo;

		public:
			Base_Solver(
						const std::string &output_file_name,
						const constant_numerics &constants,
						std::vector<Develop_System::System<dim>> &equation_info,
                        ExactSolution::Base_ExactSolution<dim> *exact_solution,
                        const std::vector<int> &nEqn,
                        const std::vector<int> &nBC);



            DeclException2 (ExcCellCenter, double, double,
                        << "Cell Center = " << arg1 << "Neighbor Center = " << arg2 << "Mesh not lexiographical ");

            DeclException2 (ExcNoElementInSparsity, size_t, size_t,
                        << "Dof-1 " << arg1 << " Dof-2 " << arg2 << " Entry does not exist in sparsity pattern");

			const constant_numerics constants;
            const std::vector<int> nEqn;
            const std::vector<int> nBC;

            ExactSolution::Base_ExactSolution<dim> *base_exactsolution;
			std::vector<Develop_System::System<dim>> system_info;
            MatrixOpt::Base_MatrixOpt matrix_opt;

			// the following order is important or else we run into trouble
			FESystem<dim> finite_element;
			DoFHandler<dim> dof_handler;
			TrilinosWrappers::SparsityPattern sparsity_pattern;
			
			
			// A matrix in the final Ax = b
			TrilinosWrappers::SparseMatrix global_matrix;
      

			// x in final Ax = b
			Vector<double> solution;

			// The rhs in final Ax = b.
        	Vector<double> system_rhs;

        	// running routines for a ring
        	void distribute_dof_allocate_matrix();
		    void allocate_vectors(); 
            void run();

            // routines for multiple systems
            void 

        	MappingQ<dim,dim> mapping;
        	const unsigned int ngp;
        	const unsigned int ngp_face;
            unsigned int integrate_inflow = 0;

        	// routines for a periodic box
            void distribute_dof_allocate_matrix_periodic_box();
        	void run_periodic();

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

            // different implementations of the same routine
        	void assemble_system_char();
            void assemble_system_odd();
            void assemble_system_periodic_char();
            void assemble_system_periodic_odd();


            // integration for every different element
        	void integrate_cell_manuel(FullMatrix<double> &cell_matrix, Vector<double> &cell_rhs,
        								FEValuesBase<dim> &fe_v,  std::vector<double> &J,
        								std::vector<Vector<double>> &source_term_value, 
                                        const typename DoFHandler<dim>::active_cell_iterator &cell);

        	// integrate the boundary manuelly using odd boundary implementation
        	void integrate_boundary_manuel_odd(FullMatrix<double> &cell_matrix,
        									Vector<double> &cell_rhs,
        									FEValuesBase<dim> &fe_v,
        									std::vector<double> &J,
        									Vector<double> &component,
        									const typename DoFHandler<dim>::active_cell_iterator &cell,
                                            const unsigned int b_id);
        	
        	// integrate the boundary using characteristic variables
        	void integrate_boundary_manuel_char(FullMatrix<double> &cell_matrix,
        									Vector<double> &cell_rhs,
        									FEValuesBase<dim> &fe_v,
        									std::vector<double> &J,
        									Vector<double> &component,
        									const typename DoFHandler<dim>::active_cell_iterator &cell,
                                            const unsigned int b_id);

        	void integrate_face_manuel(FullMatrix<double> &u1_v1,
                                       FullMatrix<double> &u1_v2,
                                       FullMatrix<double> &u2_v1,
                                       FullMatrix<double> &u2_v2,
        								FEValuesBase<dim> &fe_v,
        								FEValuesBase<dim> &fe_v_neighbor,
        								std::vector<double> &J,
        								Vector<double> &component,
        								const typename DoFHandler<dim>::active_cell_iterator &cell);

            // Data for post proc
            ConvergenceTable convergence_table;
            std::vector<double> error_per_itr;

            // the following quantity is a measure of how good our FE solution satisfies the strong form of the equations
            Vector<double> residual;
            
	};

	template<int dim>
	Base_Solver<dim>::Base_Solver(
								  const std::string &output_file_name,
								  const constant_numerics &constants,
								  std::vector<Develop_System::System<dim>> &equation_info,
                                  ExactSolution::Base_ExactSolution<dim> *exact_solution,
                                  const std::vector<int> &nEqn,
                                  const std::vector<int> &nBC)
	:
	MeshGenerator::Base_MeshGenerator<dim>(output_file_name,constants),
	PeriodicityHandler::Base_PeriodicityHandler<dim>(constants.xl,constants.xr),
	constants(constants),
    nEqn(nEqn),
    nBC(nBC),
    base_exactsolution(exact_solution),
	system_info(equation_info),
	finite_element(FE_DGQ<dim>(constants.p),nEqn[0]),
	dof_handler(this->triangulation),
	mapping(constants.mapping_order),
	ngp(constants.p + 1),
	ngp_face(constants.p + 1)
	{
        // we have not implemented any other system till now
        AssertDimension(nEqn.size(),1);
    }

    //can be used for all the applications apart from periodic box
    template<int dim>
    void
    Base_Solver<dim>::distribute_dof_allocate_matrix()
    {
        dof_handler.distribute_dofs(finite_element);

        // the vector which stores the residual 
        DynamicSparsityPattern dsp(dof_handler.n_dofs(),dof_handler.n_dofs());


	//std::cout << "making flux sparsity pattern " << std::endl;
        DoFTools::make_flux_sparsity_pattern (dof_handler, dsp);
	
        global_matrix.reinit(dsp);   
	
    }

   template<int dim>
   void 
   Base_Solver<dim>::allocate_vectors()
   {
	solution.reinit(dof_handler.n_dofs());
	system_rhs.reinit(dof_handler.n_dofs());
    residual.reinit(dof_handler.n_dofs());
   }

  


    #include "Integrate_PerCell.h"
  	#include "AssembleSystem_Meshworker.h"
    #include "AssembleSystem_Manuel.h"
    #include "Run_System.h"     // run a particular test case(Not involving periodic boundaries)
    #include "Run_Periodic.h"   // A periodic square box
    

}
