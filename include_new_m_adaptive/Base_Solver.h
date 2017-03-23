namespace FEM_Solver
{
	using namespace dealii;

    
	template<int dim>
	class
	Base_Solver:public MeshGenerator::Base_MeshGenerator<dim>,
				//public PeriodicityHandler::Base_PeriodicityHandler<dim>,
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

            // the exact solution will correspond to a particular system
            ExactSolution::Base_ExactSolution<dim> *base_exactsolution;
			std::vector<Develop_System::System<dim>> system_info;
            MatrixOpt::Base_MatrixOpt matrix_opt;

            
			// the following order is important or else we run into trouble
            // a finite element structure for a single component
            FE_DGQ<dim> fe_dg;

            // block sizes for the finite element objects
            std::vector<int> block_finite_elements;

            // the main finite element object
			hp::FECollection<dim> finite_element;

            // the hp dof handler object
			hp::DoFHandler<dim> dof_handler;
						
			// A matrix in the final Ax = b
			TrilinosWrappers::SparseMatrix global_matrix;

			// x in final Ax = b
			Vector<double> solution;

			// The rhs in final Ax = b.
        	Vector<double> system_rhs;

        	// running routines for a ring
            void construct_block_structure(std::vector<int> &block_structure,const std::vector<int> &nEqn);
            void construct_fe_collection();
            

        	void distribute_dof_allocate_matrix(Vector<double> &error_per_cell, const double tolerance,
                                                const unsigned int present_cycle, const unsigned int total_cycles);

            void allocate_fe_index(Vector<double> &error_per_cell, const double tolerance,
                                 const unsigned int present_cycle, const unsigned int total_cycles);

		    void allocate_vectors(); 
            void run();

        	MappingQ<dim,dim> mapping_basic;
            hp::MappingCollection<dim,dim> mapping;
        	const unsigned int ngp;
        	const unsigned int ngp_face;
            unsigned int integrate_inflow = 0;
            std::vector<unsigned int> dofs_per_cell;
            unsigned int dofs_per_component;

        	// routines for a periodic box
         //    void distribute_dof_allocate_matrix_periodic_box();
        	// void run_periodic();

        	// assembling routines for meshworker
    		// variables and functions for meshworker
        	// void assemble_rhs();

        	// void integrate_cell_term (DoFInfo &dinfo,
        	// 	CellInfo &info);
        	// void integrate_boundary_term_odd (DoFInfo &dinfo,
        	// 								CellInfo &info);

        	// // integrate boundary term using characteristic variables
        	// void integrate_boundary_term_char(DoFInfo &dinfo,
        	// 								CellInfo &info);

        	// void integrate_face_term (DoFInfo &dinfo1,
        	// 	DoFInfo &dinfo2,
        	// 	CellInfo &info1,
        	// 	CellInfo &info2);


        	//void assemble_system_meshworker();

            // different implementations of the same routine
        	//void assemble_system_char();
            void assemble_system_odd();
            //void assemble_system_periodic_char();
            //void assemble_system_periodic_odd();


            // integration for every different element
        	void integrate_cell_manuel(FullMatrix<double> &cell_matrix,
                                        Vector<double> &cell_rhs,
                                        const FEValuesBase<dim> &fe_v,
                                        std::vector<double> &J,
                                        std::vector<Vector<double>> &source_term_value,
                                        std::vector<Vector<double>> &component_to_system,
                                        const unsigned int fe_index);

        	// // integrate the boundary manuelly using odd boundary implementation
        	void integrate_boundary_manuel_odd(FullMatrix<double> &cell_matrix,
                                               Vector<double> &cell_rhs,
                                                const FEValuesBase<dim> &fe_v,
                                                std::vector<double> &J,
                                                std::vector<Vector<double>> &component_to_system,
                                                const unsigned int b_id,
                                                const unsigned int fe_index);
        	
        	// // integrate the boundary using characteristic variables
        	// void integrate_boundary_manuel_char(FullMatrix<double> &cell_matrix,
        	// 								Vector<double> &cell_rhs,
        	// 								FEValuesBase<dim> &fe_v,
        	// 								std::vector<double> &J,
        	// 								Vector<double> &component,
        	// 								const typename DoFHandler<dim>::active_cell_iterator &cell,
         //                                    const unsigned int b_id);

        	void integrate_face_manuel(FullMatrix<double> &u1_v1,
                        FullMatrix<double> &u1_v2,
                        FullMatrix<double> &u2_v1,
                        FullMatrix<double> &u2_v2,
                        const FEValuesBase<dim> &fe_v,
                        const FEValuesBase<dim> &fe_v_neighbor,
                        std::vector<double> &J,
                        const unsigned int fe_index1,
                        const unsigned int fe_index2);

            // Data for post proc
            ConvergenceTable convergence_table;
            std::vector<double> error_per_itr;

         //    // the following quantity is a measure of how good our FE solution satisfies the strong form of the equations
             Vector<double> residual;
             Vector<double> error_per_cell_VelocitySpace;
            
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
	//PeriodicityHandler::Base_PeriodicityHandler<dim>(constants.xl,constants.xr),
	constants(constants),
    nEqn(nEqn),
    nBC(nBC),
    base_exactsolution(exact_solution),
	system_info(equation_info),
	fe_dg(constants.p),
	dof_handler(this->triangulation),
	mapping_basic(constants.mapping_order),
	ngp(constants.p + 1),
	ngp_face(constants.p + 1)
	{
        // we have not implemented any other system till now
        AssertDimension(nEqn.size(),3);

        // we construct the block structure for the finite element object and develop the finite element objects 
        // which will be needed for the present problem
        construct_fe_collection();

        // the total number of finite element objects must be equal to the total number of systems sent to the solver
        AssertDimension(nEqn.size(),finite_element.size());

        // now we initialize the hp mapping object which simply contain the same mapping object repeated a several times
        for (unsigned long int i = 0 ;i < nEqn.size(); i++)
            mapping.push_back(mapping_basic);

        // we now initialize the dofs per cell and the numer of dofs per component. We can do this since we are looking for 
        // the same polynomial degree for every moment system
              // for every moment system, the number of dofs in a given cell are different
      for (unsigned long int i = 0 ;i < nEqn.size() ; i++)
        dofs_per_cell.push_back(finite_element[i].dofs_per_cell);


      // these are the number of degrees of freedom used for one component. They remain the same for all the equations
      dofs_per_component = dofs_per_cell[0]/nEqn[0];


      // the vector dofs per cell should be sorted since the vector number of equations is sorted
      Assert(std::is_sorted(std::begin(dofs_per_cell),std::end(dofs_per_cell)),ExcMessage("The number of dofs are not sorted"));


    }

    //can be used for all the applications apart from periodic box
    template<int dim>
    void
    Base_Solver<dim>::distribute_dof_allocate_matrix(Vector<double> &error_per_cell, const double tolerance,
                                                const unsigned int present_cycle, const unsigned int total_cycles)
    {
        // first we need to allocate the fe_index for all the cells
        allocate_fe_index(error_per_cell, tolerance,
                          present_cycle,total_cycles);

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
    error_per_cell_VelocitySpace.reinit(this->triangulation.n_active_cells());
   }

  


    #include "Integrate_PerCell.h"
  	//#include "AssembleSystem_Meshworker.h"
    #include "AssembleSystem_Manuel.h"
    #include "Run_System.h"     // run a particular test case(Not involving periodic boundaries)
    //#include "Run_Periodic.h"   // A periodic square box
    #include "dof_handler_m_adaptivity.h"
    

}
