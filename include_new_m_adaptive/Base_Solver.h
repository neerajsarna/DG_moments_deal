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

            // the data structure which store the constants corresponding to numerics.
			const constant_numerics constants;

            // a vector containing all the number of equations in the ssytems
            const std::vector<int> nEqn;

            // a vector containing the number of boundary conditions corresponding to different systems begin 
            // considered.
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

        	// to construct our fe collection we need the block structure of the elementary finite element
            // objects.
            void construct_block_structure(std::vector<int> &block_structure,const std::vector<int> &nEqn);

            // using the block structure constructed in the previous routine, in the following routine we develop 
            // the fe collection object.
            void construct_fe_collection();
            

            // distribute the dofs and allocate the matrices
        	  void distribute_dof_allocate_matrix(const unsigned int present_cycle,const unsigned int total_cycles);

            // void distribute_dof_residual_computation();
            // void develop_solution_max_moments();
            // void compute_residual();   // the following routine computes the residual 
            void compute_equilibrium_deviation();

            // sets the tolerance band for difference moment systems
            void set_tolerance_bands();


            void allocate_fe_index(const unsigned int present_cycle,const unsigned int total_cycles);

		        void allocate_vectors(); 
            void run();

        	MappingQ<dim,dim> mapping_basic;
            hp::MappingCollection<dim,dim> mapping;
        	const unsigned int ngp;
        	const unsigned int ngp_face;
            unsigned int integrate_inflow = 0;

            // we store the number of degrees of freedom in every cell. The assumption is that all the components 
            // of all the moment system are solved with the same polynomial degree.
            std::vector<unsigned int> dofs_per_cell;

            // fe_index of the maximum moment system which has to be solver in the computation
            unsigned int max_fe_index;

            // tolerance for the residual, if we find that the residual in a particular cell is greater than this value
            // then we refine it. 
            std::vector<double> VelocitySpace_error_tolerance;

            // the total number of degrees of freedom per component will remain the same for every cell since we do not
            // have p-adaptivity presently
            unsigned int dofs_per_component;

            void assemble_system_odd();


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

         // the following quantity is a measure of how good our FE solution satisfies the weak form of the equations
             Vector<double> residual;

         // Data objects needed for refinement in the velocity space. In the residual we compute the 
         // dofs of the error and in the error_per_cell object we store the error per cell generated by these
         // degress of freedom.
             Vector<double> VelocitySpace_residual;
             Vector<double> VelocitySpace_error_per_cell;



                 
            
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
	ngp_face(constants.p + 1),
  // maximum fe index possible in the system
  max_fe_index(nEqn.size()-1)
	{

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
    Base_Solver<dim>::distribute_dof_allocate_matrix(const unsigned int present_cycle,const unsigned int total_cycles)
    {
        // first we need to allocate the fe_index for all the cells
        allocate_fe_index(present_cycle,total_cycles);

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
    VelocitySpace_error_per_cell.reinit(this->triangulation.n_active_cells());
   }

  


    #include "Integrate_PerCell.h"
  	//#include "AssembleSystem_Meshworker.h"
    #include "AssembleSystem_Manuel.h"
    #include "Run_System.h"     // run a particular test case(Not involving periodic boundaries)
    //#include "Run_Periodic.h"   // A periodic square box
    #include "dof_handler_m_adaptivity.h"

    //computes the dofs of a residual through numerical integration
//    #include "compute_residual.h"    
    // this is the deviation of the distribution function from the equilibrium
    #include "compute_equilibrium_deviation.h"

}
