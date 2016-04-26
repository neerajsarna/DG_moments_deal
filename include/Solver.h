namespace SolverDG
{
  using namespace dealii;
  using namespace std;
  using namespace Mesh_Handler;
  using namespace EquationGenerator;
  using namespace PostProc;
  using namespace ExactSolution;

  template<int dim> class Solver_DG:protected generate_systemA<dim>,
                                    protected mesh_generation<dim>,
                                    protected Base_PostProc<dim>
  {

    typedef MeshWorker::DoFInfo<dim> DoFInfo;
    typedef MeshWorker::IntegrationInfo<dim> CellInfo;


    public:
        enum Refinement
          {global,adaptive,adaptive_kelly};

        Refinement refinement;

        Solver_DG(const unsigned int p,const unsigned int mapping_order,
                  const enum Refinement refinement,const Base_ExactSolution<dim> *exact_solution);

        void run(const string mesh_to_read,const unsigned int refine_cycles);

    private:
        SphericalManifold<dim> boundary;
        Triangulation<dim> triangulation;
        FESystem<dim> finite_element;
        DoFHandler<dim> dof_handler;

        SparsityPattern sparsity_pattern;
        TrilinosWrappers::SparseMatrix global_matrix;
     
        Vector<double> solution;
        Vector<double> system_rhs;

        const MappingQ<dim> mapping;
        const unsigned int ngp;
        const unsigned int ngp_face;

        void distribute_dof_allocate_matrix(); 
        void assemble_system();
        void solve();
        void h_adapt();

  
    // variables and functions for meshworker
      void assemble_rhs();

      void integrate_cell_term (DoFInfo &dinfo,
                                   CellInfo &info);
      void integrate_boundary_term (DoFInfo &dinfo,
                                       CellInfo &info);
      void integrate_face_term (DoFInfo &dinfo1,
                                   DoFInfo &dinfo2,
                                   CellInfo &info1,
                                   CellInfo &info2);

     void assemble_system_meshworker();


    // variable for post processing
        const Base_ExactSolution<dim> *exact_solution;
        ConvergenceTable convergence_table;
        virtual void print_convergence_table(const string filename);
        virtual void error_evaluation(const Vector<double> solution);
        virtual void output_solution_details(const Triangulation<dim> &triangulation,const string file_solution
                                                ,const string file_exact,const string file_error)const ;

	
      struct output_files
      {
        string file_for_convergence_tables;
        string file_for_grid;
        string file_for_num_solution;
        string file_for_exact_solution;
        string file_for_error;
      };

      output_files output_file_names;

      void prescribe_filenames(output_files &output_file_names,const unsigned int p);

        
  };

  template<int dim> void Solver_DG<dim>::prescribe_filenames(output_files &output_file_names,const unsigned int p)
  {
    switch(refinement)
    {
        case global:
        {
          output_file_names.file_for_convergence_tables = this->sub_directory_names[2] + "/convergence_table_global_degree_"
                                                          + to_string(p);

          output_file_names.file_for_num_solution = this->sub_directory_names[1] + "/numerical_solution_global_degree_"
                                                          + to_string(p)+"_DOF_"+to_string(dof_handler.n_dofs());


          output_file_names.file_for_exact_solution = this->sub_directory_names[1] + "/exact_solution_global_degree_"
                                                          + to_string(p)+"_DOF_"+to_string(dof_handler.n_dofs());


          output_file_names.file_for_error = this->sub_directory_names[1] + "/error_global_degree_"
                                                          + to_string(p)+"_DOF_"+to_string(dof_handler.n_dofs());

          break;                
        }

        case adaptive:
        {
          output_file_names.file_for_convergence_tables = this->sub_directory_names[2] + "/convergence_table_adaptive_degree_"
                                                          + to_string(p);

          output_file_names.file_for_num_solution = this->sub_directory_names[1] + "/numerical_solution_adaptive_degree_"
                                                          + to_string(p)+"_DOF_"+to_string(dof_handler.n_dofs());

          output_file_names.file_for_exact_solution = this->sub_directory_names[1] + "/exact_solution_adaptive_degree_"
                                                          + to_string(p)+"_DOF_"+to_string(dof_handler.n_dofs());


          output_file_names.file_for_error = this->sub_directory_names[1] + "/error_adaptive_degree_"
                                                          + to_string(p)+"_DOF_"+to_string(dof_handler.n_dofs());

          break;
        }

        case adaptive_kelly:
          {
          output_file_names.file_for_convergence_tables = this->sub_directory_names[2] + "/convergence_table_adaptive_kelly_degree_"
                                                          + to_string(p);
          
          output_file_names.file_for_num_solution = this->sub_directory_names[1] + "/numerical_solution_adaptive_kelly_degree_"
                                                          + to_string(p)+"_DOF_"+to_string(dof_handler.n_dofs());

          output_file_names.file_for_exact_solution = this->sub_directory_names[1] + "/exact_solution_adaptive_kelly_degree_"
                                                          + to_string(p)+"_DOF_"+to_string(dof_handler.n_dofs());

          output_file_names.file_for_error = this->sub_directory_names[1] + "/error_adaptive_kelly_degree_"
                                                          + to_string(p)+"_DOF_"+to_string(dof_handler.n_dofs());

          break;
        }      
    }
    
  }

  template<int dim> Solver_DG<dim>::Solver_DG(const unsigned int p,const unsigned int mapping_order,
                                              const enum Refinement refinement,const Base_ExactSolution<dim> *exact_solution)
  :
  generate_systemA<dim>(generate_systemA<dim>::Upwind),
  finite_element(FE_DGQ<dim>(p),this->nEqn),
  dof_handler(triangulation),
  ngp(p+1),
  ngp_face(p+1),
  refinement(refinement),
  mapping(mapping_order),
  mesh_generation<dim>(mesh_generation<dim>::generate_internal),
  exact_solution(exact_solution)
  {
  }

  template<int dim> void Solver_DG<dim>::run(const string mesh_to_read,const unsigned int refine_cycles)
  {

    TimerOutput timer (std::cout, TimerOutput::summary,
                   TimerOutput::wall_times);
    
    cout << "Solving for: " << this->nEqn << " equations " << endl;
          
     for (unsigned int i = 0 ; i < refine_cycles; i++)
     {
       cout << "Start of Cycle: " << i << endl;
       if(i == 0)
	{
	  timer.enter_subsection("mesh_generation");
          mesh_generation<dim>::generate_mesh(triangulation,boundary,mesh_to_read);
	  cout << "no of cells in the initial mesh" << triangulation.n_cells() << endl;  
     	  timer.leave_subsection();
	}
	else
          h_adapt();
	
	  timer.enter_subsection("distributing dof");
          distribute_dof_allocate_matrix();
	  timer.leave_subsection();

	cout << "Total #DOF: " << dof_handler.n_dofs() << endl;

	  cout << "assembling the matrix...." << endl;
	  timer.enter_subsection("assemblation");
	  //assemble_system();
          assemble_system_meshworker();
	  timer.leave_subsection();
          cout << "assemblation completed...." << endl;

          cout << "solving the system...." << endl;
	  timer.enter_subsection("solving the system");
          solve();
	  timer.leave_subsection();
          cout << "done solving";

	  timer.enter_subsection("error evaluation");
          error_evaluation(solution);
	  timer.leave_subsection();

       if (i == refine_cycles - 1)
       {

        prescribe_filenames(output_file_names,finite_element.degree);
        string file_for_grid;
        file_for_grid = this->sub_directory_names[0] + "/grid_"+"_DOF_" + to_string(dof_handler.n_dofs());
        //mesh_generation<dim>::print_grid(triangulation,file_for_grid);
        print_convergence_table(output_file_names.file_for_convergence_tables);
        //output_solution_details(triangulation,output_file_names.file_for_num_solution,
        //                        output_file_names.file_for_exact_solution,
        //                        output_file_names.file_for_error);
       }
     }

  }

  template<int dim> void Solver_DG<dim>::distribute_dof_allocate_matrix()
  {
    dof_handler.distribute_dofs(finite_element);

    DynamicSparsityPattern dsp(dof_handler.n_dofs(),dof_handler.n_dofs());
    DoFTools::make_flux_sparsity_pattern (dof_handler, dsp);
 
    sparsity_pattern.copy_from(dsp);
    
    global_matrix.reinit (sparsity_pattern);
    
    solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());


  }

  #include "solve.h"
  #include "h_adaptivity.h"
  #include "error_evaluator.h"
  #include "print_convergence_table.h"
//  #include "assemble.h"
  #include "assemble_system_meshworker.h"
  #include "output_routines.h"
}
