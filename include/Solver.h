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
     
        Solver_DG(const unsigned int p,const unsigned int mapping_order,const enum Refinement refinement);
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


        ConvergenceTable convergence_table;
        virtual void print_convergence_table(const string filename);
        virtual void error_evaluation(const Vector<double> solution);
          /*void output_solution_at_vertex(const Triangulation<dim> triangulation,const string filename)const ;
          void output_exact_solution_at_vertex(const Function<dim> exact_solution,const string filename)const ;
          void output_error_at_vertext(const string filename)const;*/

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

        
  };

  template<int dim> void Solver_DG<dim>::prescribe_filenames(output_files &output_file_names,const unsigned int p)
  {
    switch(refinement)
    {
        case global:
        {
          output_file_names.file_for_convergence_tables = this->sub_directory_names[2] + "/convergence_table_global_degree_"
                                                          + to_string(p);
          break;                
        }

        case adaptive:
        {
          output_file_names.file_for_convergence_tables = this->sub_directory_names[2] + "/convergence_table_adaptive_degree_"
                                                          + to_string(p);
          break;
        }

        case adaptive_kelly:
          {
          output_file_names.file_for_convergence_tables = this->sub_directory_names[2] + "/convergence_table_adaptive_kelly_degree_"
                                                          + to_string(p);
          break;
        }      
    }
    
  }

  template<int dim> Solver_DG<dim>::Solver_DG(const unsigned int p,const unsigned int mapping_order,const enum Refinement refinement)
  :
  generate_systemA<dim>(generate_systemA<dim>::LLF),
  finite_element(FE_DGQ<dim>(p),this->nEqn),
  dof_handler(triangulation),
  ngp(p+1),
  ngp_face(p+1),
  refinement(refinement),
  mapping(mapping_order),
  mesh_generation<dim>(mesh_generation<dim>::generate_internal)
  {
    prescribe_filenames(output_file_names,p);
  }

  template<int dim> void Solver_DG<dim>::run(const string mesh_to_read,const unsigned int refine_cycles)
  {

     for (unsigned int i = 0 ; i < refine_cycles; i++)
     {
       cout << "Start of Cycle: " << i << endl;
       if(i == 0)
          mesh_generation<dim>::generate_mesh(triangulation,boundary,mesh_to_read);
       else
          h_adapt();

          distribute_dof_allocate_matrix();

          cout << "assembling the matrix...." << endl;
          assemble_system_meshworker();
          cout << "assemblation completed...." << endl;

          cout << "solving the system...." << endl;
          solve();
          cout << "done solving";

          error_evaluation(solution);

       if (i == refine_cycles - 1)
       {
        string file_for_grid;
        file_for_grid = this->sub_directory_names[0] + "/grid_"+"_DOF_" + to_string(dof_handler.n_dofs());
        //mesh_generation<dim>::print_grid(triangulation,file_for_grid);
        print_convergence_table(output_file_names.file_for_convergence_tables);
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

  #include "assemble.h"
  #include "solve.h"
  #include "h_adaptivity.h"
  #include "error_evaluator.h"
  #include "print_convergence_table.h"
  #include "assemble_system_meshworker.h"
}
