namespace SolverDG
{
  using namespace dealii;
  using namespace std;
  using namespace Mesh_Handler;
  using namespace PostProc;
  
  /*void pardiso(void *, MKL_INT *, MKL_INT *,MKL_INT *, MKL_INT *,
               MKL_INT *,void *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *,
               MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *);*/
 
  template<int num_flux,int dim> class Solver_DG:protected mesh_generation<dim>,
                                                                protected Base_PostProc<dim>,
                                                                protected Base_Basics
  {

    typedef MeshWorker::DoFInfo<dim> DoFInfo;
    typedef MeshWorker::IntegrationInfo<dim> CellInfo;


    public:

      const mesh_data mesh_info;
	
      unsigned int count_meshworker;
      unsigned int count_manuel;

      enum Solver_Type
          {Trilinos_Direct,Trilinos_GMRES,Pardiso};

        Refinement refinement;

        BC_type bc_type;

        Assembly_Type assembly_type;

        Solver_DG(const numerical_data numerical_constants,
                   EquationGenerator::Base_EquationGenerator<num_flux,dim> *system_of_equations,
                  const ExactSolution::Base_ExactSolution<dim> *exact_solution,
                  const nEqn_data num_equations,
                  physical_data &physical_constants,
                  string &output_dir,
                  mesh_data &mesh_info);


        Solver_DG(const numerical_data numerical_constants,
                  EquationGenerator::Base_EquationGenerator<num_flux,dim> *system_of_equations,
                  const nEqn_data num_equations,
                  physical_data &physical_constants,
                  string &output_dir,
                  mesh_data &mesh_info);

        void run(const unsigned int refine_cycles);

    private:
        EquationGenerator::Base_EquationGenerator<num_flux,dim> *equation_system_data;
        const ExactSolution::Base_ExactSolution<dim> *exact_solution;

        const unsigned int solve_system;
        const unsigned int no_of_BC;
        const unsigned int nEqn;

        SphericalManifold<dim> boundary;
        Triangulation<dim> triangulation;
         GridIn<dim> gridin;
        FESystem<dim> finite_element;
        DoFHandler<dim> dof_handler;

        SparsityPattern sparsity_pattern;
        TrilinosWrappers::SparseMatrix global_matrix;

        // data for periodicity
        vector<GridTools::PeriodicFacePair<typename DoFHandler<dim>::cell_iterator> > periodicity_vector;

        // stores the y coord of the cell center and cell iterator
        map<double, typename DoFHandler<dim>::cell_iterator> set_xl;
        map<double, typename DoFHandler<dim>::cell_iterator> set_xr;

        // stores the y coord of the cell center and the local face number
        map<double, unsigned int> set_xl_face;
        map<double, unsigned int> set_xr_face;

        Vector<double> solution;
        Vector<double> system_rhs;

        const MappingQ<dim> mapping;
        const unsigned int ngp;
        const unsigned int ngp_face;

        // functions for adding periodicity
        void add_periodic_sparsity(DynamicSparsityPattern &dsp,
              vector<GridTools::PeriodicFacePair<typename DoFHandler<dim>::cell_iterator>> &periodicity_vector);

        void divide_periodicity(vector<GridTools::PeriodicFacePair<typename DoFHandler<dim>::cell_iterator> > &periodicity_vector,
             map<double, typename DoFHandler<dim>::cell_iterator> &xminus1_set,
             map<double, typename DoFHandler<dim>::cell_iterator> &xplus1_set,
             map<double, unsigned int> &set_xl_face,
             map<double, unsigned int> &set_xr_face,
             const double xl,
             const double xr);

        typename DoFHandler<dim>::cell_iterator  get_periodic_neighbor(const double xcord_face,
                                                                       const double ycord_center) const;

        unsigned int get_periodic_neighbor_face(const double xcord_face,
                                       const double ycord_center) const;

        void distribute_dof_allocate_matrix(); 
        void solve(const enum Solver_Type solver_type);

        void PardisoSolve(MKL_INT *ia,MKL_INT *ja,double *a,double *b,double *x,MKL_INT n);
        
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

     // functions for manuel assemblation
     void assemble_system();

     void integrate_cell_manuel(FullMatrix<double> &cell_matrix, Vector<double> &cell_rhs,
                                FEValuesBase<dim> &fe_v,  vector<double> &J,
                                vector<Vector<double>> &source_term_value, const typename DoFHandler<dim>::active_cell_iterator &cell);

     void integrate_boundary_manuel(FullMatrix<double> &cell_matrix,
                                    Vector<double> &cell_rhs,
                                    FEValuesBase<dim> &fe_v,
                                    vector<double> &J,
                                    Vector<double> &component,
                                    const typename DoFHandler<dim>::active_cell_iterator &cell);

     void integrate_face_manuel(FullMatrix<double> &u1_v1,
                                FullMatrix<double> &u1_v2,
                                FullMatrix<double> &u2_v1,
                                FullMatrix<double> &u2_v2,
                                FEValuesBase<dim> &fe_v,
                                FEValuesBase<dim> &fe_v_neighbor,
                                vector<double> &J,
                                Vector<double> &component,
                                const typename DoFHandler<dim>::active_cell_iterator &cell);


    // variable for post processing
        ConvergenceTable convergence_table;
        virtual void print_convergence_table(const string filename);
        virtual void error_evaluation(const Vector<double> &solution);
        virtual void evaluate_norms(const Vector<double> &solution);
        virtual void output_solution_details(const Triangulation<dim> &triangulation,const string file_solution
                                                ,const string file_exact,const string file_error)const ;

       virtual void output_solution_details(const Triangulation<dim> &triangulation,const string file_solution)const;	
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

  template<int num_flux,int dim> void 
  Solver_DG<num_flux,dim>
  ::prescribe_filenames(output_files &output_file_names,const unsigned int p)
  {
    switch(refinement)
    {
        case global:
        {
          output_file_names.file_for_convergence_tables = sub_directory_names[2] + "/convergence_table_global_degree_"
                                                          + to_string(p);

          output_file_names.file_for_num_solution = sub_directory_names[1] + "/numerical_solution_global_degree_"
                                                          + to_string(p)+"_DOF_"+to_string(dof_handler.n_dofs());


          output_file_names.file_for_exact_solution = sub_directory_names[1] + "/exact_solution_global_degree_"
                                                          + to_string(p)+"_DOF_"+to_string(dof_handler.n_dofs());


          output_file_names.file_for_error = sub_directory_names[1] + "/error_global_degree_"
                                                          + to_string(p)+"_DOF_"+to_string(dof_handler.n_dofs());

          break;                
        }

        case adaptive:
        {
          output_file_names.file_for_convergence_tables = sub_directory_names[2] + "/convergence_table_adaptive_degree_"
                                                          + to_string(p);

          output_file_names.file_for_num_solution = sub_directory_names[1] + "/numerical_solution_adaptive_degree_"
                                                          + to_string(p)+"_DOF_"+to_string(dof_handler.n_dofs());

          output_file_names.file_for_exact_solution = sub_directory_names[1] + "/exact_solution_adaptive_degree_"
                                                          + to_string(p)+"_DOF_"+to_string(dof_handler.n_dofs());


          output_file_names.file_for_error = sub_directory_names[1] + "/error_adaptive_degree_"
                                                          + to_string(p)+"_DOF_"+to_string(dof_handler.n_dofs());

          break;
        }

        case adaptive_kelly:
          {
          output_file_names.file_for_convergence_tables = sub_directory_names[2] + "/convergence_table_adaptive_kelly_degree_"
                                                          + to_string(p);
          
          output_file_names.file_for_num_solution = sub_directory_names[1] + "/numerical_solution_adaptive_kelly_degree_"
                                                          + to_string(p)+"_DOF_"+to_string(dof_handler.n_dofs());

          output_file_names.file_for_exact_solution = sub_directory_names[1] + "/exact_solution_adaptive_kelly_degree_"
                                                          + to_string(p)+"_DOF_"+to_string(dof_handler.n_dofs());

          output_file_names.file_for_error = sub_directory_names[1] + "/error_adaptive_kelly_degree_"
                                                          + to_string(p)+"_DOF_"+to_string(dof_handler.n_dofs());

          break;
        }      
    }
    
  }

  template<int num_flux,int dim> 
  Solver_DG<num_flux,dim>::Solver_DG(const numerical_data numerical_constants,
                                                  EquationGenerator::Base_EquationGenerator<num_flux,dim> *system_of_equations,
                                                  const ExactSolution::Base_ExactSolution<dim> *exact_solution,
                                                  const nEqn_data num_equations,
                                                  physical_data &physical_constants,
                                                  string &output_dir,
                                                  mesh_data &mesh_info)
  :
  mesh_generation<dim>(mesh_info),
  Base_Basics(physical_constants,output_dir),
  mesh_info(mesh_info),
  refinement(numerical_constants.refinement),
  bc_type(num_equations.bc_type),
  assembly_type(numerical_constants.assembly_type),
  exact_solution(exact_solution),
  solve_system(num_equations.system_to_solve),
  no_of_BC(num_equations.nBC[num_equations.system_to_solve]),
  nEqn(num_equations.total_nEqn[num_equations.system_to_solve]),
  finite_element(FE_DGQ<dim>(numerical_constants.p),num_equations.total_nEqn[num_equations.system_to_solve]),
  dof_handler(triangulation),
  mapping(numerical_constants.mapping_order),
  ngp(numerical_constants.p+1),
  ngp_face(numerical_constants.p+1)
  {
    equation_system_data = system_of_equations;
  }

  // constructor for implementations where the exact solution is not known
  template<int num_flux,int dim>
  Solver_DG<num_flux,dim>::Solver_DG(const numerical_data numerical_constants,
                                                  EquationGenerator::Base_EquationGenerator<num_flux,dim> *system_of_equations,
                                                  const nEqn_data num_equations,
                                                  physical_data &physical_constants,
                                                  string &output_dir,
                                                  mesh_data &mesh_info)
  :
  mesh_generation<dim>(mesh_info),
  Base_Basics(physical_constants,output_dir),
  mesh_info(mesh_info),
  refinement(numerical_constants.refinement),
  bc_type(num_equations.bc_type),
  assembly_type(numerical_constants.assembly_type),
  solve_system(num_equations.system_to_solve),
  no_of_BC(num_equations.nBC[num_equations.system_to_solve]),
  nEqn(num_equations.total_nEqn[num_equations.system_to_solve]),
  finite_element(FE_DGQ<dim>(numerical_constants.p),num_equations.total_nEqn[num_equations.system_to_solve]),
  dof_handler(triangulation),
  mapping(numerical_constants.mapping_order),
  ngp(numerical_constants.p+1),
  ngp_face(numerical_constants.p+1)
  {
    equation_system_data = system_of_equations;
  }


  template<int num_flux,int dim> 
  void 
  Solver_DG<num_flux,dim>::run(const unsigned int refine_cycles)
                                                      
  {
    switch(mesh_info.mesh_type)
    {
      case ring:
        break;

      case periodic_square:
      {
        Assert(refine_cycles == 1,ExcNotImplemented());
        break;
      }

      default:
        Assert(1 == 0,ExcNotImplemented());
    }


    TimerOutput timer (std::cout, TimerOutput::summary,
                   TimerOutput::wall_times);

    cout << "Solving for: " << nEqn << " equations " << endl;
          

    for (unsigned int i = 0 ; i < refine_cycles ; i ++)
    {

      if (i == 0)
      {
        timer.enter_subsection("mesh_generation");
        mesh_generation<dim>::generate_mesh(triangulation,boundary,gridin);
        cout << "no of cells in the initial mesh" << triangulation.n_cells() << endl;  
        timer.leave_subsection();
      }
      else
        h_adapt();


      switch(mesh_info.mesh_type)
      {
        case ring:
          break;

        case periodic_square:
          {
                    // set a periodic link between bid 2 and bid 0
            GridTools::collect_periodic_faces(dof_handler, 2,
                                             0,
                                             0,
                                             periodicity_vector);

        // develops a map bewteen the boundary location and the respective cells
            divide_periodicity(periodicity_vector,
                               set_xl, 
                               set_xr,
                               set_xl_face,
                               set_xr_face,
                               mesh_info.xl,
                               mesh_info.xr);
            break;
          }
      }


      timer.enter_subsection("distributing dof");
      distribute_dof_allocate_matrix();
      timer.leave_subsection();

      cout << "Total #DOF: " << dof_handler.n_dofs() << endl;

      switch(assembly_type)
      {
        case meshworker:
        {
          cout << "assembling the matrix using meshworker...." << endl;
          // meshworker cannot presently handle periodic boundary conditions
          Assert(mesh_info.mesh_type == periodic_square, ExcNotImplemented());
          timer.enter_subsection("assemblation");
          assemble_system_meshworker();
          timer.leave_subsection();
          cout << "assemblation completed...." << endl;
          break;
        }
        case manuel:
        {
          cout << "assembling the matrix manually " << endl;
          assemble_system();
	        cout << "finished assembling the matrix " << endl;
          break;
        }
      }

      //fflush(stdout);
      //Assert(1 == 0, ExcMessage("debugging from Solver"));
     // print_dealii_sparse(global_matrix);

      cout << "solving the system...." << endl;
      timer.enter_subsection("solving the system");
      solve(Pardiso);
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
        mesh_generation<dim>::print_grid(triangulation,file_for_grid);
        //    print_convergence_table(output_file_names.file_for_convergence_tables);  
            output_solution_details(triangulation,output_file_names.file_for_num_solution,
                                    output_file_names.file_for_exact_solution,
                                    output_file_names.file_for_error);


      }
    }

  }



  template<int num_flux,int dim> 
  void 
  Solver_DG<num_flux,dim>::distribute_dof_allocate_matrix()
  {
    dof_handler.distribute_dofs(finite_element);

    DynamicSparsityPattern dsp(dof_handler.n_dofs(),dof_handler.n_dofs());


    DoFTools::make_flux_sparsity_pattern (dof_handler, dsp);
  
    switch(mesh_info.mesh_type)
    {
        case ring:
          break;

        // add coupling between the periodic cells
        case periodic_square:
        {
          add_periodic_sparsity(dsp,periodicity_vector);
          break;
        }
        

    }

    sparsity_pattern.copy_from(dsp);
 
    global_matrix.reinit(sparsity_pattern);   
 
    solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());


  }

  #include "periodicity.h"
  #include "solve.h"
  #include "h_adaptivity.h"
  #include "error_evaluator.h"
  #include "print_convergence_table.h"
  #include "assemble.h"
  #include "assemble_system_meshworker.h"
  #include "output_routines.h"
}
