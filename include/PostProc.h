namespace PostProc
{
	using namespace dealii;
	using namespace std;
	using namespace Basics;

	// an abstract base class for post proc routines;
	// finally implemented in the Solver namespace
	template<int dim> class Base_PostProc
	{
		public:
			Base_PostProc();
			virtual void print_convergence_table(const string filename) = 0;
			virtual void error_evaluation(const Vector<double> &solution) = 0;
			virtual void output_solution_details(const Triangulation<dim> &triangulation,const string file_solution,
												const string file_exact,const string file_error)const = 0;
			virtual void evaluate_norms(const Vector<double> &solution) = 0;

	};

	template<int dim> Base_PostProc<dim>::Base_PostProc()
	{
	}


}
