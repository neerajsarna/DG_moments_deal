namespace PostProc
{
	using namespace dealii;
	using namespace std;
	using namespace Basics;

	template<int dim> class Base_PostProc
	{
		public:
			Base_PostProc();
			virtual void print_convergence_table(const string filename) = 0;
			virtual void error_evaluation(const Vector<double> solution) = 0;
			/*void output_solution_at_vertex(const Triangulation<dim> triangulation,const string filename)const = 0;
			void output_exact_solution_at_vertex(const Function<dim> exact_solution,const string filename)const = 0;
			void output_error_at_vertext(const string filename)const = 0;*/
	};

	template<int dim> Base_PostProc<dim>::Base_PostProc()
	{
	};


}