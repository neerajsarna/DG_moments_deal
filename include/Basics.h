// the following function generates the directories where the output from the program will be stored
namespace Basics
{

	struct stat st = { .st_dev = 0 };
	using namespace std;
	using namespace dealii;
	class Base_Basics
	{
		// physical constants
	public:
		Base_Basics(const physical_data &physical_constants,
					const string &output_dir);
		double tau;
		double zeta;
		double chi;
		double theta0;
		double theta1;
		double uW;
		double A0;
		double A1;
		double A2;
		double kappa;

		string parameters;
		string main_output_dir;
		vector<string> sub_directory_names;
		string generate_filename_to_write(const string system_dir,string filename);
		void Sparse_matrix_dot_Vector(const  system_matrix matrix_info,
			Vector<double> &x) const;

		// the following function prints the matrix given to it 
		void print_dealii_matrix(const FullMatrix<double> &matrix,
					string matrix_name);

		void print_eigen_sparse(Sparse_matrix &sparse_matrix, 
								string matrix_name);

		void print_eigen_full(Full_matrix &full_matrix,
							  string matrix_name);

		void print_dealii_sparse(const TrilinosWrappers::SparseMatrix &matrix);
	};

	Base_Basics::Base_Basics(const physical_data &physical_constants,
						     const string &output_dir)
	{

		tau = physical_constants.tau;
		zeta = physical_constants.zeta;
		chi = physical_constants.chi;
		theta0 = physical_constants.theta0;
		theta1 = physical_constants.theta1;
		uW = physical_constants.uW;
		A0 = physical_constants.A0;
		A1 = physical_constants.A1;
		A2 = physical_constants.A2;
		kappa = physical_constants.kappa;		
		
		main_output_dir = output_dir;
		
		unsigned int num_outputs = 5;
		sub_directory_names.resize(num_outputs);

		for (unsigned int i = 0 ; i < num_outputs ; i++)
		{
			if (i == 0)
				sub_directory_names[0] = main_output_dir + "/grids";

			if (i == 1)
				sub_directory_names[1] = main_output_dir + "/solution";

			if (i == 2)
				sub_directory_names[2] = main_output_dir + "/convergence_tables";

			if (i == 3)
				sub_directory_names[3] = main_output_dir + "/sparsity_patterns";

			if ( i == 4)
				sub_directory_names[4] = main_output_dir + "/matrix_read";			
		}


		if (stat(main_output_dir.c_str(),&st) == -1)
			mkdir(main_output_dir.c_str(),0777);

		for (unsigned int i = 0 ; i < num_outputs ; i ++)
			if (stat(sub_directory_names[i].c_str(),&st) == -1)
				mkdir(sub_directory_names[i].c_str(),0777);

		}

		string Base_Basics
		::generate_filename_to_write(const string system_dir,string filename)
		{
			string file_to_write = sub_directory_names[4];
			const unsigned int position_to_erase = filename.find(system_dir);
			const unsigned int len_to_erase = system_dir.length();

			file_to_write += "/" + filename.erase(position_to_erase,len_to_erase);

			return(file_to_write);
		}

		void Base_Basics
		::Sparse_matrix_dot_Vector(const  system_matrix matrix_info,
			Vector<double> &x) const
		{
			Assert(x.size() != 0 || 
				matrix_info.matrix.rows() != 0 ||
				matrix_info.matrix.cols(),ExcNotInitialized());

			Assert(matrix_info.matrix.IsRowMajor,ExcInternalError());

			Vector<double> result(x.size());

			for (unsigned int m = 0 ; m < matrix_info.matrix.outerSize(); m++)
			{
				result(m) = 0;
				for (Sparse_matrix::InnerIterator n(matrix_info.matrix,m); n ; ++n)
					result(m) += n.value() * x(n.col());
			}

			for (unsigned int i = 0 ; i < result.size() ; i ++)
				x(i) = result(i);
		}

	 // print the matrix to standard input output
	  void  Base_Basics::print_dealii_matrix(const FullMatrix<double> &matrix,
					   string matrix_name) 
	  {
		const unsigned int n_rows = matrix.m();
		const unsigned int n_cols = matrix.n();

		printf("%s \n",matrix_name.c_str());

		for (unsigned int i = 0 ; i < n_rows ; i ++)
		{
			for (unsigned int j = 0 ; j < n_cols ; j ++)
				cout << matrix(i,j) << " " ;

			cout << endl;
		}

	  }

	  void Base_Basics::print_eigen_sparse(Sparse_matrix &sparse_matrix,
	  									   string matrix_name)
	  {
	  	const unsigned int n_rows = sparse_matrix.rows();
	  	const unsigned int n_cols = sparse_matrix.cols();

	  	printf("%s\n",matrix_name.c_str());
	  	for (unsigned int i = 0 ;  i < n_rows ; i ++)
	  	{
	  		for (unsigned int j = 0 ; j < n_cols ; j ++)
	  			cout << " " << sparse_matrix.coeffRef(i,j) ; 

	  		cout << endl;
	  	}

	  }

	  void Base_Basics::print_eigen_full(Full_matrix &full_matrix,
	  									 string matrix_name)
	  {
	  	const unsigned int n_rows = full_matrix.rows();
	  	const unsigned int n_cols = full_matrix.cols();

	  	printf("%s\n",matrix_name.c_str());
	  	for (unsigned int i = 0 ; i < n_rows; i ++)
	  	{
	  		for (unsigned int j = 0 ; j < n_cols ; j++)
	  			cout << " " << full_matrix(i,j) ;

	  		cout << endl;
	  	}
	  }

	  void Base_Basics::print_dealii_sparse(const TrilinosWrappers::SparseMatrix &matrix)
	  {
	  	FILE *fp;
	  	fp = fopen("printed_matrices/global_matrix","w+");

	  	const unsigned int num_non_zero = matrix.n_nonzero_elements();
	  	typename TrilinosWrappers::SparseMatrix::const_iterator it = matrix.begin();
	  	const typename TrilinosWrappers::SparseMatrix::const_iterator it_end = matrix.end();

	  	for (; it != it_end; it++)
	  		fprintf(fp,"%u %u %lf \n",it->row(),it->column(),it->value());
	  	

	  	fclose(fp);


	  }

	}
