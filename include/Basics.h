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
		Base_Basics();
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
	};

	Base_Basics::Base_Basics()
	{

		// variables corresponding to exact solution for system-A
			tau = 0.10;
			zeta = 1.0;					//zeta
			chi = 1.0;
			theta0 = 2.0;
			theta1 = 1.0; 
			uW = 0.1;
			A0 = 0.0;
			A1 = 0.2;
			A2 = 0.1;

		// variables corresponding to exact solution for system-B
	/*		tau = 0.10;
			zeta = 1.0;					//zeta
			kappa = 0.0;
			chi = 1.0;
			theta0 = 2.0;
			theta1 = 1.0; 
			uW = 0.0;

			A0 = 0.0;
			A1 = 0.0;
			A2 = 0.1;*/


			 parameters =
			 " tau " + to_string(tau) +
			 " zeta " + to_string (zeta) +
			 " chi " + to_string(chi) +
			 " theta0 " + to_string(theta0) +
			 " theta1 " + to_string(theta1) +
			 " uW " + to_string(uW) +
			 " A0 " + to_string(A0) +
			 " A1 " + to_string(A1) +	
			 " A2 " + to_string(A2);

			 main_output_dir = "outputs_trial";
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

}
