// the following function generates the directories where the output from the program will be stored
namespace Basics
{

	struct stat st = {0};
	using namespace std;
	class Base_Basics
	{
		// physical constants
	public:
		Base_Basics();
		string main_output_dir;
		vector<string> sub_directory_names;
		string generate_filename_to_write(const string folder_name,const string filename) const;
		double tau;
		double zeta;					//zeta
		double chi;
		double theta0;
		double theta1; 
		double kappa;
		double uW;
		double A0;
		double A1;
		double A2;


	};

	Base_Basics::Base_Basics()
	{
				kappa = 0;
				tau = 0.10;
				zeta = 1.0;					//zeta
				chi = 1.0;
				theta0 = 2.0;
				theta1 = 1.0; 
				uW = 0.0;
				A0 = 0.0;
				A1 = 0.0;
				A2 = 0.1;

			 main_output_dir = "outputs";
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

	string Base_Basics::generate_filename_to_write(const string folder_name,const string filename) const
	{

			string filename_out;
			filename_out = sub_directory_names[4] + "/" + filename.substr(folder_name.size());

			return filename_out;
	}


}
