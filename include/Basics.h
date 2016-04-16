// the following function generates the directories where the output from the program will be stored
namespace Basics
{

	struct stat st = {0};
	using namespace std;
	class Base_Basics
	{
		// physical constants
	public:
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
		unsigned int nEqn;

		string parameters;
		string main_output_dir;
		vector<string> sub_directory_names;
	};

	class constants_SystemA:public Base_Basics
	{
		public:
			constants_SystemA();

	};

	constants_SystemA::constants_SystemA()
	{
				tau = 0.10;
				zeta = 1.0;					//zeta
				chi = 1.0;
				theta0 = 2.0;
				theta1 = 1.0; 
				uW = 0.1;
				A0 = 0.0;
				A1 = 0.2;
				A2 = 0.1;
				nEqn = 6;

			 parameters =
			 " nEqn " + to_string(nEqn) +
			 " tau " + to_string(tau) +
			 " zeta " + to_string (zeta) +
			 " chi " + to_string(chi) +
			 " theta0 " + to_string(theta0) +
			 " theta1 " + to_string(theta1) +
			 " uW " + to_string(uW) +
			 " A0 " + to_string(A0) +
			 " A1 " + to_string(A1) +	
			 " A2 " + to_string(A2);

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


	class constants_SystemB:public Base_Basics
	{
		public:
			constants_SystemB();

	};


	constants_SystemB::constants_SystemB()
	{
		tau = 0.1;
		zeta = 1.0;
		kappa = 2.0;
		chi = 1.0;
  		theta0 = 2.0; // angular: 0.0, radial: 2.0
  		theta1 = 1.0; // angular: 0.0, radial: 1.0
  		uW = 0.0; // angular: 5.0, radial: 0.0
  
  		A0 = 0.0; // angular: 0.0, radial: 0.0
  		A1 = 0.0; // angular: -0.4, radial: 0.0
  		A2 = 0.1; // angular: 0.0, radial: 0.1
		nEqn = 10;  

		parameters =
		" nEqn " + to_string(nEqn) +
		" tau " + to_string(tau) +
		" kappa " + to_string(kappa) + 
		" zeta " + to_string (zeta) +
		" chi " + to_string(chi) +
		" theta0 " + to_string(theta0) +
		" theta1 " + to_string(theta1) +
		" uW " + to_string(uW) +
		" A0 " + to_string(A0) +
		" A1 " + to_string(A1) +	
		" A2 " + to_string(A2);

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



	};
}