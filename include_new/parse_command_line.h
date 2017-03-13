class parse_command_line
{
	public:
		parse_command_line(int argc,char **argv);
		void read_input_filename(std::vector<std::string> &output_dir);

	private:
		std::list<std::string> args;
};

parse_command_line
::parse_command_line(int argc,char **argv)
{
	// first arguments is the name of the executable
	// second is the number of threads to be used for meshworker
	// thirst is -p
	// fourth is the name of the input file for the system
	// again -p
	// fifth is the name of the input file for the numerical experiment
	if (argc != 6)
	{
		std::cout << "command line input not in proper format" << std::endl;
		exit(1);
	}

	// collect all the input arguments in a string for ease of manipulations
	for ( int i = 2 ; i < argc ; i++)
		args.push_back(argv[i]);
}

void
parse_command_line
::read_input_filename(std::vector<std::string> &output_dir)
{
	unsigned int read_filename = 0;
	while(args.size())
	{


		if (args.front() == std::string("-p"))
		{
			if (args.size() == 1)
				{
					std::cout << "-p must be followed by input file name" << std::endl;
					exit(1);
				}

			args.pop_front();
			output_dir[read_filename] = args.front();
			args.pop_front();
			read_filename++;
		}
	}
	AssertDimension(read_filename,2);

	printf("Using %s and %s  as input file \n",output_dir[0].c_str(),output_dir[1].c_str());
}
