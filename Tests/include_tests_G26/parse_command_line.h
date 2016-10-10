class parse_command_line
{
	public:
		parse_command_line(int argc,char **argv);
		void read_input_filename(std::string &output_dir);

	private:
		std::list<std::string> args;
};

parse_command_line
::parse_command_line(int argc,char **argv)
{
	if (argc != 4)
	{
		std::cout << "command line input not in proper format" << std::endl;
		exit(1);
	}

	for ( int i = 2 ; i < argc ; i++)
		args.push_back(argv[i]);
}

void
parse_command_line
::read_input_filename(std::string &output_dir)
{
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
			output_dir = args.front();
			args.pop_front();
		}
	}

	printf("Using %s as input file",output_dir.c_str());
}