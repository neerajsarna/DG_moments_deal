struct nEqn_data
{
	unsigned int no_of_total_systems;
	std::vector<unsigned int> total_nEqn;
	std::map<unsigned int,char> system_id_nEqn;			// maps the system_id with the number of equations
};

enum System_Type
{symmetric,un_symmetric};

enum Num_Flux
{Upwind, LLF};