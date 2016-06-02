
enum System_Type
{symmetric,un_symmetric};

enum Num_Flux
{Upwind, LLF};

enum Force_Type
{type1,type2};

enum Refinement
    {global,adaptive,adaptive_kelly};

enum BC_type
	{characteristic,odd};

struct numerical_data
{
	int p;
	int mapping_order;
	int refine_cycles;
	Refinement refinement;
};

struct nEqn_data
{
	unsigned int no_of_total_systems;
	unsigned int system_to_solve;

	std::vector<unsigned int> total_nEqn;
	std::vector<unsigned int> nBC;
	std::map<unsigned int,char> system_id_nEqn;			// maps the system_id with the number of equations
	
	Force_Type force_type;
	System_Type system_type;			
	BC_type bc_type;
};

struct physical_data
{
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
};