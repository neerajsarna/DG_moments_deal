#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <map>

using namespace std;


void allocate_map(	std::vector<std::map<std::string,unsigned int> > &variable_map)
{
	std::map<std::string,unsigned int> map_1;

	map_1["neeraj"] = 0;
	variable_map.push_back(map_1);	
	std::cout << "size " << variable_map.size() << std::endl;
}

int main(int argc,char **argv)
{
	std::vector<std::map<std::string,unsigned int> > variable_map;

	allocate_map(variable_map);

	std::cout << variable_map[0].find("neeraj")->second << std::endl;
}