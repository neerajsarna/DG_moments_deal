#include <typeinfo>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <list>
#include <algorithm>

int main(int argc,char **argv)
{
	FILE *fp;


	for (unsigned int i = 5 ; i <= 50; i++)
	{
		std::string filename = "input_N" + std::to_string(i) + ".in";
		int nBC = i/2;

	fp = fopen(filename.c_str(),"w+");

	fprintf(fp, "%s\n","subsection System Properties");
	fprintf(fp , "\t %s %d\n","set total_systems = ",1 );
	fprintf(fp , "\t %s %d\n","set Ntensors0 = ",i );
	fprintf(fp , "\t %s %d\n","set system_id0 = ",i );
	fprintf(fp , "\t %s %d\n","set num_equations0 = ",i );
	fprintf(fp , "\t %s %d\n","set nBC0 = ",nBC );
	fprintf(fp , "%s","end");
	fclose(fp);
	}



	



}