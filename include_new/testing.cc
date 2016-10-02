#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

using namespace std;

class Base
{
	public:
		enum solver_type
		{Trilinos,Pardiso};

		Base(solver_type &solver)
		{std::cout << "solver value is " << solver << std::endl;}
};

int main(int argc,char **argv)
{
	Base base(enum Base::solver_type Trilinos);
}