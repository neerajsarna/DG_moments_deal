#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

using namespace std;
class Base
{
	public:
		Base();

	struct test
	{
		double a;
		double b;
	};
};

class Child:public Base
{
	public:
		Child();
		test omega;
};

int main(int argc,char **argv)
{
	Child child;
	child.omega.a = 10.0;

}