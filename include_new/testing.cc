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
		double test_fn(const double a);
};

double Child::tenst_fn(const double a)
{
	std::cout << "Value of a is " << a << std::endl;
}

int main(int argc,char **argv)
{
	Child child;
	Base *base  = &child;
	base.test_fn(10);

}