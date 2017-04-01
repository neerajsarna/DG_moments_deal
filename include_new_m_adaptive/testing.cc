#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>

using namespace std;


struct test
{	
	double a;
	void print();
};

void test::print()
{
	std::cout << "value of a is " << a << std::endl;
}

int main(int argc,char **argv)
{
	test neeraj;
	neeraj.a = 10;
	
	neeraj.print();

}