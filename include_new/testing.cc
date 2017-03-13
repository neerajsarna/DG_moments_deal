#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>

using namespace std;

class Base
{
	public:
		Base(){};
		virtual void print_variable() = 0;
};

class Child:public Base
{
	public:	
		Child(){};
		virtual void print_variable() {std::cout << "Neeraj " << std::endl;};

};

class Base_Testing
{

	public:
		Base_Testing(){};
		Base *base;
		
};

class Child_Testing:public Base_Testing
{
	public:
		Child child;
		int a;
		Child_Testing(int a):
		a(a)
		{this->base = &child;};
		
};

int main(int argc,char **argv)
{
	std::vector<Child_Testing> test;

	for (int i = 0 ; i < 4 ; i++)
		test.push_back(Child_Testing(i));

	test[0].base->print_variable();

}