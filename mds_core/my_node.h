#ifndef MY_NODE
#define MY_NODE

#include <iostream>

using namespace std;

class my_node
{
public:
	my_node()
	{
		index = 0;
		value = 0;
		visit = false;
	}
	my_node(int index, double value, bool visit)
	{
		this->index = index;
		this->value = value;
		this->visit = visit;
	}
	int index;
	double value;
	bool visit;
};

struct cmp
{
	bool operator()(const my_node& n1, const my_node& n2)
	{
	return n1.value > n2.value;
	}
};

#endif
