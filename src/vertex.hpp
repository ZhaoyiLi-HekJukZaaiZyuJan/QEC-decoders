#pragma once

#include "coord.hpp"
#include <vector>
using namespace std;

struct vertex {
	vertex(){}
	vertex(const int&, const coord&);
	int c;
	vector<int> partial;
};