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

struct subvertex {
	subvertex(){}
	subvertex(const int& hash, const subcoord&);
	int x;
	int y;
	vector<int> partial;
};
