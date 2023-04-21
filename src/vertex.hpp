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

struct subvertex_x{
	subvertex_x(){};
	subvertex_x(const int&, const subcoord&);
	
	int x;
	int y;
	vector<int> partial;
};

struct subvertex_z{
	subvertex_z(){};
	subvertex_z(const int&, const subcoord&);
	
	int x;
	int y;
	vector<int> partial;
};
