#pragma once
# include <vector>
# include <iostream>

using namespace std;

int divmod (int a, int b); //modular division

struct coord{
	int x,y,z,l;
	coord(){}
	coord(const int&, const int&, const int&, const int&);
	coord(const int&, const coord&);
	int hash(const coord&);
	void operator=(const coord&);
	//dual lattice (cubes)
	int getFaceQubits(const coord&, const int&, const int&, const int&);
	int getFaceQubits(const coord&, const int&);
};

struct subcoord{
	int x,y,l;
	subcoord(){}
	subcoord(const int&, const int&, const int&);
//		subcoord(const int&);
	
	subcoord(const int& c, const subcoord&);
	int hash(const subcoord&);
	void operator=(const subcoord&);
};

//print coord
ostream& operator<<(ostream&, const coord&);