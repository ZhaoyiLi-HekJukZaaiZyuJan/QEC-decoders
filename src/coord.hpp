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

//print coord
ostream& operator<<(ostream&, const coord&);

//coord comparison
bool operator<(const coord&, const coord&);

bool operator==(const coord&, const coord&);