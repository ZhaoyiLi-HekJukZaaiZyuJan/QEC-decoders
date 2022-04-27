# include <iostream>
# include <random>
# include <vector>
# include <cstdlib>
# include <string>
# include <algorithm>
# include <functional>
using namespace std;
#include "functions.hpp"

#include "PerfectMatching.h"
#include "GEOM/GeomPerfectMatching.h"


int divmod (int a, int b) {return  (a%b+b) % b;} //modular division
istream& operator>> (istream& is, surfacetype& aSurf){
	string str;
	is >> str;
	if (str == "PLANE") aSurf = PLANE;
	else;
	return is;
}
string to_string(surfacetype& surf){
	if (surf == PLANE) return "PLANE";
	return "";
}

istream& operator>> (istream& is, lossmodel& N){
	string str;
	is >> str;
	if (str == "INDEP_loss") N = INDEP_loss;
	else;
	return is;
}
string to_string(lossmodel& N){
	if (N == INDEP_loss) return "INDEP_loss";
	return "";
}

istream& operator>> (istream& is, noisemodel& N){
	string str;
	is >> str;
	if (str == "INDEP") N = INDEP;
	else if (str == "EM2") N = EM2;
	else if (str == "GATE") N = GATE;
	else if (str == "INDEP_211") N = INDEP_211;
	else;
	return is;
}
string to_string(noisemodel& N){
	if (N == INDEP) return "INDEP";
	else if (N == EM2) return "EM2";
	else if (N == GATE) return "GATE";
	else if (N == INDEP_211) return "INDEP_211";
	return "";
}

//################################## vector printing ##############################
template <class T>
ostream& operator<<(ostream& os, const vector<T> vec) {
	os << "[";
	for (int i = 0; i < vec.size(); i++) {
		os << vec[i] << ", ";
	}
	os << "]";
	return os;
}

template <class T, class S>
ostream& operator<<(ostream& os, const pair<T,S> aPair) {
	os << "(" << aPair.first << "," << aPair.second << ")";
	return os;
}

template <class T>
void printMatrix(const vector<T> vec) {
	for (int i=0; i<vec.size(); i++) {
		cout << vec[i] << endl;
	}
}

coord::coord(const int& x, const int& y, const int& z, const int& l){
	this->x = x;
	this->y = y;
	this->z = z;
	this->l = l;
};
	
coord::coord(const int&c, const coord& S){
	int L = S.x;
	int M = S.y;
	int N = S.z;
	*this = {c % (M*L*N) % (M*L) % L, c % (M*L*N) % (M*L) / L, c % (M*L*N) / (M*L), c / (M*L*N)};
}
int coord::hash(const coord& S){
	int L = S.x;
	int M = S.y;
	int N = S.z;
	return l*L*M*N + z*L*M + y*L + x;
}
void coord::operator=(const coord& c) {
	x = c.x;
	y = c.y;
	z = c.z;
	l = c.l;
}


//dual lattice (cubes)
int coord::getFaceQubits(const coord& S, const int& face, const int& direction, const int& k){
	vector<int> helper {0,1,2,3};
	if (direction == 1) {
		iter_swap(helper.begin() + 0, helper.begin() + 3);
		iter_swap(helper.begin() + 1, helper.begin() + 2);
	}
	if (face == 0) { // x-y plane
		if (S.z % 1 == 0) {
			iter_swap(helper.begin() + 0, helper.begin() + 1);
			iter_swap(helper.begin() + 2, helper.begin() + 3);
		}
		if (k == helper[0]) {
			return coord(x, divmod(y+1, S.y), z, 0).hash(S);//0
		} else if (k == helper[1]) {
			return hash(S);//1
		} else if (k == helper[2]) {
			return coord(x, y, z, 1).hash(S);//2
		} else{
			return coord(divmod(x + 1, S.x), y, z, 1).hash(S);//3
		}
	} else if (face == 1) { // z-y plane
		if (S.x % 1 == 0) {
			iter_swap(helper.begin() + 0, helper.begin() + 1);
			iter_swap(helper.begin() + 2, helper.begin() + 3);
		}
		if (k == helper[0]) {
			return coord(x, y, z, 2).hash(S);//4m
		} else if (k == helper[1]) {
			return coord(divmod(S.x+1, S.x), y, z, 2).hash(S);//6
		} else if (k == helper[2]) {
			return coord(x, divmod(y+1, S.y), divmod(z+1, S.z), 0).hash(S);//8
		} else{
			return coord(x, divmod(y+1, S.y), z, 0).hash(S);//0
		}
	} else { //x-z plane
		if (S.y % 1 == 0) {
			iter_swap(helper.begin() + 0, helper.begin() + 1);
			iter_swap(helper.begin() + 2, helper.begin() + 3);
		}
		if (k == helper[0]) {
			return coord(x, y, z, 1).hash(S);//2
		} else if (k == helper[1]) {
			return coord(x, y, divmod(z+1, S.z), 1).hash(S);//10
		} else if (k == helper[2]) {
			return coord(x, divmod(y+1, S.y), z, 2).hash(S);//5
		} else{
			return coord(x, y, z, 2).hash(S);//4
		}
	}
}

//print coord
ostream& operator<<(ostream& os, const coord& c) {
	os << "(" << c.x << "," << c.y << "," << c.z << "," << c.l << ")";
	return os;
}

//coord comparison
bool operator<(const coord& c1, const coord& c2) {
	if (c1.x < c2.x) {
		return true;
	} else if (c1.x > c2.x) {
		return false;
	} else {
		if (c1.y < c2.y){
			return true;
		} else if (c1.y > c2.y){
			return false;
		} else {
			if (c1.z < c2.z) {
				return true;
			} else {
				return false;
			}
		}
	}
}

bool operator==(const coord& c1, const coord& c2) {
	if (c1.x == c2.x && c1.y == c2.y && c1.z == c2.z && c1.l == c2.l) return true;
	else return false;
}

vertex::vertex(const int& c, const coord& S){
	this->c = c;
	coord C(c,S);
	partial ={
		coord(divmod(C.x-1, S.x), C.y, C.z, 0).hash(S), coord(C.x, C.y, C.z, 0).hash(S),
		coord(C.x, divmod(C.y-1, S.y), C.z, 1).hash(S), coord(C.x, C.y, C.z, 1).hash(S),
		coord(C.x, C.y, divmod(C.z-1, S.z), 2).hash(S), coord(C.x, C.y, C.z, 2).hash(S),
	};
	//each vertex contains 6 physical_qubits, left (Mx-), right(Mx+), up (Ly-), down (Ly+), less (Nz-), more (Nz+)
};

int getTaxicabDistance(const coord& S, const int& c1, const int& c2){  // compute taxicab distance between two coords
	coord C1(c1,S);
	coord C2(c2,S);
	return abs(C1.x - C2.x) + abs(C1.y - C2.y) + abs(C1.z - C2.z);
}
//
coord getTaxicabDisplacement(const coord& S, const int& c1, const int& c2){  //
	coord C1(c1,S);
	coord C2(c2,S);
	return {C2.x - C1.x, C2.y - C1.y, C2.z - C1.z, 0};
}

int getTaxicabDistance(const coord& S, const vector<int>& chunk1, const vector<int>& chunk2){  // compute taxicab distance between two coords
	int min_dist = INT_MAX;
	for (int i = 0; i < chunk1.size(); i++) {
		for (int j = 0; j < chunk2.size(); j++) {
			int dist = getTaxicabDistance(S, chunk1[i], chunk2[j]);
			min_dist = (dist < min_dist) ? dist : min_dist;
		}
	}
	return min_dist;
}
