# include <iostream>
# include <random>
# include <vector>
# include <cstdlib>
# include <string>
# include <algorithm>
# include <functional>

# include "functions.hpp"

# include "libs/blossom5-v2.05.src/PerfectMatching.h"
# include "libs/blossom5-v2.05.src/GEOM/GeomPerfectMatching.h"

using namespace std;
//===================================================================//

istream& operator>> (istream& is, surfacetype& aSurf){
	string str;
	is >> str;
	if (str == "PLANE") aSurf = PLANE;
	else if (str == "TORUS") aSurf = TORUS;
	else;
	return is;
}
string To_string(surfacetype& surf){
	if (surf == PLANE) return "PLANE";
	else if (surf == TORUS) return "TORUS";
	return "";
}

ostream& operator<< (ostream& os,surfacetype & surf){
	os << To_string(surf);
	return os;
}

//===================================================================//

istream& operator>> (istream& is, concattype& concat){
	string str;
	is >> str;
	if (str == "concat311a") concat = concat311a;
	else if (str == "concat311b") concat = concat311b;
	else;
	return is;
}
string To_string(concattype& concat){
	if (concat == concat311a) return "concat311a";
	else if (concat == concat311b) return "concat311b";
	else;
	return "";
}

ostream& operator<< (ostream& os,concattype & concat){
	os << To_string(concat);
	return os;
}

//===================================================================//

istream& operator>> (istream& is, noisemodel& N){
	string str;
	is >> str;
	if (str == "INDEP") N = INDEP;
	else if (str == "INDEP_211") N = INDEP_211;
	else if (str == "EM2") N = EM2;
	else if (str == "EM2p0") N = EM2p0;
	else if (str == "EM2_full") N = EM2_full;
	else if (str == "GATE") N = GATE;
	else if (str == "GATE_full") N = GATE_full;
	else if (str == "GATE_biased") N = GATE_biased;
	else if (str == "DEPOL1") N = DEPOL1;
	else if (str == "DEPOL2") N = DEPOL2;
	else if (str == "LOSS") N = LOSS;
	else if (str == "MANUAL") N = MANUAL;
	else if (str == "TEST") N = TEST;
	else if (str == "OFF") N = OFF;
	else;
	return is;
}
string To_string(noisemodel& N){
	if (N == INDEP) return "INDEP";
	else if (N == INDEP_211) return "INDEP_211";
	else if (N == EM2) return "EM2";
	else if (N == EM2p0) return "EM2p0";
	else if (N == EM2_full) return "EM2_full";
	else if (N == GATE) return "GATE";
	else if (N == GATE_full) return "GATE_full";
	else if (N == GATE_biased) return "GATE_biased";
	else if (N == DEPOL1) return "DEPOL1";
	else if (N == DEPOL2) return "DEPOL2";
	else if (N == TEST) return "TEST";
	else if (N == LOSS) return "LOSS";
	else if (N == MANUAL) return "MANUAL";
	else if (N == OFF) return "OFF";
	else return "";
}

ostream& operator<< (ostream& os, noisemodel& N){
	os << To_string(N);
	return os;
}

map<lossmodel, lossmodelfunc> LOSSMODELMAP = {
	{INDEP_loss, [](const double& p,const double& q)->pair<double,double>{return {p,0};}}
};

//===================================================================//

istream& operator>> (istream& is, lossmodel& L){
	string str;
	is >> str;
	if (str == "INDEP_loss") L = INDEP_loss;
	else if (str == "OFF_loss") L = OFF_loss;
	else if (str == "toLoss") L = toLoss;
	else;
	return is;
}

string To_string(lossmodel& L){
	if (L == INDEP_loss) return "INDEP_loss";
	if (L == OFF_loss) return "OFF_loss";
	if (L == toLoss) return "toLoss";
	return "";
}


ostream& operator<< (ostream& os, lossmodel& L){
	os << To_string(L);
	return os;
}

map<noisemodel, noisemodelfunc> NOISEMODELMAP = {
	{EM2p0, [](const double& p,const double& q)->pair<double,double>{return {p*32/15, p*4/15};}}, //to first order
	{EM2, [](const double& p,const double& q)->pair<double,double>{return {p*32/15 + 2*p, p*4/15};}}, //to first order
	{INDEP, [](const double& p,const double& q)->pair<double,double>{return {p, q};}},
	{TEST, [](const double& p,const double& q)->pair<double,double>{return {16/15*p, 4/15*p};}},
	{EM2_full, [](const double& p,const double& q)->pair<double,double>{return {p*62/15 - pow(p,2)*1096/75 + pow(p,3)*32224/1125 - pow(p,4)*567296 /16875 + pow(p,5)*1196032/50625 - pow(p,6)*4194304/455625 + pow(p,7)*2097152/1366875, 1/2-sqrt(1/4-p*4/15)};}}
};


//################################## vector printing ##############################

typedef function<pair<double, double>(double, double)> noisemodelfunc;

//################################## other functions ##############################

//
coord getTaxicabDisplacement(const coord& S, const int& c1, const int& c2, const surfacetype& surf){  //
	coord C1(c1,S);
	coord C2(c2,S);
	int L = S.x, M = S.y, x1 = C1.x, x2 = C2.x, y1 = C1.y, y2 = C2.y, z1 = C1.z, z2 = C2.z;
	if (surf == TORUS) { //TORUS geometry
			int x_dist, y_dist, x_relative_pos, y_relative_pos;
			x_dist = min(M-abs(x1 - x2), abs(x1 - x2));
			y_dist = min(L-abs(y1 - y2), abs(y1 - y2));

			//replace this part with displacement function
			if (x_dist == 0) { //two vertices on same x
				x_relative_pos = 0;

			} else if ((x_dist == abs(x1 - x2) && x1 < x2) ||
					(x_dist == M - abs(x1 - x2) && x1 > x2)){
				x_relative_pos = 1;
			} else {
				x_relative_pos = -1;
			}

			if (y_dist == 0) { //two vertices on same y
				y_relative_pos = 0;
			} else if ((y_dist == abs(y1 - y2) && y1 < y2) ||
					(y_dist == L - abs(y1 - y2) && y1 > y2)){
				y_relative_pos = 1;

			} else {
				y_relative_pos = -1;
			}
			return {x_relative_pos * x_dist, y_relative_pos * y_dist, z2 - z1, 0};
	} else if (surf == PLANE){
		return {x2 - x1, y2 - y1, z2 - z1, 0};
	} else {
		return {0,0,0,0};
	}
}

int getTaxicabDistance(const coord& S, const int& c1, const int& c2, const surfacetype& surf){  // compute taxicab distance between two coords
	coord C1(c1,S);
	coord C2(c2,S);
	int L = S.x, M = S.y;
	if (surf == TORUS) {
		return min(abs(C1.x - C2.x), M-abs(C1.x - C2.x)) + min(abs(C1.y - C2.y), L-abs(C1.y - C2.y)) + abs(C1.z - C2.z);
	} else if (surf == PLANE){
		coord C1(c1,S);
		coord C2(c2,S);
		return abs(C1.x - C2.x) + abs(C1.y - C2.y) + abs(C1.z - C2.z);
	} else return INT_MAX;
}

int getTaxicabDistance(const coord& S, const coord& C1, const coord& C2, const surfacetype& surf){  // compute taxicab distance between two coords
	int L = S.x, M = S.y;
	if (surf == TORUS) {
		return min(abs(C1.x - C2.x), M-abs(C1.x - C2.x)) + min(abs(C1.y - C2.y), L-abs(C1.y - C2.y)) + abs(C1.z - C2.z);
	} else if (surf == PLANE){
		return abs(C1.x - C2.x) + abs(C1.y - C2.y) + abs(C1.z - C2.z);
	} else return INT_MAX; 
}

int getTaxicabDistanceChunks(const coord& S, const vector<int>& chunk1, const vector<int>& chunk2, const surfacetype& surf){  // compute taxicab distance between two coords
	int min_dist = INT_MAX;
	for (int i = 0; i < chunk1.size(); i++) {
		for (int j = 0; j < chunk2.size(); j++) {
			int dist = getTaxicabDistance(S, chunk1[i], chunk2[j], surf);
			min_dist = (dist < min_dist) ? dist : min_dist;
		}
	}
	return min_dist;
}
