#pragma once

# include <algorithm>
# include <cstdlib>
# include <fstream>
# include <functional>
# include <iostream>
# include <map>
# include <random>
# include <stdio.h>
# include <string>
# include <vector>
# include "coord.hpp"
# include "types.hpp"

using namespace std;
//===================================================================//
istream& operator>> (istream&, surfacetype&);
string To_string(surfacetype&);
ostream& operator<< (ostream&, surfacetype&);

//===================================================================//
istream& operator>> (istream&, subsurfacetype&);
string To_string(subsurfacetype&);
ostream& operator<< (ostream&, subsurfacetype&);

//===================================================================//
istream& operator>> (istream&, concattype&);
string To_string(concattype&);
ostream& operator<< (ostream&, concattype&);

//===================================================================//

istream& operator>> (istream&, lossmodel&);
string To_string(lossmodel&);
ostream& operator<< (ostream&, lossmodel&);

typedef function<pair<double, double>(double,double)> lossmodelfunc;

extern std::map<lossmodel, lossmodelfunc> LOSSMODELMAP;

//===================================================================//
//INDEP: p1 = p, p2 = q
//DEPOL1: simulated time depol
//DEPOL2: Pseudo depol noise with primal-dual lattice identification
//LOSS: invoke loss model and loss decoder
//EM2p0: Raussendorf model with p_M = p_S = p_P = 0; p_th = 0.96% 
//EM2:
//GATE:
// TEST:

istream& operator>> (istream&, noisemodel&);
string To_string(noisemodel&);
ostream& operator<< (ostream&, noisemodel&);

typedef function<pair<double, double>(double,double)> noisemodelfunc;

extern int myint;
extern std::map<noisemodel, noisemodelfunc> NOISEMODELMAP;


//################################## vector printing ##############################
template <class T>
ostream& operator<<(ostream&, const vector<T>);

template <class T, class S>
ostream& operator<<(ostream&, const pair<T,S>);

template <class T>
void printMatrix(const vector<T>);

#include "functions.tpp"

//======================== 3D Distance =========================//

coord getTaxicabDisplacement(const coord&, const int&, const int&, const surfacetype&);

int getTaxicabDistance(const coord&, const coord&, const coord&, const surfacetype&);

int getTaxicabDistanceChunks(const coord&, const vector<int>&, const vector<int>&, const surfacetype&);

//======================== 2D Distance =========================//

subcoord getTaxicabDisplacement(const subcoord&, const int&, const int&, const subsurfacetype&);

int getTaxicabDistance(const subcoord&, const int&, const int&, const subsurfacetype&);

void testDecoding(const int L, const int M, const double p, subsurfacetype);

void loopDecoding(const int L_range, const int trials, const string, subsurfacetype, bool);

//======================== Coordinate comparison =========================//
bool operator<(const coord&, const coord&);

bool operator==(const coord&, const coord&);

bool operator<(const subcoord&, const subcoord&);

bool operator==(const subcoord&, const subcoord&);