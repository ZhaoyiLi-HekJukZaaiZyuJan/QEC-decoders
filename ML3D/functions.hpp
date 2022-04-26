//===================================================================//
//===================  ML Assisted decoding (3D)=====================//
//=====================  functions.hpp      =========================//
//===================================================================//

# pragma once
# include <stdio.h>
# include <iostream>
# include <fstream>
# include <random>
# include <vector>
# include <cstdlib>
# include <string>
# include <algorithm>
# include <functional>
# include <tensorflow/c/c_api.h>
# include "cppflow/cppflow.h"
#include  "cppflow/ops.h"
#include  "cppflow/model.h"

using namespace std;

int divmod (int, int); //modular division
//===================================================================//

enum surfacetype {TORUS, PLANE};
istream& operator>> (istream&, surfacetype&);
string to_string(surfacetype&);

//===================================================================//
enum lossmodel {OFF_loss, INDEP_loss, toLoss};


//===================================================================//
enum noisemodel {OFF, INDEP, EM2, EM2_full, GATE, DEPOL1, DEPOL2, LOSS, TEST, MANUAL};
//INDEP: p1 = p, p2 = q
//DEPOL1: simulated time depol
//DEPOL2: Pseudo depol noise with primal-dual lattice identification
//LOSS: invoke loss model and loss decoder
//EM2:
//GATE:
//TEST:


//===================================================================//
//################################## vector printing ##############################
template <class T>
ostream& operator<<(ostream&, const vector<T>);

template <class T>
void printMatrix(const vector<T>);

//===================================================================//



