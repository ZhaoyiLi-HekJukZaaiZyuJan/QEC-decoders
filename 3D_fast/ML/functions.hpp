//===================================================================//
//===================  ML Assisted decoding     =====================//
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
enum subsurfacetype {TORUS, PLANE};
istream& operator>> (istream&, subsurfacetype&);
//string to_string(subsurfacetype&);

//========================== Noise Model =============================//
enum noisemodel {INDEP, DEPOL};
istream& operator>> (istream&, noisemodel&);
string to_string(noisemodel&);

//===================================================================//
//################################## vector printing ##############################
template <class T>
ostream& operator<<(ostream&, const vector<T>);

template <class T>
void printMatrix(const vector<T>);

//===================================================================//



