# pragma once
# include <iostream>
# include <fstream>
# include <random>
# include <vector>
# include <cstdlib>
# include <string>
# include <algorithm>
# include <functional>

using namespace std;

int divmod (int, int); //modular division
//===================================================================//

enum subsurfacetype {TORUS, PLANE};
istream& operator>> (istream&, subsurfacetype&);

//string to_string(subsurfacetype&);

//===================================================================//
//################################## vector printing ##############################
template <class T>
ostream& operator<<(ostream&, const vector<T>);

template <class T>
void printMatrix(const vector<T>);

//===================================================================//



