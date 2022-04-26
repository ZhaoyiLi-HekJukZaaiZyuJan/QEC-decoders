//===================================================================//
//===================  ML Assisted decoding     =====================//
//=====================  functions.cpp      =========================//
//===================================================================//

#include "functions.hpp"
using namespace std;

int divmod (int a, int b) {return  (a%b+b) % b;} //modular division
//===================================================================//

istream& operator>> (istream& is, subsurfacetype& aSurf){
	string str;
	is >> str;
	if (str == "TORUS") aSurf = TORUS;
	if (str == "PLANE") aSurf = PLANE;
	else;
	return is;
}
//===================================================================//
//################################## vector printing ##############################
template <class T>
ostream& operator<<(ostream& os, const vector<T> vec) {
	for (int i = 0; i < vec.size(); i++) {
		os << vec[i] << ", ";
	}
	return os;
}

template <class T>
void printMatrix(const vector<T> vec) {
	for (int i=0; i<vec.size(); i++) {
		cout << vec[i] << endl;
	}
}

//========================== Noise Model =============================//

istream& operator>> (istream& is, noisemodel& N){
	string str;
	is >> str;
	if (str == "INDEP") N = INDEP;
	else if (str == "DEPOL") N = DEPOL;
	else;
	return is;
}
string tostring(const noisemodel& N){
	if (N == INDEP) return "INDEP";
	else if (N == DEPOL) return "DEPOL";
	return "";
}