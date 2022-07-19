#include "subfunctions.hpp"
using namespace std;

int divmod (int a, int b) {return  (a%b+b) % b;} //modular division
//===================================================================//

istream& operator>> (istream& is, subsurfacetype& aSurf){
	string str;
	is >> str;
	if (str == "subTORUS") aSurf = subTORUS;
	if (str == "subPLANE") aSurf = subPLANE;
	else;
	return is;
}
//===================================================================//
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

template <class T>
void printMatrix(const vector<T> vec) {
	for (int i=0; i<vec.size(); i++) {
		cout << vec[i] << endl;
	}
}

//===================================================================//
