//===================================================================//
//===================  ML Assisted decoding (3D)=====================//
//=====================  functions.hpp      =========================//
//===================================================================//


#include "functions.hpp"
#include <map>
 
using namespace std;
int divmod (int a, int b) {return  (a%b+b) % b;} //modular division
//===================================================================//

istream& operator>> (istream& is, surfacetype& aSurf){
	string str;
	is >> str;
	if (str == "PLANE") aSurf = PLANE;
	if (str == "TORUS") aSurf = TORUS;
	else;
	return is;
}
string to_string(surfacetype& surf){
	if (surf == PLANE) return "PLANE";
	if (surf == TORUS) return "TORUS";
	else return "";
}

istream& operator>> (istream& is, noisemodel& N){
	string str;
	is >> str;
	if (str == "INDEP") N = INDEP;
	else if (str == "EM2") N = EM2;
	else if (str == "EM2_full") N = EM2_full;
	else if (str == "GATE") N = GATE;
	else if (str == "DEPOL1") N = DEPOL1;
	else if (str == "DEPOL2") N = DEPOL2;
	else if (str == "LOSS") N = LOSS;
	else if (str == "MANUAL") N = MANUAL;
	else if (str == "TEST") N = TEST;
	else if (str == "OFF") N = OFF;
	else;
	return is;
}
string to_string(noisemodel& N){
	if (N == INDEP) return "INDEP";
	else if (N == EM2) return "EM2";
	else if (N == DEPOL1) return "DEPOL1";
	else if (N == DEPOL2) return "DEPOL2";
	else if (N == TEST) return "TEST";
	else if (N == LOSS) return "LOSS";
	else if (N == EM2_full) return "EM2_full";
	else if (N == GATE) return "GATE";
	else if (N == MANUAL) return "MANUAL";
	else if (N == OFF) return "OFF";
	else return "";
}

typedef function<pair<double, double>(double, double)> noisemodelfunc;

map<noisemodel, noisemodelfunc> NOISEMODELMAP {
	{EM2, [](const double& p,const double& q)->pair<double,double>{return {p*32/15,p*4/15};}}, //to first order
	{INDEP, [](const double& p,const double& q)->pair<double,double>{return {p, q};}},
	{TEST, [](const double& p,const double& q)->pair<double,double>{return {16/15*p, 4/15*p};}},
	{EM2_full, [](const double& p,const double& q)->pair<double,double>{return {p*62/15 - pow(p,2)*1096/75 + pow(p,3)*32224/1125 - pow(p,4)*567296 /16875 + pow(p,5)*1196032/50625 - pow(p,6)*4194304/455625 + pow(p,7)*2097152/1366875, 1/2-sqrt(1/4-p*4/15)};}}
};

istream& operator>> (istream& is, lossmodel& L){
	string str;
	is >> str;
	if (str == "INDEP_loss") L = INDEP_loss;
	else if (str == "OFF_loss") L = OFF_loss;
	else;
	return is;
}
string to_string(lossmodel& L){
	if (L == INDEP_loss) return "INDEP_loss";
	if (L == OFF_loss) return "OFF_loss";
	return "";
}

//===================================================================//
typedef function<pair<double, double>(double)> lossmodelfunc;


map<lossmodel, lossmodelfunc> LOSSMODELMAP {
	{INDEP_loss, [](const double& p)->pair<double,double>{return {p,0};}}
};

//===================================================================//
//################################## vector printing ##############################
template <class T>
ostream& operator<<(ostream& os, const vector<T> vec) {
	for (int i = 0; i < vec.size(); i++) {
		os << vec[i] << "|";
	}
	return os;
}

template <class T>
void printMatrix(const vector<T> vec) {
	for (int i=0; i<vec.size(); i++) {
		cout << vec[i] << endl;
	}
}


//===================================================================//
// ./simulate -s TORUS --pmin 0 -N TEST --pmax 0.007 --Np 10 --Nq 1 -n 10000 --Lmin 3 --Lmax 17 -v 1 -d "/scratch/users/ladmon/" -m 'model,L=3(7),layer=5x512,epochs=1000,p=' 
