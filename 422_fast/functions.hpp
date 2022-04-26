# include <iostream>
# include <random>
# include <vector>
# include <cstdlib>
# include <string>
# include <algorithm>
# include <functional>

using namespace std;

int divmod (int, int); //modular division
enum surfacetype {PLANE};
istream& operator>> (istream&, surfacetype&);
//string to_string(surfacetype&);

enum lossmodel {INDEP_loss};
istream& operator>> (istream&, lossmodel&);
//string to_string(lossmodel&);

typedef function<pair<double, double>(double)> lossmodelfunc;

map<lossmodel, lossmodelfunc> LOSSMODELMAP {
	{INDEP_loss, [](const double& p)->pair<double,double>{return {p,0};}}
};

enum noisemodel {INDEP, EM2, GATE};
istream& operator>> (istream&, noisemodel&);
//string to_string(noisemodel&);
typedef function<pair<double, double>(double)> noisemodelfunc;

map<noisemodel, noisemodelfunc> NOISEMODELMAP {
	{EM2, [](const double& p)->pair<double,double>{return {p*(32/15+2),p*4/15};}},
	{INDEP, [](const double& p)->pair<double,double>{return {p,0};}}
};

//################################## vector printing ##############################
template <class T>
ostream& operator<<(ostream&, const vector<T>);

template <class T, class S>
ostream& operator<<(ostream&, const pair<T,S>);

template <class T>
void printMatrix(const vector<T>);

struct coord{
	int x,y,z,l;
	coord(){}
	coord(const int&, const int&, const int&, const int&);
	coord(const int&, const coord&);
	int hash(const coord&);
	void operator=(const coord&);
	
	//dual lattice (cubes)
	int getFaceQubits(const coord&, const int&, const int&);
};

//print coord
ostream& operator<<(ostream&, const coord&);

//coord comparison
bool operator<(const coord&, const coord&);

bool operator==(const coord&, const coord&);

struct vertex {
	vertex(){}
	vertex(const int&, const coord&);
	int c;
	vector<int> partial;
};

int getTaxicabDistance(const coord&, const int&, const int&);

coord getTaxicabDisplacement(const coord&, const int&, const int&);

int getTaxicabDistance(const coord&, const vector<int>&, const vector<int>&);
