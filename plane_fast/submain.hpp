# include <assert.h>
# include "libs/cxxopts/cxxopts.hpp"
# include "subfunctions.hpp"
# include "libs/blossom5-v2.05.src/PerfectMatching.h"
# include "libs/blossom5-v2.05.src/GEOM/GeomPerfectMatching.h"

using namespace std;

struct subcluster {
	static int cc;
	struct subvertex;
	struct subcoord{
		int x,y,l;
		subcoord(){}
		subcoord(const int&, const int&, const int&);
//		subcoord(const int&);
		
		subcoord(const int& c);
		int hash();
		void operator=(const subcluster::subcoord&);
	};
	subcluster(const subcoord&, const subsurfacetype&);

	void addNoise(const double&);

	void printQubit();

	void getx_vec();
	void getXVec();
	int decodeWithMWPM(int);

	//	vector<int> surf;///decode
	void getQubit();///decode
	//	void printSurf();///decode
	//	void surfaceCorrect(PerfectMatching*, const int&, const vector<int>&, const vector<int>&, const int&);///decode
	static subcoord S;
	vector<int> z_error_pos;
private:
	//high level structures
	
	vector<int> x_vec;
	subsurfacetype this_surf = subPLANE;
	
};

int subcluster::cc = 1; // definition

//===================================================================//

subcluster::subcoord getTaxicabDisplacement(const int&, const int&, const subsurfacetype&);

int getTaxicabDistance(const subcluster::subcoord&, const int&, const int&, const subsurfacetype&);

void testDecoding(const int L, const int M, const double p, subsurfacetype);

void loopDecoding(const int L_range, const int trials, const string, subsurfacetype, bool);
