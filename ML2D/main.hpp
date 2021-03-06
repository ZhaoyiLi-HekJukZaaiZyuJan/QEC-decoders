//===================================================================//
//===================  ML Assisted decoding     =====================//
//=====================       main.hpp      =========================//
//===================================================================//

# include "functions.cpp"
# include "../src/libs/cxxopts/cxxopts.hpp"
# include <assert.h>
# include "../src/libs/blossom5-v2.05.src/PerfectMatching.h"
# include "../src/libs/blossom5-v2.05.src/GEOM/GeomPerfectMatching.h"
# include "cppflow/ops.h"
# include "cppflow/model.h"

using namespace std;

struct subcluster {
	struct subcoord{
		int x,y,l;
		subcoord(){}
		subcoord(const int&, const int&, const int&);
		//		subcoord(const int&);
		
		subcoord(const int& c);
		subcoord(const int& c, const int&);
		int hash();
		void operator=(const subcluster::subcoord&);
	};
	static subcoord S;
	
	struct subvertex_x;
	struct subvertex_z;
	subcluster(const subcoord&, const subsurfacetype&);

	void addNoise(const double&, const noisemodel, const int &);
	void addError();

	void printQubit(int);
	vector<int> decodeWithMWPM(int, bool, bool);
	void decodeWithNN(cppflow::model, bool, int, double);
	void Predecode(bool, int);

	//	vector<int> surf;///decode
	void getQubit();///decode
	//	void printSurf();///decode
	//	void surfaceCorrect(PerfectMatching*, const int&, const vector<int>&, const vector<int>&, const int&);///decode
	void getx_measurements();
	void getz_measurements();
	vector<float> getWindow(const subcoord&);
	vector<int> z_error_pos;
	vector<int> x_error_pos;
	vector<int> stabs;
private:

	//high level structures
	subsurfacetype this_surf = PLANE;
};

//===================================================================//

subcluster::subcoord getTaxicabDisplacement(const int&, const int&, const subsurfacetype&);

int getTaxicabDistance(const subcluster::subcoord&, const int&, const int&, const subsurfacetype&);

void loopDecoding(const int L_range, const int trials, const string, subsurfacetype, bool);
