//########################## Faster Main function ##########################
//########################## 3D Fast ##########################
//# labeling with single integer
//Surf: PLANE; noisemodel: EM2, GATE (two seperate probabilities)

# include "functions.cpp"
# include "libs/cxxopts/cxxopts.hpp"
# include <assert.h>
# include "libs/blossom5-v2.05.src/PerfectMatching.h"
# include "libs/blossom5-v2.05.src/GEOM/GeomPerfectMatching.h"
# include "cppflow/ops.h"
# include "cppflow/model.h"

using namespace std;

struct coord{
	int x,y,z,l;
	coord(){}
	coord(const int& x, const int& y, const int& z, const int& l){
		this->x = x;
		this->y = y;
		this->z = z;
		this->l = l;
	}
	coord(const int&c, const coord& S){
		int L = S.x;
		int M = S.y;
		int N = S.z;
		*this = {c % (M*L*N) % (M*L) % L, c % (M*L*N) % (M*L) / L, c % (M*L*N) / (M*L), c / (M*L*N)};
	}
	int hash(const coord& S){
		int L = S.x;
		int M = S.y;
		int N = S.z;
		return l*L*M*N + z*L*M + y*L + x;
	}
	void operator=(const coord& c) {
		x = c.x;
		y = c.y;
		z = c.z;
		l = c.l;
	}
	//dual lattice (cubes)
	int getFaceQubits(const coord& S, const int& k){ //new version
		switch(k) {
			case 10:
				return coord(divmod(x + 1, S.x), y, z, 1).hash(S); //10
			case 8:
				return coord(x, y, z, 0).hash(S);//8
			case 9:
				return coord(x, divmod(y + 1, S.y), z, 0).hash(S); //9
			case 0:
				return coord(x, y, divmod(z+1, S.z), 0).hash(S);//0
			case 6:
				return coord(x, y, z, 2).hash(S);//6		
			case 4:
				return coord(divmod(x + 1, S.x), y, z, 2).hash(S);//4
			case 7:
				return coord(x, divmod(y + 1, S.x), z, 2).hash(S);//７
			case 3:
				return coord(x, y, divmod(z+1, S.z), 1).hash(S);//３
			case 11:
				return coord(x, y, z, 1).hash(S);//11
		}
		return 0;
	}
};

class cluster {		
		coord S;
		surfacetype this_surf;
		vector<int> surf;///debug


	public:
		vector<int> stabs;
		vector<int> c_error_pos;
		vector<int> d_error_pos;
		vector<int> c_loss_pos;

		vector<vector<int>> super_chunks;
		vector<int> left_boundary_chunks;
		vector<int> right_boundary_chunks;
		vector<int> boundary_info;

		cluster(const coord&, const surfacetype&);
		void addNoise(const double&, const double&, const noisemodel, const int&);
		void addGateNoise(const double&, const noisemodel, const int&);
		void addLoss(const double&, const lossmodel, const int&);
		void addError(); // manually add error

		void printQubit(surfacetype);

		void getZMeasurements();
		void getXMeasurements();

		int decodeWithMWPM(int, bool, bool);

		// loss decoder helper functions
		int findParent(int);
		void getSuperChunks();
		int getLeftDistance(const coord&, const vector<int>&);
		int getRightDistance(const coord&, const vector<int>&);
		int decodeWithMWPMloss(int, bool, bool);

		// NN decoder helper functions
		vector<float> getWindow(coord&, const bool&);
		void decodeWithNN(cppflow::model, int);
		
		
		void getSurf();//debugging
		void printSurf();//debugging
		void printSuperChunks();

		void surfaceCorrect(PerfectMatching*, const int&, const vector<int>&, const vector<int>&, const int&);///debug
		
		int checkMeasurementOutcomeX1();
		int checkMeasurementOutcomeX2();
};

//===================================================================//

coord getTaxicabDisplacement(const int&, const int&, const surfacetype&);

int getTaxicabDistance(const coord& S, const coord&, const coord&, const surfacetype&);

int getTaxicabDistance(const coord& S, const vector<int>& chunk1, const vector<int>& chunk2, const surfacetype&);

void testDecoding(const int L, const int M, const double p, surfacetype);

void loopDecoding(const int L_range, const int trials, const string, surfacetype, bool);
