//defines the cluster and subcluster(2D) classes

#pragma once

# include "coord.hpp"
# include "types.hpp"
# include "../src/libs/blossom5-v2.05.src/PerfectMatching.h"
# include "../src/libs/cppflow-master/include/cppflow/ops.h"
# include "../src/libs/cppflow-master/include/cppflow/model.h"

class cluster {
	public:
	vector<int> c_error_pos;
	vector<int> d_error_pos;
	vector<int> stabs;
	//loss, superchecks, and superchunks
	
	vector<int> c_loss_pos;
	vector<int> d_loss_pos;
	
	vector<vector<int>> super_chunks;
	vector<int> left_boundary_chunks;
	vector<int> right_boundary_chunks;
	vector<int> boundary_info;
	
	//concatenation structures
	vector<vector<int>> c_error_pos_211;

	coord S;
	surfacetype this_surf = PLANE;

	cluster(const coord&, const surfacetype&);
	void addNoise(const double&, const double&, const noisemodel, const lossmodel, const int& seed =0);
	void addPauli(const double&, const double&, const int& seed = 0);
	//add partially probability combined noise
	void addGateNoise(const double&, const int&);
	//add noise respectively for each gate applications
	void addFullGateNoise(const double&, const int&);
	
	void addBiasedGateNoise(const double&, const double&, const int&);
	void addError(); // manually add error
	void printPrimal(surfacetype);
	
	void addLoss(const double&, const double&, const int&);
	void getSuperChunks();
	void printSuperChunks();
	int findParent(int);
	int getLeftDistance(const coord&, const vector<int>&);
	int getRightDistance(const coord&, const vector<int>&);
	
	void getStabs();
	int decodeWithMWPM(int, bool, surfacetype);
	int decodeWithMWPMFull(int verbosity = 0, bool dir = 0, bool make_corrections = 0);
	int decodeWithMWPMLoss(int, bool, surfacetype);
	
	vector<int> surf;///decode
	void getSurf();///decode
	void printSurf();///decode
	void surfaceCorrect(PerfectMatching*, const int&, const vector<int>&, const vector<int>&, const int& verbose, const surfacetype&);
	void surfaceCorrectLoss(PerfectMatching*, const int&, const vector<int>&, const vector<int>&, const int&, const surfacetype&);///decode
	
	int checkMeasurementOutcomeX1();
	int checkMeasurementOutcomeX2();
};


class subcluster {
	public:
	subcoord S;
	subcluster(const subcoord&, const subsurfacetype&);

	void addNoise(const double&);
	void addNoiseWithX(const double&, const noisemodel , const int& seed=0);
	void clearNoise();
	void addNoiseManually(const int&);

	void printQubit();
	void printQubitWithWindow(int);
	void addError();

	void getx_measurements();
	void getz_measurements();
	int decodeWithMWPM(int);
	vector<int> decodeWithMWPMFull(int verbose = 0, bool dir = 0, bool make_corrections = 0);
	void decodeWithNN(cppflow::model model, bool binary_output, int verbose = 0, double cutoff = 1);
	void Predecode(bool binary_output, int verbose = 0);

	//NN Decoder functions
	vector<float> getWindow(const subcoord&);

	//	vector<int> surf;///decode
	void getQubit();///decode
	//	void printSurf();///decode
	//	void surfaceCorrect(PerfectMatching*, const int&, const vector<int>&, const vector<int>&, const int&);///decode
	//high level structures
	vector<int> z_error_pos;
	vector<int> x_error_pos;
	vector<int> stabs;
	subsurfacetype this_surf = subPLANE;

	
	
};