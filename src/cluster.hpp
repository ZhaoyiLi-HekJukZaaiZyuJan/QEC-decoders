#pragma once

# include "coord.hpp"
# include "types.hpp"
# include "../src/libs/blossom5-v2.05.src/PerfectMatching.h"

class cluster {
	public:
	vector<int> c_error_pos;
	vector<int> stabs;
	//loss, superchecks, and superchunks
	vector<int> c_loss_pos;
	
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
	void addGateNoise(const double&, const int&);
	
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
	int decodeWithMWPMLoss(int, int, bool, surfacetype);
	
	vector<int> surf;///decode
	void getSurf();///decode
	void printSurf();///decode
	void surfaceCorrect(PerfectMatching*, const int&, const vector<int>&, const vector<int>&, const int&, const surfacetype&);///decode
	
	int checkMeasurementOutcomeX1();
	int checkMeasurementOutcomeX2();

};