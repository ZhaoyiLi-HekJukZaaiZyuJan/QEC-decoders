//########################## [[3,1,1]] Toric Code ##########################
//########################## main.cpp ##########################

//system libraries
#include <fstream>
#include <functional>
#include <algorithm>
#include <utility>
#include <iomanip>
#include <assert.h>
#include <time.h>

//personal libraries
#include <thread_for.hpp>
#include <cxxopts/cxxopts.hpp>
#include <PerfectMatching.h>
#include <GEOM/GeomPerfectMatching.h>

//project libraries
#include <vertex.hpp> //directory added through -I in compiler
#include <functions.hpp>
#include <cluster.hpp>


using namespace std;

clock_t start_t = clock();
auto start = chrono::steady_clock::now();


class Cluster : public cluster {
	vector<vector<int>> c_error_pos_311;
	public:
	void addPauli311(const double&, const double&, const int&);
	void addNoise(const double&, const double&, const noisemodel, const lossmodel, const int& seed = 0);
	void addGateNoise(const double&, const int& seed=0);
	void addBiasedGateNoise1(const double &, const double &, const int& seed=0);
	void addBiasedGateNoise2(const double &, const double &, const int& seed=0);
	Cluster(const coord&, const surfacetype&);
	void decode311(int);
	void printQubit(surfacetype);
};

Cluster::Cluster(const coord& S, const surfacetype& this_surf) : cluster(S, this_surf){
	vector<vector<int>> c_error_pos_311(3*S.x*S.y*S.z, vector<int>(3, 1));
	this->c_error_pos_311 = c_error_pos_311;
}

void Cluster::addPauli311(const double & p, const double & q, const int& seed){
	//Total heuristic Probabilities
	random_device rd;
	mt19937 engine{rd()};
	if (seed != 0) {
		engine.seed(seed);
	}
	uniform_real_distribution<> dist(0.0, 1.0);
	//Initialization of error operator
	//Error Model 1: uncorrelated error distribution for all physical c_error_poss
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		for (int k = 0; k < 3; k++) {
			if(dist(engine) < p) {
				c_error_pos_311[c][k] = -1;
			} else {
				c_error_pos_311[c][k] = 1;
			}
		}
	}
	//	//Error Model 2: Correlation Error on the Same Edge (check correctness)
	//	for (int c = 0; c < S.x*S.y*S.z; c++) {
	//		vertex aVertex(c, S); //get vertices as when getting x_vec's
	//		for (int pos = 0; pos < 6; pos = pos + 2) {
	//			if(dist(engine) < NOISEMODELMAP[N](p).second) {
	//				c_error_pos[aVertex.partial[pos]] *= -1;
	//				c_error_pos[aVertex.partial[pos + 1]] *= -1;
	//			}
	//		}
	//	}
}

void Cluster::addGateNoise(const double & p, const int& seed){
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		c_error_pos_311[c][0] = 1;
		c_error_pos_311[c][1] = 1;
	}
	//Total heuristic Probabilities
	random_device rd;
	mt19937 engine{rd()};
	if (seed != 0) {
		engine.seed(seed);
	}
	uniform_real_distribution<> dist(0.0, 1.0);

//	//Initialization of error operator
//	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
//		double p1 = dist(engine); //+1 error
//		double p2 = dist(engine); //+2 error
//		double pCX = dist(engine); //CZ error
//
//		if(p1 < p*2/3) {
//			c_error_pos_311[c][1] *= -1;
//		}
//		if(p2 < p*2/3) {
//			c_error_pos_311[c][0] *= -1;
//			c_error_pos_311[c][1] *= -1;
//		}
//		if(pCX < p*4/15) {
//			c_error_pos_311[c][1] *= -1;
//		} else if(p*4/15 < pCX && pCX < p*8/15) {
//			c_error_pos_311[c][0] *= -1;
//		} else if(p*8/15 < pCX && pCX < p*12/15) {
//			c_error_pos_311[c][0] *= -1;
//			c_error_pos_311[c][1] *= -1;
//		}
//	}
	//initialization, measurement and storage errors
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		for (int i = 0; i < 2; i++) {// 311
//			if(dist(engine) < p*2-pow(p,2)*4/3) {
			if(dist(engine) < p*4/3) {
				c_error_pos_311[c][i] = -1;
			} else {
				c_error_pos_311[c][i] = 1;
			}
		}
	}
	
	for (int c = 0; c < S.x*S.y*S.z; c++) {
		coord C(c,S);
		if (C.z == S.z - 1 || C.x == 0) {//remove boundary cubes
			continue;
		}
		for (int i = 0; i < 2; i++) {//311
			for (int face = 0; face < 3; face ++) {
				for (int direction = 0; direction < 2; direction++) {//horizontal/vertical
				
					double p1 = dist(engine); //process 1 (black)
					double p2 = dist(engine); //process 2 (pink)
					double p3 = dist(engine); //process 3 (rose)
					double p4 = dist(engine); //process 4 (purple)
					

					if (p1 < p*8/15) { //black
						c_error_pos_311[C.getFaceQubits(S, face, direction, 3)][i] *= -1;
					}
					if (p4 < p*4/15) { //purple
						c_error_pos_311[C.getFaceQubits(S, face, direction, 0)][i] *= -1;
					} else if (p*4/15 < p4 && p4 < p*8/15) { //purple
						c_error_pos_311[C.getFaceQubits(S, face, direction, 1)][i] *= -1;
						c_error_pos_311[C.getFaceQubits(S, face, direction, 2)][i] *= -1;
						c_error_pos_311[C.getFaceQubits(S, face, direction, 3)][i] *= -1;
					} else if (p*8/15 < p4 && p4 < p*12/15) { //purple
						c_error_pos_311[C.getFaceQubits(S, face, direction, 0)][i] *= -1;
						c_error_pos_311[C.getFaceQubits(S, face, direction, 1)][i] *= -1;
						c_error_pos_311[C.getFaceQubits(S, face, direction, 2)][i] *= -1;
						c_error_pos_311[C.getFaceQubits(S, face, direction, 3)][i] *= -1;
					}
					
					
					if (p2 < p*4/15) { //pink
						c_error_pos_311[C.getFaceQubits(S, face, direction, 2)][i] *= -1;
					} else if (p*4/15 < p2 && p2 < p*8/15) {
						c_error_pos_311[C.getFaceQubits(S, face, direction, 3)][i] *= -1;
					} else if (p*8/15 < p2 && p2 < p*12/15) {
						c_error_pos_311[C.getFaceQubits(S, face, direction, 2)][i] *= -1;
						c_error_pos_311[C.getFaceQubits(S, face, direction, 3)][i] *= -1;
					}

					if (p3 < p*4/15) { //rose
						c_error_pos_311[C.getFaceQubits(S, face, direction, 1)][i] *= -1;
					} else if (p*4/15 < p3 && p3 < p*8/15) {
						c_error_pos_311[C.getFaceQubits(S, face, direction, 1)][i] *= -1;
						c_error_pos_311[C.getFaceQubits(S, face, direction, 2)][i] *= -1;
						c_error_pos_311[C.getFaceQubits(S, face, direction, 3)][i] *= -1;
					} else if (p*8/15 < p3 && p3 < p*12/15) {
						c_error_pos_311[C.getFaceQubits(S, face, direction, 2)][i] *= -1;
						c_error_pos_311[C.getFaceQubits(S, face, direction, 3)][i] *= -1;
					}
					
				}
			}
		}
	}
}


void Cluster::addBiasedGateNoise1(const double & p, const double & B, const int& seed){
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		c_error_pos_311[c][0] = 1;
		c_error_pos_311[c][1] = 1;
		c_error_pos_311[c][2] = 1;
	}
	//Total heuristic Probabilities
	random_device rd;
	mt19937 engine{rd()};
	if (seed != 0) {
		engine.seed(seed);
	}
	uniform_real_distribution<> dist(0.0, 1.0);

//	//Initialization of error operator
//	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
//		double p1 = dist(engine); //+1 error
//		double p2 = dist(engine); //+2 error
//		double pCX = dist(engine); //CZ error
//
//		if(p1 < p*2/3) {
//			c_error_pos_311[c][1] *= -1;
//		}
//		if(p2 < p*2/3) {
//			c_error_pos_311[c][0] *= -1;
//			c_error_pos_311[c][1] *= -1;
//		}
//		if(pCX < p*4/15) {
//			c_error_pos_311[c][1] *= -1;
//		} else if(p*4/15 < pCX && pCX < p*8/15) {
//			c_error_pos_311[c][0] *= -1;
//		} else if(p*8/15 < pCX && pCX < p*12/15) {
//			c_error_pos_311[c][0] *= -1;
//			c_error_pos_311[c][1] *= -1;
//		}
//	}
	//initialization, measurement and storage errors
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		for (int i = 0; i < 2; i++) {// 311
//			if(dist(engine) < p*2-pow(p,2)*4/3) {
			if(dist(engine) < p*4/3) {
				c_error_pos_311[c][i] = -1;
			} else {
				c_error_pos_311[c][i] = 1;
			}
		}
	}
	
	for (int c = 0; c < S.x*S.y*S.z; c++) {
		coord C(c,S);
		if (C.z == S.z - 1 || C.x == 0) {//remove boundary cubes
			continue;
		}
		for (int i = 0; i < 2; i++) {//311
			for (int face = 0; face < 3; face ++) {
				
				double p1x = dist(engine); //process 1 (black)
				double p1z = dist(engine); //process 1 (black)
				double p2x = dist(engine); //process 2 (pink)
				double p2z = dist(engine); //process 2 (pink)
				double p3x = dist(engine); //process 3 (rose)
				double p3z = dist(engine); //process 3 (rose)
				double p4x = dist(engine); //process 4 (purple)
				double p4z = dist(engine); //process 4 (purple)
				

				if (p1z < p*2/3) {
					c_error_pos_311[C.getFaceQubits(S, face, 0, 3)][i] *= -1;					
				}
				if (p1x < p*2/3/B) { //black
					c_error_pos_311[C.getFaceQubits(S, face, 0, 3)][i] *= -1;					
				}	

				if (p2z < p*2/3) { //pink
					c_error_pos_311[C.getFaceQubits(S, face, 0, 2)][i] *= -1;
				}
				if (p2x < p*2/3/B) {
					c_error_pos_311[C.getFaceQubits(S, face, 0, 2)][i] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 3)][i] *= -1;
				}

				if (p3z < p*2/3) { //rose
					c_error_pos_311[C.getFaceQubits(S, face, 0, 1)][i] *= -1;
				}
				if (p3x < p*2/3/B) {
					c_error_pos_311[C.getFaceQubits(S, face, 0, 1)][i] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 2)][i] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 3)][i] *= -1;
				} 

				if (p4z < p*2/3) { //purple
					c_error_pos_311[C.getFaceQubits(S, face, 0, 0)][i] *= -1;
				}
				if (p4x < p*2/3/B) { //purple
					c_error_pos_311[C.getFaceQubits(S, face, 0, 0)][i] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 1)][i] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 2)][i] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 3)][i] *= -1;
				}			
					


				p1x = dist(engine); //process 1 (black)
				p1z = dist(engine); //process 1 (black)
				p2x = dist(engine); //process 2 (pink)
				p2z = dist(engine); //process 2 (pink)
				p3x = dist(engine); //process 3 (rose)
				p3z = dist(engine); //process 3 (rose)
				p4z = dist(engine); //process 4 (purple)
				

				if (p1z < p*2/3) {//black
					c_error_pos_311[C.getFaceQubits(S, face, 0, 3)][i] *= -1;					
				}
				if (p1x < p*2/3/B) { 
					c_error_pos_311[C.getFaceQubits(S, face, 0, 0)][i] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 1)][i] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 2)][i] *= -1;			
				}	

				if (p2z < p*2/3) { //pink
					c_error_pos_311[C.getFaceQubits(S, face, 0, 2)][i] *= -1;
				}
				if (p2x < p*2/3/B) {
					c_error_pos_311[C.getFaceQubits(S, face, 0, 0)][i] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 1)][i] *= -1;
				}

				if (p3z < p*2/3) { //rose
					c_error_pos_311[C.getFaceQubits(S, face, 0, 1)][i] *= -1;
				}
				if (p3x < p*2/3/B) {
					c_error_pos_311[C.getFaceQubits(S, face, 0, 0)][i] *= -1;
				} 

				if (p4z < p*2/3) { //purple
					c_error_pos_311[C.getFaceQubits(S, face, 0, 0)][i] *= -1;
				}				
			}
		}
	}
}


//gate sequence 2 (top layer (logical gate) order)
void Cluster::addBiasedGateNoise2(const double & p, const double & B, const int& seed){
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		c_error_pos_311[c][0] = 1;
		c_error_pos_311[c][1] = 1;
		c_error_pos_311[c][2] = 1;
	}
	//Total heuristic Probabilities
	random_device rd;
	mt19937 engine{rd()};
	if (seed != 0) {
		engine.seed(seed);
	}
	uniform_real_distribution<> dist(0.0, 1.0);

//	//Initialization of error operator
//	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
//		double p1 = dist(engine); //+1 error
//		double p2 = dist(engine); //+2 error
//		double pCX = dist(engine); //CZ error
//
//		if(p1 < p*2/3) {
//			c_error_pos_311[c][1] *= -1;
//		}
//		if(p2 < p*2/3) {
//			c_error_pos_311[c][0] *= -1;
//			c_error_pos_311[c][1] *= -1;
//		}
//		if(pCX < p*4/15) {
//			c_error_pos_311[c][1] *= -1;
//		} else if(p*4/15 < pCX && pCX < p*8/15) {
//			c_error_pos_311[c][0] *= -1;
//		} else if(p*8/15 < pCX && pCX < p*12/15) {
//			c_error_pos_311[c][0] *= -1;
//			c_error_pos_311[c][1] *= -1;
//		}
//	}
	//initialization, measurement and storage errors
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		for (int i = 0; i < 2; i++) {// 311
//			if(dist(engine) < p*2-pow(p,2)*4/3) {
			if(dist(engine) < p*4/3) {
				c_error_pos_311[c][i] = -1;
			} else {
				c_error_pos_311[c][i] = 1;
			}
		}
	}
	
	for (int c = 0; c < S.x*S.y*S.z; c++) {
		coord C(c,S);
		if (C.z == S.z - 1 || C.x == 0) {//remove boundary cubes
			continue;
		}
		for (int i = 0; i < 3; i++) {//311
			for (int face = 0; face < 3; face ++) {
				
				double p1x = dist(engine); //process 1 (black)
				double p1z = dist(engine); //process 1 (black)
				double p2x = dist(engine); //process 2 (pink)
				double p2z = dist(engine); //process 2 (pink)
				double p3x = dist(engine); //process 3 (rose)
				double p3z = dist(engine); //process 3 (rose)
				double p4x = dist(engine); //process 4 (purple)
				double p4z = dist(engine); //process 4 (purple)
				

				//1,3,5,7 1<=>2 swapped as in ref

				if (p1z < p*2/3) {
					c_error_pos_311[C.getFaceQubits(S, face, 0, 3)][divmod(i,3)] *= -1;					
				}
				if (p1x < p*2/3/B) { //black
					c_error_pos_311[C.getFaceQubits(S, face, 0, 3)][divmod(i,3)] *= -1;					
				}	

				if (p2z < p*2/3) { //pink
					c_error_pos_311[C.getFaceQubits(S, face, 0, 2)][divmod(i,3)] *= -1;
				}
				if (p2x < p*2/3/B) {
					cout << "here" << endl;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 2)][divmod(i,3)] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 3)][divmod(i,3)] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 3)][divmod(i+1,3)] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 3)][divmod(i+2,3)] *= -1;
				}

				if (p3z < p*2/3) { //rose
					c_error_pos_311[C.getFaceQubits(S, face, 0, 1)][divmod(i,3)] *= -1;
				}
				if (p3x < p*2/3/B) {
					cout << "here" << endl;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 1)][divmod(i+1,3)] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 1)][divmod(i+2,3)] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 0)][divmod(i,3)] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 0)][divmod(i+1,3)] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 0)][divmod(i+2,3)] *= -1;
				} 

				if (p4z < p*2/3) { //purple
					c_error_pos_311[C.getFaceQubits(S, face, 0, 0)][divmod(i,3)] *= -1;
				}
				if (p4x < p*2/3/B) { //purple
					cout << "here" << endl;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 0)][divmod(i+1,3)] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 0)][divmod(i+2,3)] *= -1;
				}			
					
				//2,4,6,8

				p1x = dist(engine); //process 1 (black)
				p1z = dist(engine); //process 1 (black)
				p2x = dist(engine); //process 2 (pink)
				p2z = dist(engine); //process 2 (pink)
				p3x = dist(engine); //process 3 (rose)
				p3z = dist(engine); //process 3 (rose)
				p4x = dist(engine); //process 4 (purple)
				p4z = dist(engine); //process 4 (purple)
				

				if (p1z < p*2/3) {//black
					c_error_pos_311[C.getFaceQubits(S, face, 0, 3)][divmod(i+1,3)] *= -1;					
				}
				if (p1x < p*2/3/B) { 
					c_error_pos_311[C.getFaceQubits(S, face, 0, 3)][divmod(i,3)] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 3)][divmod(i+1,3)] *= -1;		
				}	

				if (p2z < p*2/3) { //pink
					c_error_pos_311[C.getFaceQubits(S, face, 0, 2)][divmod(i+1,3)] *= -1;
				}
				if (p2x < p*2/3/B) {
					c_error_pos_311[C.getFaceQubits(S, face, 0, 2)][divmod(i,3)] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 2)][divmod(i+1,3)] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 3)][divmod(i,3)] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 3)][divmod(i+1,3)] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 3)][divmod(i+2,3)] *= -1;
				}

				if (p3z < p*2/3) { //rose
					c_error_pos_311[C.getFaceQubits(S, face, 0, 1)][divmod(i+1,3)] *= -1;
				}
				if (p3x < p*2/3/B) {
					c_error_pos_311[C.getFaceQubits(S, face, 0, 1)][divmod(i+2,3)] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 0)][divmod(i,3)] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 0)][divmod(i+1,3)] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 0)][divmod(i+2,3)] *= -1;
				} 

				if (p4z < p*2/3) { //purple
					c_error_pos_311[C.getFaceQubits(S, face, 0, 0)][divmod(i+1,3)] *= -1;
				}	
				if (p4x < p*2/3) { //purple
					c_error_pos_311[C.getFaceQubits(S, face, 0, 0)][divmod(i+2,3)] *= -1;
				}

				p1x = dist(engine); //process 1 (black)
				p1z = dist(engine); //process 1 (black)
				p2x = dist(engine); //process 2 (pink)
				p2z = dist(engine); //process 2 (pink)
				p3x = dist(engine); //process 3 (rose)
				p3z = dist(engine); //process 3 (rose)
				p4z = dist(engine); //process 4 (purple)

				if (p1z < p*2/3) {//black
					c_error_pos_311[C.getFaceQubits(S, face, 0, 3)][divmod(i+2,3)] *= -1;					
				}
				if (p1x < p*2/3/B) { 
					c_error_pos_311[C.getFaceQubits(S, face, 0, 3)][divmod(i,3)] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 3)][divmod(i+1,3)] *= -1;	
					c_error_pos_311[C.getFaceQubits(S, face, 0, 3)][divmod(i+2,3)] *= -1;	
				}	

				if (p2z < p*2/3) { //pink
					c_error_pos_311[C.getFaceQubits(S, face, 0, 2)][divmod(i+2,3)] *= -1;
				}
				if (p2x < p*2/3/B) {
					c_error_pos_311[C.getFaceQubits(S, face, 0, 2)][divmod(i,3)] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 2)][divmod(i+1,3)] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 2)][divmod(i+2,3)] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 3)][divmod(i,3)] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 3)][divmod(i+1,3)] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 3)][divmod(i+2,3)] *= -1;
				}

				if (p3z < p*2/3) { //rose
					c_error_pos_311[C.getFaceQubits(S, face, 0, 1)][divmod(i+2,3)] *= -1;
				}
				if (p3x < p*2/3/B) {
					c_error_pos_311[C.getFaceQubits(S, face, 0, 0)][divmod(i,3)] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 0)][divmod(i+1,3)] *= -1;
					c_error_pos_311[C.getFaceQubits(S, face, 0, 0)][divmod(i+2,3)] *= -1;
				} 

				if (p4z < p*2/3) { //purple
					c_error_pos_311[C.getFaceQubits(S, face, 0, 0)][divmod(i+2,3)] *= -1;
				}			
			}
		}
	}
}



//No seed used for seed = 0
void Cluster::addNoise(const double & p, const double & q, const noisemodel N, const lossmodel L, const int& seed){
	if (N == GATE) { //p:error probability, //q: bias
		addGateNoise(p, seed);
	}
	if (N == GATE_biased) { //p:error probability, //q: bias
		addBiasedGateNoise2(p, q, seed);
	}
	if (N != GATE && N != GATE_biased) {
		addPauli311(NOISEMODELMAP[N](p,0).first, NOISEMODELMAP[N](p,0).second,  seed);
	} 

	// Loss part
	if (L == OFF_loss){}
	if (L == toLoss) { //convert measured syndrom to loss errors
		getStabs();
		for (int c = 0; c < S.x*S.y*S.z; c++){
			if (stabs[c] == -1) {
				vertex Vertex(c, S);
				for (int pos = 0; pos < 6; pos ++) {
					c_loss_pos[Vertex.partial[pos]] = -1 ? 1 : -1;
				}
			}
		}
	}
	if (L != OFF_loss && L != toLoss){
		addLoss(LOSSMODELMAP[L](q,0).first, LOSSMODELMAP[L](q,0).second, seed);
	}
}

void Cluster::decode311(int verbosity = 0){
	//initialization loss and z_error
	
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		//stabalizer measurement (xx)
		int xxi = c_error_pos_311[c][0]*c_error_pos_311[c][1]; //stabalizer operator
		int ixx = c_error_pos_311[c][1]*c_error_pos_311[c][2]; //stabalizer operator
		
		if (xxi == -1 && ixx == -1) {//correction scheme
			// c_error_pos_311[c][1] *= -1;
		} else if (xxi == -1) {//correction scheme
			c_error_pos_311[c][0] *= -1;
		} else if (ixx == -1){//correction scheme
			// c_error_pos_311[c][2] *= -1;
		}
		
		//pauli measurement (xii)
		c_error_pos[c] = c_error_pos_311[c][0];//erronous logical qubit
	}
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		coord C(c, S);
		if(this_surf == PLANE && ((C.l == 1 && (C.y == S.y - 1 || C.x == 0)) ||(C.l == 2 && (C.z == S.z - 1 || C.x == 0)))){
			c_error_pos[c] = 1;
			c_loss_pos[c] = 1;
		}//	for planar code, remove errors on left, more, and lower boundaries to create edges.
	}
}

void Cluster::printQubit(surfacetype s = PLANE){
	getStabs();
	for(int k = 0; k < 2*S.z; k++){
		vector<vector<string>> print_out;
		if(s == TORUS){
			for(int i = 0; i < 2*S.x; i++){
				vector<string> a_line;
				for(int j = 0; j < 2*S.y; j++){
					a_line.push_back(" ");
				}
				print_out.push_back(a_line);
			}
		} else if(s == PLANE){
			for(int i = 0; i < 2*S.y - 1; i++){
				vector<string> a_line;
				for(int j = 0; j < 2*S.x; j++){
					a_line.push_back(" ");
				}
				print_out.push_back(a_line);
			}
		}	
		if (k % 2 == 0) {//even layer
			for (int x = 0; x < S.x; x++) {
				for (int y = 0; y < S.y; y++) {
					if (c_error_pos[coord(x, y, k/2, 0).hash(S)] < 0){
						print_out[2 * y][2 * x + 1] = "\033[1;96mZ\033[0m";
					} else if (c_error_pos[coord(x, y, k/2, 0).hash(S)] == 2){
						print_out[2 * y][2 * x + 1]  = "\033[1;33mC\033[0m";
					}
					if (c_error_pos[coord(x, y, k/2, 1).hash(S)] < 0){
						print_out[2 * y + 1][2 * x] = "\033[1;96mZ\033[0m";
					} else if (c_error_pos[coord(x, y, k/2, 1).hash(S)] == 2){
						print_out[2 * y + 1][2 * x]  = "\033[1;33mC\033[0m";
					}
					if (s == PLANE && x == 0){
						print_out[2 * y][2 * x] = " ";
					} else if (stabs[coord(x, y, k/2, 0).hash(S)] > 2){
						print_out[2 * y][2 * x] = to_string(stabs[coord(x, y, k/2, 0).hash(S)]);
					} else if (stabs[coord(x, y, k/2, 0).hash(S)] > 0){
						print_out[2 * y][2 * x] = "◯";
					} else {
						print_out[2 * y][2 * x] = "\033[1;92m⊕\033[0m";
					}		
					
				}
			}
		} else {//odd layer
			if (s == PLANE && k == 2*S.z-1){
				continue;
			}
			for (int x = 0; x < S.x; x++) {
				for (int y = 0; y < S.y; y++) {
					if (c_error_pos[coord(x, y, k/2, 2).hash(S)] < 0){
						print_out[2 * y][2 * x] = "\033[1;96mZ\033[0m";
					} else if (c_error_pos[coord(x, y, k/2, 2).hash(S)] == 2){
						print_out[2 * y][2 * x]   = "\033[1;33mC\033[0m";
					}
					// if (s == PLANE && y == S.y-1){
					// } else if (stabs[coord(x, y, k/2, 0).hash(S) + S.x*S.y*S.z] > 2){
					// 	print_out[2*y+1][2*x+1] = to_string(stabs[coord(x, y, k/2, 0).hash(S) + S.x*S.y*S.z]);
					// } else if (stabs[coord(x, y, k/2, 0).hash(S) + S.x*S.y*S.z] > 0){
					// 	print_out[2*y+1][2*x+1] = "◯";
					// } else {
					// 	print_out[2*y+1][2*x+1] = "\033[1;92m⊕\033[0m";
					// }		
				}
			}
		}
		cout << "layer:" << k << endl;
		printMatrix(print_out);
	}
	
}

void testDecoding(Cluster& test_cluster, const double& p, const double& q, const int& seed, surfacetype surf, noisemodel N, lossmodel L, int verbosity = 0){
	start_t = clock();
	//add noise
	test_cluster.addNoise(p, q, N, L, seed);
	test_cluster.decode311(0);
	// test_cluster.addError();
	
	cout << "t(Generation):"<< double(clock()-start_t)/CLOCKS_PER_SEC << endl;//timing
	start_t = clock();
	try {
		test_cluster.getSuperChunks();
	} catch (...) {
		cout << "failure" << endl;
		return;
	}
	
	if (verbosity >= 1){
		if (verbosity == 2) {
			test_cluster.printSuperChunks();
			test_cluster.printQubit(surf);
			cout << "here" << endl;
			test_cluster.getSurf();
			cout << "here2" << endl;
			test_cluster.printSurf();
		}
		cout << "t(getSuperChunks):" << double(clock()-start_t)/CLOCKS_PER_SEC << endl;//timing
		start_t = clock();
	}
	cout << "there" <<endl;
	int parity = test_cluster.decodeWithMWPMLoss(verbosity, 1, 1, surf);
	
	if (verbosity >= 1){
		if (verbosity == 2) {
			cout << "here2" <<endl;
			test_cluster.printQubit();
			test_cluster.printSurf();
		}
		cout << "t(decodeWithMWPMLoss):"<< double(clock()-start_t)/CLOCKS_PER_SEC << endl;//timing
		start_t = clock();
	}
	
	cout << "X1 check:" << parity <<endl;
}

int loopDecoding(const int lmin, const int lmax, const int trials, const double pmin, const double pmax, const int Np, const double qmin, const double qmax, const int Nq, const string fname, surfacetype surf, noisemodel N, lossmodel L, bool out, int verbosity = 0, int thread = 0, bool make_corrections = 0){
	cout << "hardware_concurrency" << thread::hardware_concurrency() <<endl;
	ofstream outfile;
	if(out != 0){
		outfile.open(fname);
		outfile << "L,p_error,num_success\n";
	}
	
	int binsp = 0 == Np ? 1:Np;
	int binsq = 0 == Nq ? 1:Nq;
	int num_correct;
	
	for (int j = 0; j <= Nq; j++) {
		double q = qmin + (qmax-qmin)/binsq * j;
		for (int l = lmin; l <= lmax; l = l+2) {
			cout << endl;
			Cluster test_cluster({l,l,l,0}, surf);
			for (int i = 0; i <= Np; i++) {
				double p = pmin + (pmax-pmin)/binsp * i;
				if (!thread) {
					Cluster test_cluster({l,l,l,0}, surf);
				}
				
				//run simulation
				num_correct = 0;
				if (thread) {
					parallel_for(trials, [&](int start, int end){
						Cluster test_cluster({l,l,l,0}, surf);
						for(int k = start; k < end; ++k){
							test_cluster.addNoise(p,q, N,L, 0);
							test_cluster.decode311(0);
							try {
								test_cluster.getSuperChunks();
							} catch (...) {
								continue;
							}
							if (surf == PLANE && test_cluster.decodeWithMWPMLoss(verbosity,make_corrections,0,surf) == 1) {
								num_correct ++; //correction successful
							}
						}
					});
				} else {
					for(int k = 0; k < trials; ++k){
						test_cluster.addNoise(p,q, N,L, 0);
						test_cluster.decode311(0);
						try {
							test_cluster.getSuperChunks();
						} catch (...) {
							continue;
						}
						if (surf == PLANE && test_cluster.decodeWithMWPMLoss(verbosity,make_corrections,0,surf) == 1) {
							num_correct ++; //correction successful
						}
					}
				}
				
				//printout/outfile
				if (verbosity >= 2) {
					cout << "t=" << double(clock()-start_t)/CLOCKS_PER_SEC;
					start_t = clock();
				} if (verbosity >= 1){
					cout << l << "," << p << "," << q << "," << num_correct << endl;
				} else {
					cout << ".";
					cout.flush();
				}

				if (out != 0){
					outfile << l << "," << p << "," << q << "," << num_correct << "\n";
					outfile.close();
				}			
			}
		}
	}
	auto end = chrono::steady_clock::now();
	auto diff = end - start;
	cout << "total CPU time:" << double(clock()-start_t)/CLOCKS_PER_SEC <<";total time:"<< chrono::duration <double, milli> (diff).count() << "ms" << endl;
	return num_correct;
}

int main(int argc, const char *argv[]) {
	cout << "hardware_concurrency" << thread::hardware_concurrency() <<endl;
	string fname;
	surfacetype s;
	noisemodel N;
	lossmodel L;
	bool test, use_env, thread, make_corrections, out;
	int verbosity, return_value = 0;
	int lmin, lmax, n, Np, Nq, seed, times;
	float pmin, pmax, qmin, qmax;
	
	//getting options
	cxxopts::Options options(*argv,
							 "Simulator for fault-tolerant measurement-based quantum "
							 "computation on foliated surface code Cluster states"
							 );
	options.add_options()
	("fname", "filename", cxxopts::value(fname)->default_value(""))
	("out", "output in this directory as an .out file", cxxopts::value(out)->default_value("0"))
	("s", "surface type", cxxopts::value(s)->default_value("PLANE"))
	("lmin", "Minimal size of mesh", cxxopts::value(lmin)->default_value("3"))
	("lmax", "Maximal size of mesh", cxxopts::value(lmax)->default_value("17"))
	("n", "Number of trials", cxxopts::value(n)->default_value("100"))
	
	("N", "noise model", cxxopts::value(N)->default_value("GATE"))
	("Np", "z error p Points", cxxopts::value(Np)->default_value("10"))
	("pmin", "Minimal z error probability", cxxopts::value(pmin)->default_value("0.001"))
	("pmax", "Maximal z error probability", cxxopts::value(pmax)->default_value("0.008"))
	
	("L", "loss model", cxxopts::value(L)->default_value("OFF_loss"))
	("Nq", "loss p points", cxxopts::value(Nq)->default_value("1"))
	("qmin", "Minimal loss probability", cxxopts::value(qmin)->default_value("0"))
	("qmax", "Maximal loss probability", cxxopts::value(qmax)->default_value("0"))
	
	("v, verbosity", "verbosity switch", cxxopts::value(verbosity)->default_value("0"))
	("seed", "seed switch", cxxopts::value(seed)->default_value("0"))
	("thread", "thread switch", cxxopts::value(thread)->default_value("0"))
	("use_env", "use environment variables", cxxopts::value(use_env)->default_value("0"))
	("times", "test times switch", cxxopts::value(times)->default_value("1"))
	("make_corrections", "correct qubit", cxxopts::value(make_corrections)->default_value("0"))
	("test", "test switch", cxxopts::value(test)->default_value("0"));
	options.parse(argc, argv);
	
	cout << "hardware_concurrency" << thread::hardware_concurrency() << endl;
	if (use_env) {
		cout << "SLURM_CPUS_PER_TASK:" << atoi(getenv("SLURM_CPUS_PER_TASK")) << endl;
	}

	//outputing options
	cout << "lmin:" << lmin << ";lmax:" << lmax << endl;
	cout << "Pmin:" << pmin << ";Pmax:" << pmax << ";nP" << Np << endl;
	cout << "qmin:" << qmin << ";qmax:" << qmax << ";Nq" << Nq << endl;
	cout << "n:" << n << ";seed:" << seed << ";thread" << thread << endl;
	cout << "N:" << N << ";verbosity:" << verbosity << "L;" << L << endl;
	
	
	if (fname == "") {
		fname = "L=" + to_string(lmin) + ",P=(" + to_string(pmin).substr(3,2) + "," + to_string(pmax).substr(3,2) + "),n=" +to_string(n) + to_string(s) + "," + to_string(N) + ".out";
	}
	if (test) {
		Cluster test_cluster({lmin,lmin,lmin,0}, s);
		int i = 0;
		do {
			testDecoding(test_cluster, pmin, qmin, seed, s, N, L, verbosity);
			i++;
		}
		while (i < times);
	} else{
		return_value = loopDecoding(lmin, lmax, n, pmin, pmax, Np-1, qmin, qmax, Nq-1, fname, s, N, L, out, verbosity, thread, make_corrections);
	}
	return return_value;
}
///################ [[2,1,1]] ################

///######## 1D run ########

//### INDEP run:
	//(INDEP)  ###(p_th ~0.103)
		///./simulate -N INDEP -n 10000 --pmin 0.06 --pmax 0.12 --lmin 3 --Np 30 -v 1
	//(GATE) ###(p_th ~0.008?)
		///./simulate -N GATE -n 1000 --pmin 0.09 --pmax 0.012 --lmin 3 -v 1
	//(GATE_biased (1)) run \beta = 1000 (*p_ref = 1.37% ArXiv 1308.4776 *)
		///./simulate --qmin 1000 --pmin 0.01 --pmax 0.016 --Np 30 --Nq 1 -n 10000 --lmin 3 --lmax 21 -v 1 -N GATE_biased

//######## test ########
	///./simulate -N INDEP -n 10000 --pmin 0.06 --pmax 0.09 --lmin 3 -v 1 --test ###(pE ~0.8)
//######## timing test ########
///./simulate -s PLANE --pmin 0.01 --pmax 0.05  --Np 10 --Nq 1 -n 500 --lmin 3 -v 1

//######## test large ########
	//(simple time)
	///./simulate -s PLANE --pmin 0.01 --pmax 0.05 --Np 10 --Nq 1 -n 500 --test 1 --lmin 20 -v 2

	//(long time)
	//./simulate -s PLANE --qmin 0.05 --qmax 0.05 --pmin 0.004 --pmax 0.008  --Np 10 --Nq 1 -n 500 --lmin 30 -v 1 --test 1
