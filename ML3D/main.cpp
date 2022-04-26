
//===================================================================//
//===================  ML Assisted decoding (3D)=====================//
//=====================       main.hpp      =========================//
//===================================================================//

//# based on 3D fast
//Surf: PLANE; noisemodel: EM2, GATE (two seperate probabilities)

#include "main.hpp"
#include "thread_for.hpp"
#include <filesystem>
#include <iostream>

using namespace std;
int Sw = 5;
clock_t start_t = clock();
auto start = chrono::steady_clock::now();	
//for chunk graphics
vector<string> colorvector= { "\033[1;31m◯\033[0m","\033[1;32m◯\033[0m","\033[1;33m◯\033[0m","\033[1;34m◯\033[0m","\033[1;35m◯\033[0m","\033[1;36m◯\033[0m","\033[1;37m◯\033[0m"};

//print coord
ostream& operator<<(ostream& os, const coord& c) {
	os << "(" << c.x << "," << c.y << "," << c.z << "," << c.l << ")";
	return os;
}

//coord comparison
bool operator<(const coord& c1, const coord& c2) {
	if (c1.x < c2.x) {
		return true;
	} else if (c1.x > c2.x) {
		return false;
	} else {
		if (c1.y < c2.y){
			return true;
		} else if (c1.y > c2.y){
			return false;
		} else {
			if (c1.z < c2.z) {
				return true;
			} else {
				return false;
			}
		}
	}
}

bool operator==(const coord& c1, const coord& c2) {
	if (c1.x == c2.x && c1.y == c2.y && c1.z == c2.z && c1.l == c2.l) return true;
	else return false;
}

struct vertex {
	vertex(){}
	vertex(const int& h, const coord& S){
		coord c(h,S);
		this->c = c;
		physical_qubits ={
			coord(divmod(c.x-1, S.x), c.y, c.z, 0).hash(S), coord(c.x, c.y, c.z, 0).hash(S),
			coord(c.x, divmod(c.y-1, S.y), c.z, 1).hash(S), coord(c.x, c.y, c.z, 1).hash(S),
			coord(c.x, c.y, divmod(c.z-1, S.z), 2).hash(S), coord(c.x, c.y, c.z, 2).hash(S),
		};
		//each vertex contains 6 physical_qubits, left (Mx-), right(Mx+), up (Ly-), down (Ly+), less (Nz-), more (Nz+)
	};
	coord c;
	vector<int> physical_qubits;
};


cluster::cluster(const coord& S, const surfacetype& this_surf){
	this->S = S;
	this->this_surf = this_surf;
	
	vector<int> surf(2*S.x*S.y*S.z,1);///debug
	this->surf = surf;///debug
	
	vector<int> d_error_pos(3*S.x*S.y*S.z,1);
	this->d_error_pos = d_error_pos;
	vector<int> c_error_pos(3*S.x*S.y*S.z,1);
	this->c_error_pos = c_error_pos;
	vector<int> c_loss_pos(3*S.x*S.y*S.z,1);
	this->c_loss_pos = c_loss_pos;
	vector<int> stabs(2*S.x*S.y*S.z,1);
	this->stabs = stabs;
}

//differentiate correlated and uncorrelated errors
void cluster::addNoise(const double & p, const double & q, const noisemodel N, const int& seed=0){
	//Total heuristic Probabilities
	random_device rd;
	mt19937 engine{rd()};
	if (seed != 0) {
		engine.seed(seed);
	}
	uniform_real_distribution<> dist(0.0, 1.0);

	if (N == DEPOL1){ //2D correspondence
		for (int c = 0; c < 2*S.x*S.y*S.z; c++) {
			double r = dist(engine);
			if(r < p*2/3) {
				c_error_pos[c] = -1;
			} else {
				c_error_pos[c] = 1;
			}
			coord C(c,S);
			if(r > p*1/3 && r < p) {
				if (C.l == 0){
					d_error_pos[coord(C.x,divmod(C.y-1,S.y),C.z,1).hash(S)] = -1;
				} if (C.l == 1){
					d_error_pos[coord(divmod(C.x-1,S.x),C.y,C.z,0).hash(S)] = -1;
				}
			}else {
				if (C.l == 0){
					d_error_pos[coord(C.x,divmod(C.y-1,S.y),C.z,1).hash(S)] = 1;
				} if (C.l == 1){
					d_error_pos[coord(divmod(C.x-1,S.x),C.y,C.z,0).hash(S)] = 1;
				}
			}
		}
		for (int c = 2*S.x*S.y*S.z; c < 3*S.x*S.y*S.z; c++) { //depol noise
			double r = dist(engine);
			if(r < p*2/3) {
				c_error_pos[c] = -1;
			} else {
				c_error_pos[c] = 1;
			}
			r = dist(engine);
			if(r > p*1/3 && r < p) {
				d_error_pos[c] = -1;
			} else {
				d_error_pos[c] = 1;
			}
		}
	} else if (N == DEPOL2){
		for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
			double r = dist(engine);
			if(r < p*2/3) {
				c_error_pos[c] = -1;
			} else {
				c_error_pos[c] = 1;
			}
			if(r > p*1/3 && r < p) {
				d_error_pos[c] = -1;
			}else {
				d_error_pos[c] = 1;
			}
		}
	} else {
		//Initialization of error operator
		//q1: uncorrelated error distribution for all physical c_error_poss
		for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
			//c_error_pos[c] = 1;
			if(dist(engine) < NOISEMODELMAP[N](p, q).first) {
			c_error_pos[c] = -1;
			} else {
			c_error_pos[c] = 1;
			}
		}

		for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
			//d_error_pos[c] = 1;
			if(dist(engine) < NOISEMODELMAP[N](p, q).first) {
			d_error_pos[c] = -1;
			} else {
			d_error_pos[c] = 1;
			}
		}

		//Error Model 2: Correlation Error on the Same face (opposite edges)
		//wrong 不捣乱就可
		for (int c = 0; c < S.x*S.y*S.z; c++) {
			coord C(c,S);
			for (int face = 0; face < 3; face ++) {
				double p1 = dist(engine);
				double p2 = dist(engine);
				if (face == 0){
					if (p1 < NOISEMODELMAP[N](p, q).second) { //vertical
						c_error_pos[C.getFaceQubits(S, 10)] *= -1;//3
						c_error_pos[C.getFaceQubits(S, 11)] *= -1;//2
					}
					// if (p2 < NOISEMODELMAP[N](p, q).second) { //horizontal
					// 	c_error_pos[C.getFaceQubits(S, 8)] *= -1;//1
					// 	c_error_pos[C.getFaceQubits(S, 9)] *= -1;//0
					// }
					if (p1 < NOISEMODELMAP[N](p, q).second) { //vertical
						d_error_pos[C.getFaceQubits(S, 3)] *= -1;//3
						d_error_pos[C.getFaceQubits(S, 2)] *= -1;//2
					}
					// if (p2 < NOISEMODELMAP[N](p, q).second) { //horizontal
					// 	d_error_pos[C.getFaceQubits(S, 1)] *= -1;//1
					// 	d_error_pos[C.getFaceQubits(S, 0)] *= -1;//0
					// }
				} else if (face == 1){
					if (p1 < NOISEMODELMAP[N](p, q).second) { //vertical
						c_error_pos[C.getFaceQubits(S, 10)] *= -1;//3
						c_error_pos[C.getFaceQubits(S, 6)] *= -1;//2
					}
					// if (p2 < NOISEMODELMAP[N](p, q).second) { //horizontal
					// 	c_error_pos[C.getFaceQubits(S, 3)] *= -1;//1
					// 	c_error_pos[C.getFaceQubits(S, 11)] *= -1;//0
					// }
				} else if (face == 2){
					if (p1 < NOISEMODELMAP[N](p, q).second) { //vertical
						c_error_pos[C.getFaceQubits(S, 0)] *= -1;//3
						c_error_pos[C.getFaceQubits(S, 8)] *= -1;//2
					}
					// if (p2 < NOISEMODELMAP[N](p, q).second) { //horizontal
					// 	c_error_pos[C.getFaceQubits(S, 6)] *= -1;//1
					// 	c_error_pos[C.getFaceQubits(S, 4)] *= -1;//0
					// }
				}
				
			}
		}
	}

	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		coord C(c, S);
		if(this_surf == PLANE && ((C.l == 1 && (C.y == S.y - 1 || C.x == 0)) ||(C.l == 2 && (C.z == S.z - 1 || C.x == 0)))){
			c_error_pos[c] = 1;
		}//	for planar code, remove errors on left, more, and lower boundaries to create edges.
		if(this_surf == PLANE && ((C.l == 0 && (C.y == S.y - 1 || C.x == S.x - 1 || C.z == S.z - 1)) || (C.l == 1 && C.z == S.z - 1) || (C.l == 2 && C.y == S.y-1) )) {
			d_error_pos[c] = 1;
		}
	}
}

void cluster::addError(){ //for debugging
	//Total heuristic Probabilities	
	//add errors manually here
	coord C1(3,3,3,0);
	c_error_pos[C1.hash(S)] *= -1;
	coord C2(3,2,3,0);
	c_error_pos[C2.hash(S)] *= -1;


	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		coord C(c, S);
		if(this_surf == PLANE && ((C.l == 1 && (C.y == S.y - 1 || C.x == 0)) ||(C.l == 2 && (C.z == S.z - 1 || C.x == 0)))){
			c_error_pos[c] = 1;
		}//	for planar code, remove errors on left, more, and lower boundaries to create edges.
	}
}

void cluster::addGateNoise(const double & p, const noisemodel N, const int& seed=0){
	//Total heuristic Probabilities
	random_device rd;
	mt19937 engine{rd()};
	if (seed != 0) {
		engine.seed(seed);
	}
	uniform_real_distribution<> dist(0.0, 1.0);

	//Error Model 1: uncorrelated error distribution for all physical c_error_poss
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		c_error_pos[c] = 1;
	}
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		d_error_pos[c] = 1;
	}

	//Initialization of error operator
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		if(dist(engine) < 2*p) {
			c_error_pos[c] = -1;
		} else {
			c_error_pos[c] = 1;
		}
	}
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		if(dist(engine) < 2*p) {
			d_error_pos[c] = -1;
		} else {
			d_error_pos[c] = 1;
		}
	}


	
	for (int c = 0; c < S.x*S.y*S.z; c++) {
		//probability to erro labels
		for (int order = 0; order < 12; order ++){
			if (dist(engine) > p*16/15){
				continue;
			}
			// cout << order << endl;
			double pa = dist(engine);
			double pb = dist(engine);
			bool aX = 0, aZ = 0, bX = 0, bZ = 0;
			coord C(c,S);
			
			if (pa < 0.5) aX = 1;
			if (pb < 0.5) bX = 1;
			if (0.25 < pa && pa < 0.75) aZ = 1;
			if (0.25 < pb && pb < 0.75) bZ = 1;
			if (order == 0){ 
				//XY gate 1
				if (C.z/2){ 
					if (aX) c_error_pos[C.getFaceQubits(S, 10)] *= -1; //x on blue qubit
					if (bZ) c_error_pos[C.getFaceQubits(S, 10)] *= -1; //z on red qubit
					coord D(C.x,C.y,divmod(C.z-1,S.z),3);
					if (aZ) d_error_pos[D.getFaceQubits(S, 6)] *= -1;
					if (bX) d_error_pos[D.getFaceQubits(S, 6)] *= -1;
				} else {
					if (aX) c_error_pos[C.getFaceQubits(S, 11)] *= -1;
					if (bZ) c_error_pos[C.getFaceQubits(S, 11)] *= -1;
					coord D(divmod(C.x-1,S.x),C.y,divmod(C.z-1,S.z),3);
					if (aZ) d_error_pos[D.getFaceQubits(S, 6)] *= -1;
					if (bX) d_error_pos[D.getFaceQubits(S, 6)] *= -1;
				}
			} else if (order == 1){
				//XZ
				if (C.y/2){
					if (aX) c_error_pos[C.getFaceQubits(S, 0)] *= -1;
					if (bZ) c_error_pos[C.getFaceQubits(S, 0)] *= -1;
					coord D(C.x,divmod(C.y-1,S.y),C.z,3);
					if (aZ) d_error_pos[D.getFaceQubits(S, 11)] *= -1;
					if (bX) d_error_pos[D.getFaceQubits(S, 11)] *= -1;
				} else {
					if (aX) c_error_pos[C.getFaceQubits(S, 8)] *= -1;
					if (bZ) c_error_pos[C.getFaceQubits(S, 8)] *= -1;
					coord D(C.x,divmod(C.y-1,S.y),divmod(C.z-1,S.z),3);
					if (aZ) d_error_pos[D.getFaceQubits(S, 3)] *= -1;
					if (bX) d_error_pos[D.getFaceQubits(S, 3)] *= -1;
				}
			} else if (order == 2){
				//labels to error operators
				if (C.x/2){
					if (aX) c_error_pos[C.getFaceQubits(S, 7)] *= -1;
					if (bZ) c_error_pos[C.getFaceQubits(S, 7)] *= -1;
					coord D(divmod(C.x-1,S.x),C.y,C.z,3);
					if (aZ) d_error_pos[D.getFaceQubits(S, 8)] *= -1;
					if (bX) d_error_pos[D.getFaceQubits(S, 8)] *= -1;
				} else {
					if (aX) c_error_pos[C.getFaceQubits(S, 6)] *= -1;
					if (bZ) c_error_pos[C.getFaceQubits(S, 6)] *= -1;
					coord D(divmod(C.x-1,S.x),divmod(C.y-1,S.y),C.z,3);
					if (aZ) d_error_pos[D.getFaceQubits(S, 9)] *= -1;
					if (bX) d_error_pos[D.getFaceQubits(S, 9)] *= -1;
				}
			} else if (order == 3){  //XY
				if (C.z/2){
					if (aX){
						c_error_pos[C.getFaceQubits(S, 11)] *= -1;
						c_error_pos[C.getFaceQubits(S, 10)] *= -1;
					}
					if (bZ) c_error_pos[C.getFaceQubits(S, 10)] *= -1;
					coord D(divmod(C.x-1,S.x),C.y,divmod(C.z-1,S.z),3);
					if (aZ) d_error_pos[D.getFaceQubits(S, 8)] *= -1; //4
					if (bX){ 
						d_error_pos[D.getFaceQubits(S, 4)] *= -1;
						d_error_pos[D.getFaceQubits(S, 6)] *= -1;
					}
				} else {
					if (aX){
						c_error_pos[C.getFaceQubits(S, 11)] *= -1;
						c_error_pos[C.getFaceQubits(S, 10)] *= -1;
					}
					if (bZ) c_error_pos[C.getFaceQubits(S, 11)] *= -1;
					coord D(C.x,C.y,divmod(C.z-1,S.z),3);
					if (aZ) d_error_pos[D.getFaceQubits(S, 6)] *= -1;
					if (bX){ 
						d_error_pos[D.getFaceQubits(S, 4)] *= -1;
						d_error_pos[D.getFaceQubits(S, 6)] *= -1;
					}
				}
			} else if (order == 4){  //XZ
				if (C.y/2){
					if (aX){
						c_error_pos[C.getFaceQubits(S, 8)] *= -1;
						c_error_pos[C.getFaceQubits(S, 0)] *= -1;
					}
					if (bZ) c_error_pos[C.getFaceQubits(S, 8)] *= -1;
					coord D(C.x,divmod(C.y-1,S.y),divmod(C.z-1,S.z),3);
					if (aZ) d_error_pos[D.getFaceQubits(S, 3)] *= -1;
					if (bX){ 
						d_error_pos[D.getFaceQubits(S, 3)] *= -1;
						d_error_pos[D.getFaceQubits(S, 11)] *= -1;
					}
				} else {
					if (aX){
						c_error_pos[C.getFaceQubits(S, 8)] *= -1;
						c_error_pos[C.getFaceQubits(S, 0)] *= -1;
					}
					if (bZ) c_error_pos[C.getFaceQubits(S, 0)] *= -1;
					coord D(C.x,divmod(C.y-1,S.y),C.z,3);
					if (aZ) d_error_pos[D.getFaceQubits(S, 11)] *= -1;
					if (bX){ 
						d_error_pos[D.getFaceQubits(S, 3)] *= -1;
						d_error_pos[D.getFaceQubits(S, 11)] *= -1;
					}
				}
			} else if (order == 5){  //YZ
				if (C.x/2){
					if (aX){
						c_error_pos[C.getFaceQubits(S, 6)] *= -1;
						c_error_pos[C.getFaceQubits(S, 7)] *= -1;
					}
					if (bZ) c_error_pos[C.getFaceQubits(S, 6)] *= -1;
					coord D(divmod(C.x-1,S.x),divmod(C.y-1,S.y),C.z,3);
					if (aZ) d_error_pos[D.getFaceQubits(S, 9)] *= -1;
					if (bX){ 
						d_error_pos[D.getFaceQubits(S, 8)] *= -1;
						d_error_pos[D.getFaceQubits(S, 9)] *= -1;
					}
				} else {
					if (aX){
						c_error_pos[C.getFaceQubits(S, 6)] *= -1;
						c_error_pos[C.getFaceQubits(S, 7)] *= -1;
					}
					if (bZ) c_error_pos[C.getFaceQubits(S, 7)] *= -1;
					coord D(divmod(C.x-1,S.x),C.y,C.z,3);
					if (aZ) d_error_pos[D.getFaceQubits(S, 8)] *= -1;
					if (bX){ 
						d_error_pos[D.getFaceQubits(S, 8)] *= -1;
						d_error_pos[D.getFaceQubits(S, 9)] *= -1;
					}
				}
			} else if (order == 6){  //XY
				if (C.z/2){
					if (aX) c_error_pos[C.getFaceQubits(S, 9)] *= -1;
					if (bZ) c_error_pos[C.getFaceQubits(S, 8)] *= -1;
					coord D(C.x,divmod(C.y-1,S.y),divmod(C.z-1,S.z),3);
					if (aZ) d_error_pos[D.getFaceQubits(S, 7)] *= -1;
					if (bX) d_error_pos[D.getFaceQubits(S, 6)] *= -1;
				} else {
					if (aX) c_error_pos[C.getFaceQubits(S, 8)] *= -1;
					if (bZ) c_error_pos[C.getFaceQubits(S, 9)] *= -1;
					coord D(C.x,C.y,divmod(C.z-1,S.z),3);
					if (aZ) d_error_pos[D.getFaceQubits(S, 6)] *= -1;
					if (bX) d_error_pos[D.getFaceQubits(S, 7)] *= -1;
				}
			} else if (order == 7){  //XZ
				if (C.y/2){
					if (aX) c_error_pos[C.getFaceQubits(S, 4)] *= -1;
					if (bZ) c_error_pos[C.getFaceQubits(S, 6)] *= -1;
					coord D(divmod(C.x-1,S.x),divmod(C.y-1,S.y),C.z,3);
					if (aZ) d_error_pos[D.getFaceQubits(S, 10)] *= -1;
					if (bX) d_error_pos[D.getFaceQubits(S, 11)] *= -1;
				} else {
					if (aX) c_error_pos[C.getFaceQubits(S, 6)] *= -1;
					if (bZ) c_error_pos[C.getFaceQubits(S, 4)] *= -1;
					coord  D(C.x,divmod(C.y-1,S.y),C.z,3);
					if (aZ) d_error_pos[D.getFaceQubits(S, 11)] *= -1;
					if (bX) d_error_pos[D.getFaceQubits(S, 10)] *= -1;
				}
			} else if (order == 8){  //YZ
				if (C.x/2){		
					if (aX) c_error_pos[C.getFaceQubits(S, 11)] *= -1;
					if (bZ) c_error_pos[C.getFaceQubits(S, 3)] *= -1;
					coord D(divmod(C.x-1,S.x),C.y,C.z,3);
					if (aZ) d_error_pos[D.getFaceQubits(S, 8)] *= -1;
					if (bX) d_error_pos[D.getFaceQubits(S, 0)] *= -1;
				} else {
					if (aX) c_error_pos[C.getFaceQubits(S, 3)] *= -1;
					if (bZ) c_error_pos[C.getFaceQubits(S, 11)] *= -1;
					coord  D(divmod(C.x-1,S.x),C.y,divmod(C.z-1,S.z),3);
					if (aZ) d_error_pos[D.getFaceQubits(S, 0)] *= -1;
					if (bX) d_error_pos[D.getFaceQubits(S, 8)] *= -1;
				}
			} else if (order == 9){  //XY
				if (C.z/2){
					if (bZ) c_error_pos[C.getFaceQubits(S, 9)] *= -1;
					coord D(C.x,C.y,divmod(C.z-1,S.z),3);
					if (aZ) d_error_pos[D.getFaceQubits(S, 6)] *= -1;
				} else {
					if (bZ) c_error_pos[C.getFaceQubits(S, 8)] *= -1;
					coord D(C.x,divmod(C.y-1,S.y),divmod(C.z-1,S.z),3);
					if (aZ) d_error_pos[D.getFaceQubits(S, 7)] *= -1;
				}
			}  else if (order == 10){  //XZ
				if (C.y/2){
					if (bZ) c_error_pos[C.getFaceQubits(S, 4)] *= -1;
					coord D(C.x,divmod(C.y-1,S.y),C.z,3);
					if (aZ) d_error_pos[D.getFaceQubits(S, 11)] *= -1;
				} else {
					if (bZ) c_error_pos[C.getFaceQubits(S, 6)] *= -1;
					coord D(divmod(C.x-1,S.x),divmod(C.y-1,S.y),C.z,3);
					if (aZ) d_error_pos[D.getFaceQubits(S, 10)] *= -1;
				}
			} else if (order == 11){  //YZ
				if (C.x/2){
					if (bZ) c_error_pos[C.getFaceQubits(S, 11)] *= -1;
					coord D(C.x,C.y,divmod(C.z-1,S.z),3);
					if (aZ) d_error_pos[D.getFaceQubits(S, 0)] *= -1;
				} else {
					if (bZ) c_error_pos[C.getFaceQubits(S, 3)] *= -1;
					coord D(C.x,C.y,C.z,3);
					if (aZ) d_error_pos[D.getFaceQubits(S, 8)] *= -1;
				}
			}
		}
	}
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		coord C(c, S);
		if(this_surf == PLANE && ((C.l == 1 && (C.y == S.y - 1 || C.x == 0)) ||(C.l == 2 && (C.z == S.z - 1 || C.x == 0)))){
			c_error_pos[c] = 1;
			d_error_pos[c] = 1;
		}//	for planar code, remove errors on left, more, and lower boundaries to create edges.
	}
}

void cluster::addLoss(const double & pl, const lossmodel L, const int& seed=0){
	fill(c_loss_pos.begin(), c_loss_pos.end(), 1);
	//Total heuristic Probabilities
	random_device rd;
	mt19937 engine{rd()};
	if (seed != 0) {
		engine.seed(seed+1);
	}
	uniform_real_distribution<> dist(0.0, 1.0);
	
	//Initialization of error operator
	//Error Model 1: uncorrelated error distribution for all physical c_loss_poss
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		if(dist(engine) < LOSSMODELMAP[L](pl).first) {
			c_loss_pos[c] = -1;
		} else {
			c_loss_pos[c] = 1;
		}
	}
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		coord C(c, S);
		if(this_surf == PLANE && ((C.l == 1 && (C.y == S.y - 1 || C.x == 0)) ||(C.l == 2 && (C.z == S.z - 1 || C.x == 0)))){
			c_loss_pos[c] = 1;
		}//	for planar code, remove errors on left, more, and lower boundaries to create edges.
	}
}

int getTaxicabDistance(const coord& S, const coord& c1, const coord& c2, const surfacetype& surf = TORUS){  // compute taxicab distance between two cubes, given their cube numbers
	if (surf == TORUS) {
		return min(abs(c1.x - c2.x), S.x-abs(c1.x - c2.x)) + min(abs(c1.y - c2.y), S.y-abs(c1.y - c2.y)) + min(abs(c1.z - c2.z), S.z-abs(c1.z - c2.z));
	} else if (surf == PLANE){
		return abs(c1.x - c2.x) + abs(c1.y - c2.y) + abs(c1.z - c2.z);
	} else{
		return 0;
	}
}
//
coord getTaxicabDisplacement(const coord& S, const coord& c1, const coord& c2, const surfacetype& surf = TORUS){  // compute taxicab distance between two points with given coordinates

	if (surf == TORUS) {
		int x_dist, y_dist, z_dist, x_relative_pos, y_relative_pos, z_relative_pos;
		x_dist = min(S.x-abs(c1.x - c2.x), abs(c1.x - c2.x));
		y_dist = min(S.y-abs(c1.y - c2.y), abs(c1.y - c2.y));
		z_dist = min(S.z-abs(c1.z - c2.z), abs(c1.z - c2.z));

		//replace this part with displacement function
		if (x_dist == 0) { //two vertices on same x
			x_relative_pos = 0;

		} else if ((x_dist == abs(c1.x - c2.x) && c1.x < c2.x) ||
				   (x_dist == S.x - abs(c1.x - c2.x) && c1.x > c2.x)){
			x_relative_pos = 1;
		} else {
			x_relative_pos = -1;
		}

		if (y_dist == 0) { //two vertices on same y
			y_relative_pos = 0;
		} else if ((y_dist == abs(c1.y - c2.y) && c1.y < c2.y) ||
				   (y_dist == S.y - abs(c1.y - c2.y) && c1.y > c2.y)){
			y_relative_pos = 1;

		} else {
			y_relative_pos = -1;
		}

		if (z_dist == 0) { //two vertices on same z
			z_relative_pos = 0;
		} else if ((z_dist == abs(c1.z - c2.z) && c1.z < c2.z) ||
				   (z_dist == S.z - abs(c1.z - c2.z) && c1.z > c2.z)){
			z_relative_pos = 1;

		} else {
			z_relative_pos = -1;
		}

		return {x_relative_pos * x_dist, y_relative_pos * y_dist, z_relative_pos * z_dist, 0};

	} else if (surf == PLANE){
		return {c2.x - c1.x, c2.y - c1.y, c2.z - c1.z, 0};
	} else{
		return {0, 0, 0, 0};
	}
}

int getTaxicabDistance(const coord& S, const vector<int>& chunk1, const vector<int>& chunk2, const surfacetype& surf = TORUS){  // compute taxicab distance between two coords
	int min_dist = INT_MAX;
	for (int i = 0; i < chunk1.size(); i++) {
		for (int j = 0; j < chunk2.size(); j++) {
			int dist = getTaxicabDistance(S, coord(chunk1[i],S), coord(chunk2[j],S), surf);
			min_dist = (dist < min_dist) ? dist : min_dist;
		}
	}
	return min_dist;
}


void cluster::getXMeasurements(){
	for (int c = 0; c < S.x*S.y*S.z; c++){
		stabs[c] = 1;
		//measurement of vertex operator
		vertex avertex(c, S);
		for (int pos = 0; pos < 6; pos ++) {
			if (c_error_pos[avertex.physical_qubits[pos]] == -1) {
				stabs[c] *= -1;
			}
		}
	}
}

void cluster::getZMeasurements(){
	for (int c = 0; c < S.x*S.y*S.z; c++){
		stabs[c + S.x*S.y*S.z] = 1;
		//measurement of vertex operator
		vertex avertex(c, S);
		for (int pos = 0; pos < 6; pos ++) {
			if (d_error_pos[avertex.physical_qubits[pos]] == -1) {
				stabs[c + S.x*S.y*S.z] *= -1;
			}
		}
	}
}

int cluster::getLeftDistance(const coord& S, const vector<int>& chunk){
	int min_dist = INT_MAX;
	for (int i = 0; i < chunk.size(); i++) {
		int dist = coord(chunk[i],S).x;
		min_dist = (dist < min_dist) ? dist : min_dist;
	}
	for (int i = 0; i < left_boundary_chunks.size(); i++) {
		int dist = getTaxicabDistance(S, super_chunks[left_boundary_chunks[i]], chunk);
		min_dist = (dist < min_dist) ? dist : min_dist;
	}
	return min_dist;
}

int cluster::getRightDistance(const coord& S, const vector<int>& chunk){
	int min_dist = INT_MAX;
	for (int i = 0; i < chunk.size(); i++) {
		int dist = S.x-coord(chunk[i],S).x;
		min_dist = (dist < min_dist) ? dist : min_dist;
	}
	for (int i = 0; i < right_boundary_chunks.size(); i++) {
		int dist = getTaxicabDistance(S, super_chunks[right_boundary_chunks[i]], chunk);
		min_dist = (dist < min_dist) ? dist : min_dist;
	}
	return min_dist;
}

int cluster::findParent(int cV){
	for (int i = 0; i < super_chunks.size(); i++){
		if(find(super_chunks[i].begin(), super_chunks[i].end(), cV) != super_chunks[i].end()){
			return i;
		}
	}
	return -1;
}

void cluster::getSuperChunks(){
	super_chunks = {}; // initialization
	left_boundary_chunks = {};// initialization
	right_boundary_chunks = {};// initialization
	boundary_info = {}; //initialization
	for(int cV = 0; cV < S.x*S.y*S.z; cV++) {
		if(coord(cV, S).x == 0){///for PLANAR code only: only take the vertices in the PLANE
			continue;
		}
		boundary_info.push_back(0);
		super_chunks.push_back(vector<int>(1,cV));
	}
	//merge chunks
	for(int c = 0; c < 3*S.x*S.y*S.z; c++){
		if (c_loss_pos[c] < 0) {//qubit loss encountered
			coord C(c,S);
			if(C.l == 0 && C.x == 0){//on left face (boundary)
				int cV = coord(divmod(C.x+1, S.x), C.y, C.z, 0).hash(S);
				int i = findParent(cV);
				boundary_info[i] = 1;//label marks boundary node
				left_boundary_chunks.push_back(i);
			} else if(C.l == 0 && C.x == S.x-1){// on right face (boundary)
				int cV = coord(C.x, C.y, C.z, 0).hash(S);
				int i = findParent(cV);
				boundary_info[i] = 1;//label marks boundary node
				right_boundary_chunks.push_back(i);
			} else{
				//get adjacent vertices
				int cV1 = coord(C.x, C.y, C.z, 0).hash(S);
				int cV2;
				if (coord(c,S).l == 0) {
					cV2 = coord(divmod(C.x+1, S.x), C.y, C.z, 0).hash(S);
				} else if (coord(c,S).l == 1) {
					cV2 = coord(C.x, divmod(C.y+1, S.y), C.z, 0).hash(S);
				} else{
					cV2 = coord(C.x, C.y, divmod(C.z+1, S.z), 0).hash(S);
				}
				int i = findParent(cV1);
				if(find(super_chunks[i].begin(), super_chunks[i].end(), cV2) == super_chunks[i].end()){ //check if two vertices not in the same chunk
					int j = findParent(cV2);//combine two chunks
					
					if ((find(right_boundary_chunks.begin(), right_boundary_chunks.end(), i) != right_boundary_chunks.end() && find(left_boundary_chunks.begin(), left_boundary_chunks.end(), j) != left_boundary_chunks.end())||
						(find(right_boundary_chunks.begin(), right_boundary_chunks.end(), j) != right_boundary_chunks.end()&& find(left_boundary_chunks.begin(), left_boundary_chunks.end(), i) != left_boundary_chunks.end())) {//failed by percolation
						throw(0);
					} else if (find(right_boundary_chunks.begin(), right_boundary_chunks.end(), i) != right_boundary_chunks.end() || find(left_boundary_chunks.begin(), left_boundary_chunks.end(), i) != left_boundary_chunks.end()) {//i is at boundary
						for(int k = 0; k < super_chunks[j].size(); k++){//add to chunk i
							super_chunks[i].push_back(super_chunks[j][k]);
						}
						super_chunks[j]={};//clear chunk j
					} else {//j is at boundary or neither at boundary
						for(int k = 0; k < super_chunks[i].size(); k++){//add to chunk j
							super_chunks[j].push_back(super_chunks[i][k]);
						}
						super_chunks[i]={};//clear chunk i
					}
				}
			}
		}
	}
}

void cluster::getSurf(){ //debug
	for (int i =0; i < S.x; i++) {
		for (int j =0; j < S.y; j++) {
			for (int l=0; l < 2; l++) {
				int bit = 1;
				for (int k =0; k < S.z; k ++) {
					bit *= c_error_pos[coord(i, j, k, l).hash(S)];
				}
				surf[coord(i, j, 0, l).hash(S)] = bit;
			}
		}
	}
}

void cluster::printSurf(){
	cout << "surf: " <<endl;
	vector<vector<string>> print_out;
	for(int i = 0; i < 2*S.x; i++){
		vector<string> a_line;
		for(int j = 0; j < 2*S.y; j++){
			a_line.push_back(" ");
		}
		print_out.push_back(a_line);
	}
	for (int x = 0; x < S.x; x++) {
		for (int y = 0; y < S.y; y++) {
			if (surf[coord(x, y, 0, 0).hash(S)] < 0) {
				print_out[2 * y][2 * x + 1] = "\033[1;96mZ\033[0m";
			} 
			if (surf[coord(x, y, 0, 1).hash(S)] < 0){
				print_out[2 * y +1][2 * x] = "\033[1;96mZ\033[0m";
			}
		}
	}
	cout << "surface:" << endl;
	printMatrix(print_out);
}

void cluster::printQubit(surfacetype s = TORUS){
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
					if (d_error_pos[coord(x, y, divmod(k/2 - 1, S.z), 2).hash(S)] < 0){
						print_out[2 * y + 1][2*x + 1] = "\033[1;33mZ\033[0m";
					} else if (d_error_pos[coord(x, y, divmod(k/2 - 1, S.z), 2).hash(S)] == 2){
						print_out[2 * y + 1][2*x + 1]  = "\033[1;33mC\033[0m";
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
					if (d_error_pos[coord(x, y, k/2, 0).hash(S)] < 0){
						print_out[2*y+1][divmod(2*x+2,2*S.x)] = "\033[1;33mZ\033[0m";
					} else if (d_error_pos[coord(x, y, k/2, 0).hash(S)] == 2){
						print_out[2*y+1][divmod(2*x+2,2*S.x)]  = "\033[1;33mC\033[0m";
					}
					if (d_error_pos[coord(x, y, k/2, 1).hash(S)] < 0){
						print_out[divmod(2*y+2,2*S.y)][2*x+1] = "\033[1;33mZ\033[0m";
					} else if (d_error_pos[coord(x, y, k/2, 1).hash(S)] == 2){
						print_out[divmod(2*y+2,2*S.y)][2*x+1]  = "\033[1;33mC\033[0m";
					}
					if (c_error_pos[coord(x, y, k/2, 2).hash(S)] < 0){
						print_out[2 * y][2 * x] = "\033[1;96mZ\033[0m";
					} else if (c_error_pos[coord(x, y, k/2, 2).hash(S)] == 2){
						print_out[2 * y][2 * x]   = "\033[1;33mC\033[0m";
					}
					if (s == PLANE && y == S.y-1){
					} else if (stabs[coord(x, y, k/2, 0).hash(S) + S.x*S.y*S.z] > 2){
						print_out[2*y+1][2*x+1] = to_string(stabs[coord(x, y, k/2, 0).hash(S) + S.x*S.y*S.z]);
					} else if (stabs[coord(x, y, k/2, 0).hash(S) + S.x*S.y*S.z] > 0){
						print_out[2*y+1][2*x+1] = "◯";
					} else {
						print_out[2*y+1][2*x+1] = "\033[1;92m⊕\033[0m";
					}		
				}
			}
		}
		cout << "layer:" << k << endl;
		printMatrix(print_out);
	}
	
}

void cluster::printSuperChunks(){
	for(int k = 0; k < 2*S.z; k++){
		vector<vector<string>> print_out;
		for(int i = 0; i < 2*S.x; i++){
			vector<string> a_line;
			for(int j = 0; j < 2*S.y; j++){
				a_line.push_back(" ");
			}
			print_out.push_back(a_line);
		}
		if (k % 2 == 0) {//even layer
			for (int x = 0; x < S.x; x++) {
				for (int y = 0; y < S.y; y++) {
					int c = coord(x, y, k/2, 0).hash(S);
					stringstream s1, s2, s3;
					
					int parent= findParent(c);
					if (super_chunks[parent].size() != 1 && parent != -1){
						print_out[2 * y][2 * x] = colorvector[divmod(parent,colorvector.size())];//looping through colors
					}

					if (c_loss_pos[coord(x, y, k/2, 0).hash(S)] < 0){
						print_out[2 * y][2 * x + 1] = "\033[1;96mL\033[0m";
					} else if (c_loss_pos[coord(x, y, k/2, 0).hash(S)] == 2){
						print_out[2 * y][2 * x + 1]  = "\033[1;33mC\033[0m";
					}
					if (c_loss_pos[coord(x, y, k/2, 1).hash(S)] < 0){
						print_out[2 * y + 1][2 * x] = "\033[1;96mL\033[0m";
					} else if (c_loss_pos[coord(x, y, k/2, 1).hash(S)] == 2){
						print_out[2 * y + 1][2 * x]  = "\033[1;33mC\033[0m";
					}
				}
			}
		} else {//odd layer
			for (int x = 0; x < S.x; x++) {
				for (int y = 0; y < S.y; y++) {
					stringstream s1;
					if (c_loss_pos[coord(x, y, k/2, 2).hash(S)] < 0){
						print_out[2 * y][2 * x] = "\033[1;96mL\033[0m";
					} else if (c_loss_pos[coord(x, y, k/2, 2).hash(S)] == 2){
						print_out[2 * y][2 * x]   = "\033[1;33mC\033[0m";
					}
				}
			}
		}
		cout << "layer:" << k << endl;
		printMatrix(print_out);
	}
}


vector<float> cluster::getWindow(coord& C, const bool& imprint = 0){
	if (imprint){
		if (C.l <= 2){ //z errors
			c_error_pos[C.hash(S)] = 2;
		} else {
			d_error_pos[C.hash(S) - 3*S.x*S.y*S.z] = 2;
		}
	}
	
	vector<int> window_int((Sw-1)*Sw*(2*Sw-1),0);
	if (C.l == 0) { //0 qubits
		for (int dx = 0; dx < Sw; dx++) {
			for (int dy = 0; dy < Sw - 1; dy++) {
				for (int dz = 0; dz < 2*Sw - 1; dz++) {
					if (dz % 2 == 0){
						window_int[dz/2 * (Sw - 1)*Sw + dy*Sw + dx] = -(stabs[coord(divmod(C.x+dy-Sw/2+1,S.x),divmod(C.y+dx-Sw/2,S.y),divmod(C.z+dz/2-Sw/2,S.z),0).hash(S)]-1)/2; //green
						if (imprint){//debug
							stabs[coord(divmod(C.x+dy-Sw/2+1,S.x),divmod(C.y+dx-Sw/2,S.y),divmod(C.z+dz/2-Sw/2,S.z),0).hash(S)] = 2;
						}	
					} else {
						window_int[dz/2 * (Sw - 1)*Sw + dy*Sw + dx + (Sw-1)*Sw*Sw] = (stabs[coord(divmod(C.x+dx-Sw/2,S.x),divmod(C.y+dy-Sw/2,S.y),divmod(C.z+dz/2-Sw/2-1,S.z),1).hash(S)]-1)/2; //purple
						if (imprint){//debug
							stabs[coord(divmod(C.x+dx-Sw/2,S.x),divmod(C.y+dy-Sw/2,S.y),divmod(C.z+dz/2-Sw/2,S.z),1).hash(S)] = 2;
						}
					}
				}
			}
		}
	} else if (C.l == 1) { //1 qubits
		for (int dx = 0; dx < Sw; dx++) {
			for (int dy = 0; dy < Sw - 1; dy++) {
				for (int dz = 0; dz < 2*Sw - 1; dz++) {
					if (dz % 2 == 0){
						window_int[dz/2 * (Sw - 1)*Sw + dy*Sw + dx] = (stabs[coord(divmod(C.x+dx-Sw/2,S.x),divmod(C.y+dy-Sw/2+1,S.y),divmod(C.z+dz/2-Sw/2,S.z),0).hash(S)]-1)/2; //green
						if (imprint){//debug
							stabs[coord(divmod(C.x+dx-Sw/2,S.x),divmod(C.y+dy-Sw/2+1,S.y),divmod(C.z+dz/2-Sw/2,S.z),0).hash(S)] = 3;
						}	
					} else {
						window_int[dz/2 * (Sw - 1)*Sw + dy*Sw + dx + (Sw-1)*Sw*Sw] = -(stabs[coord(divmod(C.x+dy-Sw/2,S.x),divmod(C.y+dx-Sw/2,S.y),divmod(C.z+dz/2-Sw/2,S.z),1).hash(S)]-1)/2; //purple
						if (imprint){//debug
							stabs[coord(divmod(C.x+dy-Sw/2,S.x),divmod(C.y+dx-Sw/2,S.y),divmod(C.z+dz/2-Sw/2,S.z),1).hash(S)] = 3;
						}
					}
				}
			}
		}
	} else if (C.l == 2) { //2 qubits
		for (int dx = 0; dx < Sw; dx++) {
			for (int dy = 0; dy < Sw - 1; dy++) {
				for (int dz = 0; dz < 2*Sw - 1; dz++) {
					if (dz % 2 == 0){
						window_int[dz/2 * (Sw-1)*Sw + dy*Sw + dx] = 3*(stabs[coord(divmod(C.x+dx-Sw/2,S.y),divmod(C.y+dz/2-Sw/2,S.x),divmod(C.z+dy-Sw/2+1,S.z),0).hash(S)]-1)/2; //green
						if (imprint){//debug
							stabs[coord(divmod(C.x+dx-Sw/2,S.y),divmod(C.y+dz/2-Sw/2,S.x),divmod(C.z+dy-Sw/2+1,S.z),0).hash(S)] = 4;
						}	
					} else {
						window_int[dz/2 * (Sw - 1)*Sw + dy*Sw + dx + (Sw-1)*Sw*Sw] = 3*(stabs[coord(divmod(C.x+dy-Sw/2,S.x),divmod(C.y+dz/2-Sw/2,S.y),divmod(C.z+dx-Sw/2,S.z),1).hash(S)]-1)/2; //purple
						if (imprint){//debug
							stabs[coord(divmod(C.x+dy-Sw/2,S.x),divmod(C.y+dz/2-Sw/2,S.y),divmod(C.z+dx-Sw/2,S.z),1).hash(S)] = 4;
						}
					}
				}
			}
		}
	} if (C.l == 3) { //3 qubits
		for (int dx = 0; dx < Sw; dx++) {
			for (int dy = 0; dy < Sw - 1; dy++) {
				for (int dz = 0; dz < 2*Sw - 1; dz++) {
					if (dz % 2 == 0){
						window_int[dz/2 * (Sw - 1)*Sw + dy*Sw + dx] = (stabs[coord(divmod(C.x+dy-Sw/2+1,S.x),divmod(C.y+dx-Sw/2,S.y),divmod(C.z+dz/2-Sw/2,S.z),1).hash(S)]-1)/2; //purple
						if (imprint){//debug
							stabs[coord(divmod(C.x+dy-Sw/2+1,S.x),divmod(C.y+dx-Sw/2,S.y),divmod(C.z+dz/2-Sw/2,S.z),1).hash(S)] = 3;
						}	
					} else {
						window_int[dz/2 * (Sw - 1)*Sw + dy*Sw + dx + (Sw-1)*Sw*Sw] = (stabs[coord(divmod(C.x+dx-Sw/2+1,S.x),divmod(C.y+dy-Sw/2+1,S.y),divmod(C.z+dz/2-Sw/2+1,S.z),0).hash(S)]-1)/2; //purple
						if (imprint){//debug
							stabs[coord(divmod(C.x+dx-Sw/2+1,S.x),divmod(C.y+dy-Sw/2+1,S.y),divmod(C.z+dz/2-Sw/2+1,S.z),0).hash(S)] = 3;
						}
					}
				}
			}
		}
	} else if (C.l == 4) { //4 qubits
		for (int dx = 0; dx < Sw; dx++) {
			for (int dy = 0; dy < Sw - 1; dy++) {
				for (int dz = 0; dz < 2*Sw - 1; dz++) {
					if (dz % 2 == 0){
						window_int[dz/2 * (Sw - 1)*Sw + dy*Sw + dx] = (stabs[coord(divmod(C.x+dx-Sw/2,S.x),divmod(C.y+dy-Sw/2+1,S.y),divmod(C.z+dz/2-Sw/2,S.z),1).hash(S)]-1)/2; //purple
						if (imprint){//debug
							stabs[coord(divmod(C.x+dx-Sw/2,S.x),divmod(C.y+dy-Sw/2+1,S.y),divmod(C.z+dz/2-Sw/2,S.z),1).hash(S)] = 3;
						}	
					} else {
						window_int[dz/2 * (Sw - 1)*Sw + dy*Sw + dx + (Sw-1)*Sw*Sw] = (stabs[coord(divmod(C.x+dy-Sw/2+1,S.x),divmod(C.y+dx-Sw/2+1,S.y),divmod(C.z+dz/2-Sw/2+1,S.z),0).hash(S)]-1)/2; //green
						if (imprint){//debug
							stabs[coord(divmod(C.x+dy-Sw/2+1,S.x),divmod(C.y+dx-Sw/2+1,S.y),divmod(C.z+dz/2-Sw/2+1,S.z),0).hash(S)] = 3;
						}
					}
				}
			}
		}
	} else if (C.l == 5) { //5 qubits
		for (int dx = 0; dx < Sw; dx++) {
			for (int dy = 0; dy < Sw - 1; dy++) {
				for (int dz = 0; dz < 2*Sw - 1; dz++) {
					if (dz % 2 == 0){
						window_int[dz/2 * (Sw-1)*Sw + dy*Sw + dx] = (stabs[coord(divmod(C.x+dx-Sw/2,S.y),divmod(C.y+dz/2-Sw/2,S.x),divmod(C.z+dy-Sw/2+1,S.z),1).hash(S)]-1)/2; //purple
						if (imprint){//debug
							stabs[coord(divmod(C.x+dx-Sw/2,S.y),divmod(C.y+dz/2-Sw/2,S.x),divmod(C.z+dy-Sw/2+1,S.z),1).hash(S)] = 3;
						}	
					} else {
						window_int[dz/2 * (Sw - 1)*Sw + dy*Sw + dx + (Sw-1)*Sw*Sw] = (stabs[coord(divmod(C.x+dy-Sw/2+1,S.x),divmod(C.y+dz/2-Sw/2+1,S.y),divmod(C.z+dx-Sw/2+1,S.z),0).hash(S)]-1)/2; //green
						if (imprint){//debug
							stabs[coord(divmod(C.x+dy-Sw/2+1,S.x),divmod(C.y+dz/2-Sw/2+1,S.y),divmod(C.z+dx-Sw/2+1,S.z),0).hash(S)] = 3;
						}
					}
				}
			}
		}
	}
	vector<float> window(window_int.begin(), window_int.end());
	// window.push_back(float(C.l));
	return window;
}

void cluster::surfaceCorrect(PerfectMatching* pm, const int& vertices_num, const vector<int>& vertexPosition, const vector<int>& boundary_nodes, const int& verbosity = 0){//debug
	//Find the error operator of each pair.
	vector<int> error_op_Z_pos; //correction operator not reduced
	vector<int> matchPosition; //matching vertex of vertexPosition respectively
	for (int i = 0; i < vertices_num; i++) {
		int ac = vertexPosition[i];
		coord relative_pos; // relative position between two vertices in a pair
		//get relative position (vector) of pair to determine the path in between.
		vertex avertex(ac, S), bvertex;
		if (pm->GetMatch(i) < vertices_num){
			int bc = vertexPosition[pm->GetMatch(i)];// position of vertexa's match, vertexb
			matchPosition.push_back(bc);
			if (count(matchPosition.begin(), matchPosition.end(), ac)) continue; //Prevent recounting of vertices
			bvertex = vertex(bc, S);
			relative_pos = getTaxicabDisplacement(S, coord(ac, S), coord(bc, S), this_surf);
		}else if (this_surf == PLANE){ // matched to boundary
			matchPosition.push_back(0);
			relative_pos = {boundary_nodes[i], 0, 0, 0};
		}
		int n;
		if (relative_pos.x > 0) {//a to the left of b
			n = avertex.physical_qubits[1];//use right c_error_pos
			for(int i =0; i < abs(relative_pos.x); i++){
				error_op_Z_pos.push_back(coord(divmod(coord(n,S).x + i, S.x), coord(n,S).y, 0, 0).hash(S));
			}
		} else if(relative_pos.x < 0) {//a to the right of b
			n = avertex.physical_qubits[0];//use left c_error_pos
			for(int i =0; i < abs(relative_pos.x); i++){
				error_op_Z_pos.push_back(coord(divmod(coord(n,S).x - i, S.x), coord(n,S).y, 0, 0).hash(S));
			}
		}
		if (relative_pos.y > 0) {//a above b
			n = bvertex.physical_qubits[2];//use upper c_error_pos
			for(int i =0; i < abs(relative_pos.y); i++){
				error_op_Z_pos.push_back(coord(coord(n,S).x, divmod(coord(n,S).y - i, S.y), 0, 1).hash(S));
			}
		} else if (relative_pos.y < 0) {//a below b
			n = bvertex.physical_qubits[3];//use lower c_error_pos
			for(int i =0; i < abs(relative_pos.y); i++){
				error_op_Z_pos.push_back(coord(coord(n,S).x, divmod(coord(n,S).y + i, S.y), 0, 1).hash(S));
			}
		}

	}

	for (int i = 0; i < error_op_Z_pos.size(); i++) {
		surf[error_op_Z_pos[i]] *= -1; //act local Z operator
	}
	if(verbosity == 3){
		cout << "vertex/match/bn: "<< endl;
		//printout test
		for (int i = 0 ; i < vertices_num; i++) {
			cout << coord(vertexPosition[i],S)<< "/";
			cout << coord(matchPosition[i],S)<<"/";
			if (this_surf == PLANE){
				cout << boundary_nodes[i] << "/";
			}	
			cout << getTaxicabDistance(S,coord(vertexPosition[i],S),coord(matchPosition[i],S),PLANE) << endl;
		}
		cout << "errorOps: " << endl;
		for (int i = 0 ; i <error_op_Z_pos.size(); i++) {
			cout << coord(error_op_Z_pos[i],S) << " ";
		}
		cout<<endl;
	}
}

void cluster::decodeWithNN(cppflow::model model, int verbosity = 0){
	vector<int> check_pos;
	for (int c = 0; c < S.x*S.y*S.z; c++) {
		if(stabs[c] < 0){
			vertex aVertex(c,S);
			for (int i = 0; i < 6; i ++) {
				int pos = aVertex.physical_qubits[i];
				if (count(check_pos.begin(), check_pos.end(), pos)) continue;
				check_pos.push_back(pos);
			}
		}
	}
	// for (int c = S.x*S.y*S.z; c < 2*S.x*S.y*S.z; c++) {
	// 	if(stabs[c] < 0){
	// 		vertex aVertex(c-S.x*S.y*S.z,S);
	// 		for (int i = 0; i < 6; i ++) {
	// 			int pos = aVertex.physical_qubits[i] + 3*S.x*S.y*S.z;
	// 			if (count(check_pos.begin(), check_pos.end(), pos)) continue;
	// 			check_pos.push_back(pos);
	// 		}
	// 	}
	// }

	for (int i = 0; i < check_pos.size(); i++){
		int pos = check_pos[i];
		coord C(pos,S);
		// cout << C << endl;
		vector<float> window = getWindow(C);
		
		
		auto input = cppflow::tensor(window, {1,(Sw-1)*Sw*(2*Sw-1)});
		auto output = model(input);
		
		// cout << output <<endl;
		// if (C.x == 8 && C.y == 7 && C.z == 2){
		// 	cout << window <<endl;
		// }
		if (pos < 3*S.x*S.y*S.z){
			c_error_pos[pos] *= -(2*int(round(output.get_data<float>()[0]))-1);
		}
		// } else {
		// 	d_error_pos[pos-3*S.x*S.y*S.z] *= -(2*int(round(output.get_data<float>()[0]))-1);
		// }
		
	}
	// cout <<"testing part" <<endl;
	// vector<float> test_window = {0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1, -1, 0, 0, -1, -1, 0, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1};
	// auto input = cppflow::tensor(test_window, {1,(Sw-1)*Sw*(2*Sw-1)+1});
	// auto output = model(input);
	// cout << output << endl;
}

int cluster::decodeWithMWPM(int verbosity = 0, bool dir = 0, bool make_corrections = 0){
	//PAIR MATCHING
    vector<int> vertexPosition;///actual (non boundary) vertices to be matched
	
    for (int c = 0; c < S.x*S.y*S.z; c++) {
		int c_prime = c;
		if (dir){
			c_prime +=  S.x*S.y*S.z;
		}
        if (stabs[c_prime] == -1) {
			if(this_surf == PLANE && coord(c, S).x == 0){///for PLANAR code only: only take the vertices in the PLANE
				continue;
			}
            vertexPosition.push_back(c);
        }
    }


	int vertices_num = vertexPosition.size(), matches_num, edges_num;

	if (this_surf == TORUS) {
		matches_num = vertices_num;
		edges_num = vertices_num * (vertices_num - 1)/2;
	} else if (this_surf == PLANE){
		matches_num = 2 * vertices_num;
		edges_num = vertices_num * (vertices_num - 1) + vertices_num;
	}
    PerfectMatching *pm = new PerfectMatching(matches_num, edges_num);
	pm->options.verbose = false;
	vector<int> boundary_nodes; //vector to keep track of x distance to boundary nodes

	if (this_surf == PLANE) {//add in the boundary nodes
		for (int i = 0; i < vertices_num; i++) {
			int cx = vertexPosition[i] % (S.x*S.y*S.z) % (S.x*S.y) % S.x;
			if (cx * 2 <= S.x) { // use left boundary
				boundary_nodes.push_back(-cx);
				pm->AddEdge(i, i + vertices_num, cx);
			} else{ //use right boundary
				boundary_nodes.push_back(S.x - cx);
				pm->AddEdge(i, i + vertices_num, S.x - cx);
			}
		}
	}

	for (int i = 0; i < vertices_num; i++) {
		coord c1(vertexPosition[i], S);
        for (int j = i + 1; j < vertices_num; j++) {// add in interconnections
            coord c2(vertexPosition[j],S);
			int dist = getTaxicabDistance(S, c1, c2, this_surf);
			if (this_surf == TORUS){
				pm->AddEdge(i,j,dist);
			} else if (this_surf == PLANE && dist < abs(boundary_nodes[i]) + abs(boundary_nodes[j])) {//optimization
				pm->AddEdge(i,j,dist);
				pm->AddEdge(i+vertices_num,j+vertices_num,0);//boundary interconnection			}
			}
        }
    }
	//solve the graph using MWPM decoder
	pm->Solve();
	//surface correction
	if (make_corrections) {
		surfaceCorrect(pm, vertices_num, vertexPosition, boundary_nodes, verbosity);///debug
	}

	//alternative check (quick method)
	int parity = 1;
	
	vector<int> matchPosition; //matching vertex of vertexPosition respectively
	for (int i = 0; i < vertices_num; i++) {
		if (this_surf == PLANE){
			if (pm->GetMatch(i) == i + vertices_num && boundary_nodes[i] > 0) { //matched to right boundary
				parity *= -1;
			}
		} else if (this_surf == TORUS){		
			int aC = vertexPosition[i];
			int bC = vertexPosition[pm->GetMatch(i)];// position of vertexa's match, vertexb

			matchPosition.push_back(bC);
			int x1 = aC % S.x;
			int x2 = bC % S.x;
			if (!((2*abs(x1 - x2) > S.x && x1 < x2))){ //crossed boundary
				parity *= -1;
			}
		}
	}

	for (int i = 0; i < S.y; i++){//errors on right boundary
		for (int j = 0; j < S.z; j++){
			if (!dir){
				parity *= c_error_pos[coord(S.x-1, i, j, 0).hash(S)];
			} else {
				parity *= d_error_pos[coord(S.x-1, i, j, 0).hash(S)];
			}	
		}
	}
	
	delete pm;
	return parity;
}

int cluster::decodeWithMWPMloss(int verbose = 0, bool dir = 0, bool make_corrections = 0){
    vector<int> wrong_chunks;//actual (non boundary) chunks to be matched
	//e.g. [1,3,5,7]
	//get wrong_chunks
	for (int i = 0; i < super_chunks.size(); i++) {
		if (boundary_info[i] != 0) {// remove boundary chunks
			continue;
		}
		int measurement = 1;
		for (int j = 0; j < super_chunks[i].size(); j++) {
			measurement *= stabs[super_chunks[i][j]];
		}
        if (measurement == -1) {
            wrong_chunks.push_back(i);
        }
    }
	
	
//	cout << wrong_chunks <<endl; /// debug
 	int vertices_num = wrong_chunks.size(), matches_num, edges_num;
	
	matches_num = 2 * vertices_num;
	edges_num = vertices_num * (vertices_num - 1) + vertices_num;

    PerfectMatching *pm = new PerfectMatching(matches_num, edges_num);
	pm->options.verbose = false;

	//add vertex-boundary connections;
    for (int i = 0; i < vertices_num; i++) {
		int dl = getLeftDistance(S, super_chunks[wrong_chunks[i]]);
		int dr = getRightDistance(S, super_chunks[wrong_chunks[i]]);
		int dmin1 = dl <= dr ? dl : dr;
		boundary_info[wrong_chunks[i]] = dl <= dr ? -1 : 1;
//		cout << wrong_chunks[i] << ":b:" << dmin1 << endl; /// debug
		pm->AddEdge(i, i + vertices_num, dmin1);
		//add vertex-vertex; boundary-boundary connections;
		for (int j = i + 1; j < vertices_num; j++) {// add in interconnections
			int dl = getLeftDistance(S, super_chunks[wrong_chunks[j]]);
			int dr = getRightDistance(S, super_chunks[wrong_chunks[j]]);
			int dmin2 = dl <= dr ? dl : dr;
			int dist = getTaxicabDistance(S, super_chunks[wrong_chunks[i]], super_chunks[wrong_chunks[j]]);
//			cout << wrong_chunks[i] << ":" << wrong_chunks[j]  << dist << endl; /// debug
			if (dist < dmin1 + dmin2) {//optimization
				pm->AddEdge(i,j,dist);
				pm->AddEdge(i + vertices_num, j + vertices_num, 0);//boundary interconnection
				
			}
		}
	}
	
	//solve the graph using MWPM decoder
	pm->Solve();

//	//surface correction
//	surfaceCorrect(pm, vertices_num, wrong_chunks, boundary_info, verbose);
	
	int parity = 1;
	for (int i = 0; i < vertices_num; i++) {
		if (pm->GetMatch(i) == i + vertices_num && boundary_info[wrong_chunks[i]] < 0) { //matched to left boundary
			parity *= -1;
		}
	}
	for (int i = 0; i < S.x; i++){//errors on left plane
		for (int j = 0; j < S.z; j++){
			parity *= c_error_pos[coord(0, i, j, 0).hash(S)];
		}
	}
	for (int i = 0; i < left_boundary_chunks.size(); i++) {
		for(int cV = 0; cV < super_chunks[left_boundary_chunks[i]].size(); cV++){
			parity *= stabs[super_chunks[left_boundary_chunks[i]][cV]];
		}
	}

	delete pm;
	return parity;
}


void testDecoding(cppflow::model model, cluster& testcluster, const double& p, const double& q, const int& seed, surfacetype surf, noisemodel N, lossmodel L, int verbosity = 0, bool decode_w_NN =1){
	start_t = clock();

	//noise model secsion
	if (N == OFF);
	else if (N == GATE) {
		testcluster.addGateNoise(p, N, seed);
	} else if (N == MANUAL){
		testcluster.addError();
	} else {
		testcluster.addNoise(p, 0, N, seed);
	}

	//loss decoding section
	if (L != OFF_loss){ 
		testcluster.addLoss(q, L);
		try {
			testcluster.getSuperChunks();
		} catch (...) {
			cout << "percolation failure!" << endl;
			return;
		}
		if (verbosity >=1) {
			testcluster.printSuperChunks();
		}
	}

	//printing and visualization
	testcluster.getXMeasurements();
	testcluster.getZMeasurements();

	if (verbosity >= 1){
		testcluster.getSurf();///debug
		testcluster.printSurf();///debug
		if (verbosity == 2) {
			testcluster.printQubit(surf);
		}
		cout << "t(Generation):" << double(clock()-start_t)/CLOCKS_PER_SEC << endl;//timing
		start_t = clock();
	}
	
	//NN decoding section
	if (decode_w_NN && surf != PLANE){
		testcluster.decodeWithNN(model);
		testcluster.getXMeasurements();
		testcluster.getZMeasurements();
		if (verbosity >= 1){
			testcluster.getSurf();///debug
			testcluster.printSurf();///debug
			if (verbosity == 2) {
				testcluster.printQubit(surf);
			}
			cout << "t(decodeWithNN):" << double(clock()-start_t)/CLOCKS_PER_SEC << endl;//timing
			start_t = clock();
		}
	}
	int parity;
	if (L == OFF_loss){
		parity = testcluster.decodeWithMWPM(verbosity, 0, 1);
	} else {
		parity = testcluster.decodeWithMWPMloss(verbosity, 0, 1);
	}

	testcluster.getXMeasurements();
	testcluster.getZMeasurements();

	if (verbosity >= 1){
		testcluster.printSurf();///debug
		if (verbosity == 2) {
			
		}
		cout << "t(decodeWithMWPM):"<< double(clock()-start_t)/CLOCKS_PER_SEC << endl;//timing
		start_t = clock();
	}
	
	cout << "Parity" << parity <<endl;
}

void generator(const int trials, const double p, const double q, int seed, const int Np, const string fname, surfacetype surf, noisemodel N, bool verbosity, bool thread, bool use_env){
	ofstream outfile;
	outfile.open(fname);
	
	random_device rd;
	mt19937 engine{rd()};
	if (seed != 0) {
		engine.seed(seed);
	}
	
	int i = 0;
	while(true){
		if (verbosity > 0 && i%(trials/10) == 0) {
			cout << ".";
			cout.flush();
		}
		cluster testcluster({Sw+2,Sw+2,Sw+2,0},surf);
		int num_add = 0;
		if (N == GATE) {
			testcluster.addGateNoise(p, N);
		} else {
			testcluster.addNoise(p, q, N);
		}
		testcluster.getXMeasurements();
		testcluster.getZMeasurements();
		//generate mode 1 (random choose syndrom, random choose one neighbor)
		vector<int> test_pos;
		for (int c = 0; c < (Sw+2)*(Sw+2)*(Sw+2); c++) {
			if(testcluster.stabs[c] < 0){
				test_pos.push_back(c);
			}
		}
		int size = test_pos.size();
		if (size == 0){ //no syndrome present
			continue;
		} else {
			uniform_int_distribution<int> dist(0,size-1);
			int h = test_pos[dist(engine)]; //random choose syndrom
			int c; //choose random neighbor
			uniform_int_distribution<int> dist2(0,3);

			vertex aVertex(h,{Sw+2,Sw+2,Sw+2,0});
			c = aVertex.physical_qubits[dist2(engine)];
			
			coord C(c,{Sw+2,Sw+2,Sw+2,0});
			vector<float> window = testcluster.getWindow(C);
			
			int result;

			result = -(testcluster.c_error_pos[c]-1)/2;
			
			// result = -(testcluster.d_error_pos[c-(Sw+2)*(Sw+2)*(Sw+2)]-1)/2;
			outfile << window << result << endl;
			i++;
		}
		if (i == trials) {
			return;
		}
	}
}

int loopDecoding(string directory, string model_name, const int Lmin, const int Lmax, const int trials, const double pmin, const double pmax, const int Np, const double qmin, const double qmax, const int Nq, const string fname, surfacetype surf, noisemodel N, lossmodel L, int verbosity = 0, int thread = 1, bool make_corrections = 0, bool decode_with_NN = 0, bool dir = 0){
	ofstream outfile;
	outfile.open(directory + fname);
	outfile << "L,p_error,num_success\n";
	int p_bins = 0 == Np ? 1 : Np;
	int q_bins = 0 == Nq ? 1 : Nq;
	int num_correct;

	for (int l = Lmin; L <= Lmax; l=l+2) {
		cluster testcluster({L,L,L,0}, surf);
		cout << endl;
		for (int iq = 0; iq <= Nq; iq ++) {
			double q = qmin + (qmax-qmin)/q_bins * iq;
			
			for (int ip = 1; ip <= Np; ip ++) {
				double p = pmin + (pmax-pmin)/(p_bins + 1) * ip;
				float p_rounded = round(p*10000)/10000;
				string p_name = to_string(p);
				while (p_name.back()=='0'){
					p_name.pop_back();
				}
				//run simulation
				int num_correct = 0;

				if (decode_with_NN && surf != PLANE){
					cppflow::model model(directory + "models/" + model_name + p_name);//absolute directory to my model
					if (thread) {
						parallel_for(trials, [&](int start, int end){
							cluster testcluster({l,l,l,0}, surf);
							for(int i = start; i < end; ++i){
								if (N == GATE) {
									testcluster.addGateNoise(p, N);
								} else {
									testcluster.addNoise(p, q, N);
								}
								if (L != OFF_loss){
									testcluster.addLoss(q, L);
									try {
										testcluster.getSuperChunks();
									} catch (...) {
										continue;
									}
								}
								testcluster.getXMeasurements();
								testcluster.getZMeasurements();

								testcluster.decodeWithNN(model);
								testcluster.getXMeasurements();
								testcluster.getZMeasurements();
								//classical decoding (without NN)
								if (N != LOSS){
									if (testcluster.decodeWithMWPM(verbosity, dir, 0) == 1) {
										num_correct ++; //correction successful
									}
								} else {
									if (testcluster.decodeWithMWPMloss(verbosity, dir, 0) == 1) {
										num_correct ++; //correction successful
									}
								}		
							}
						});
					} else {
						for(int i = 0; i < trials; ++i){
							if (N == GATE) {
								testcluster.addGateNoise(p, N);
							} else {
								testcluster.addNoise(p, q, N);
							}
							if (L != OFF_loss){
								testcluster.addLoss(q, L);
								try {
									testcluster.getSuperChunks();
								} catch (...) {
									continue;
								}
							}
							testcluster.getXMeasurements();
							testcluster.getZMeasurements();

							testcluster.decodeWithNN(model);
							testcluster.getXMeasurements();
							testcluster.getZMeasurements();
							//classical decoding (without NN)
							if (N != LOSS){
								if (testcluster.decodeWithMWPM(verbosity, dir, 0) == 1) {
									num_correct ++; //correction successful
								}
							} else {
								if (testcluster.decodeWithMWPMloss(verbosity, dir, 0) == 1) {
									num_correct ++; //correction successful
								}
							}	
						}
					}
				} else {
					if (thread) {
						parallel_for(trials, [&](int start, int end){
							cluster testcluster({l,l,l,0}, surf);
							for(int i = start; i < end; ++i){
								if (N == GATE) {
									testcluster.addGateNoise(p, N);
								} else {
									testcluster.addNoise(p, q, N);
								}
								if (L != OFF_loss){
									testcluster.addLoss(q, L);
									try {
										testcluster.getSuperChunks();
									} catch (...) {
										continue;
									}
								}
								testcluster.getXMeasurements();
								testcluster.getZMeasurements();
								if (testcluster.decodeWithMWPM(verbosity, dir, make_corrections) == 1) {
									num_correct ++; //correction successful
								}
							}
						});
					} else {
						for(int i = 0; i < trials; ++i){
							if (N == GATE) {
								testcluster.addGateNoise(p, N);
							} else {
								testcluster.addNoise(p, q, N);
							}
							if (L != OFF_loss){
								testcluster.addLoss(q, L);
								try {
									testcluster.getSuperChunks();
								} catch (...) {
									continue;
								}
							}
							testcluster.getXMeasurements();
							testcluster.getZMeasurements();
							if (testcluster.decodeWithMWPM(verbosity, dir, make_corrections) == 1) {
								num_correct ++; //correction successful
							}
						}
					}
				}
				
				//printout/outfile
				if (verbosity >= 2) {
					cout << "t=" << double(clock()-start_t)/CLOCKS_PER_SEC << endl;
					start_t = clock();
				} else if (verbosity >= 1){
					cout << L << "," << p << "," << q << "," << num_correct << endl;
				} else {
					cout << ".";
					cout.flush();
				}
				outfile << L << "," << p << "," << q << "," << num_correct << "\n";
			}
		}
	}
	auto end = chrono::steady_clock::now();
	auto diff = end - start;
	cout << "total CPU time:" << double(clock()-start_t)/CLOCKS_PER_SEC <<";total time:"<< chrono::duration <double, milli> (diff).count() << "ms" << endl;
	outfile.close();
	return num_correct;
}

int main(int argc, const char *argv[]) {
	string fname,directory, model_name;
	surfacetype s;
	noisemodel N;
	lossmodel L;
	bool test, use_env, thread, make_corrections, generate, decode_with_NN, dir;
	int Lmin, Lmax, n, Np, Nq, seed, times, verbosity, return_value;
	float pmin, pmax, qmin, qmax;
	cxxopts::Options options(*argv,
							 "Simulator for fault-tolerant measurement-based quantum "
							 "computation on foliated surface code cluster states"
							 );
	options.add_options()
	("fname", "filename", cxxopts::value(fname)->default_value(""))
	("d, directory", "model directory", cxxopts::value(directory)->default_value("/users/VanLadmon/OneDrive - Stanford/PHYSICS/Research/Patrick/ML projects/ML3D/"))
	("m, model", "model name", cxxopts::value(model_name)->default_value("model,L=3(7),layer=5x512,epochs=1000,p="))
	("s", "surface type", cxxopts::value(s)->default_value("PLANE"))
	("Lmin", "Minimal size of mesh", cxxopts::value(Lmin)->default_value("3"))
	("Lmax", "Maximal size of mesh", cxxopts::value(Lmax)->default_value("17"))
	("n", "Number of trials", cxxopts::value(n)->default_value("10000"))
	("dir", "decoding direction", cxxopts::value(dir)->default_value("0"))
	
	("N", "noise model", cxxopts::value(N)->default_value("OFF"))
	("Np", "single c_error_pos error Points", cxxopts::value(Np)->default_value("10"))
	("pmin", "Minimal single c_error_pos error probability/test probability", cxxopts::value(pmin)->default_value("0"))
	("pmax", "Maximal single c_error_pos error probability")
	
	("L", "loss model", cxxopts::value(L)->default_value("OFF_loss"))
	("Nq", "correlated error points", cxxopts::value(Nq)->default_value("10"))
	("qmin", "Minimal correlated/loss error probability", cxxopts::value(qmin)->default_value("0"))
	("qmax", "Maximal correlated/loss error probability")
	
	("v, verbosity", "verbosity switch", cxxopts::value(verbosity)->default_value("0"))
	("generate", "generation mode switch", cxxopts::value(generate)->default_value("0"))
	("seed", "seed switch", cxxopts::value(seed)->default_value("0"))
	("times", "test times switch", cxxopts::value(times)->default_value("0"))
	("thread", "thread switch", cxxopts::value(thread))
	("use_env", "use environment variables", cxxopts::value(use_env))
	("test", "test switch", cxxopts::value(test))
	("make_corrections", "correct qubit", cxxopts::value(make_corrections)->default_value("0"))
	("decode_with_NN", "turn on NN decoder", cxxopts::value(decode_with_NN)->default_value("0"));
	
	options.parse(argc, argv);

	if(verbosity >= 1){
		cout <<"use_env:" << use_env << "thread" << thread << "test" << test << endl;
		cout << "Lmin:" << Lmin << ";Lmax:" << Lmax << endl;
		cout << "Pmin:" << pmin << ";Pmax:" << pmax << ";nP" << Np << endl;
		cout << "Qmin:" << qmin << ";Qmax:" << qmax << ";nQ" << Nq << endl;
		cout << "n:" << n << ";seed:" << seed << ";thread:" << thread << ";NN:" << decode_with_NN << endl;
		cout << "v:" << verbosity << ";dir:" << dir << endl;
		cout << "noise model (N):" << to_string(N) << endl;
		cout << "loss model (L):" << to_string(L) << endl;
		//outputing options;
		cout << "hardware_concurrency:" << thread::hardware_concurrency() << endl;
		if (use_env) {
			cout << "SLURM_CPUS_PER_TASK:" << atoi(getenv("SLURM_CPUS_PER_TASK")) << endl;
		}
	}
	
	
	
	if (fname == "") {
		fname = "L=" + to_string(Lmin) + ",P=(" + to_string(pmin).substr(3,2) + "," + to_string(pmax).substr(3,2) + "),n=" +to_string(n) + to_string(s) + "," + to_string(N) + ".out";
	}
	if (test) {
		string p_name = to_string(pmin);
		while (p_name.back()=='0'){
			p_name.pop_back();
		}
		cppflow::model model(directory + "models/" + model_name + p_name);//absolute directory to my_model
		cluster testcluster({Lmin,Lmin,Lmin,0}, s);
		int i = 0;
		do {
			testDecoding(model, testcluster, pmin, qmin, seed + i, s, N, L, verbosity, decode_with_NN);
			i++;
		}
		while (i < times);
	} else if (generate) {
		if (fname == "") {
			fname = "train_data/train_set_L=" + to_string(Sw) + ",P=(" + to_string(pmin).substr(3,2) + "," + to_string(pmax).substr(3,2) + "),n=" +to_string(n) + ".out";
		}
		generator(n, pmin, qmin, seed, Np-1, directory + fname, s, N, verbosity, thread, use_env);
	} else{
		loopDecoding(directory, model_name, Lmin, Lmax, n, pmin, pmax, Np-1, qmin, qmax, Nq-1, fname, s, N, L, verbosity, thread, make_corrections, decode_with_NN, dir);
	}
	return return_value;
}
/////////// INDEP noise ////////////

// ./simulate -s TORUS --pmin 0 --pmax 0.04  --Np 10 --Nq 1 --Lmin 3 -v 1 -d "/scratch/users/ladmon/" -m 'model,L=3(7),layer=5x512,epochs=1000,p='
// ./simulate -s TORUS --pmin 0 --pmax 0.04 --Np 100 -n 10000 --Lmin 3 --Lmax 17 -v 1 -d "/scratch/users/ladmon/ML3D/" --decode_with_NN -m 'model,L=3(7),layer=5x512,epochs=1000,p='
// ./simulate -s PLANE --pmin 0.03 --pmax 0.04 --Np 100 -n 10000 --Lmin 10 --Lmax 17 -v 2 -d "/scratch/users/ladmon/ML3D/" --test

/////////// DEPOL noise ///////////
//test
// ./simulate -s TORUS --pmin 0.006 -N DEPOL1 --Np 20 --Nq 1 -n 1000 --Lmin 10 --Lmax 17 -v 1 --fname '/ML3D/results/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out' -d /scratch/users/ladmon/ML3D/ -m 'depol1model,L=5(7),layer=5x256,epochs=10000,p=' --decode_with_NN --test --seed 1

/////////// GATE noise/ //////////
// ./simulate -s TORUS --pmin 0 -N GATE --pmax 0.007 --Np 10 --Nq 1 -n 10000 --Lmin 3 --Lmax 17 -v 1 -d "/scratch/users/ladmon/ML3D/" -m 'model,L=3(7),layer=5x512,epochs=1000,p=' 
// ./simulate -s TORUS --pmin 0 -N GATE --pmax 0.04 --Np 100 -n 10000 --Lmin 3 --Lmax 17 -v 1 -d "/scratch/users/ladmon/ML3D/" -m 'model,L=3(7),layer=5x512,epochs=1000,p='
//test
// ./simulate -s TORUS -N GATE --pmin 0.00035 --Np 100 -n 10000 --Lmin 7 --Lmax 17 -v 2 -d "/scratch/users/ladmon/ML3D/" --decode_with_NN -m 'gatemodelS,L=3(5),layer=3x128,epochs=10000,p=' --test

//large test
///./simulate -s PLANE --pmin 0 --pmax 0.04  --qmin 0 --qmax 0.04 --Np 10  --Nq 10 -n 500 --Lmin 3 -v 1
// ./simulate -s TORUS --pmin 0 --pmax 0.04 --Np 100 -n 10000 --Lmin 3 --Lmax 17 -v 2 -d "/scratch/users/ladmon/ML3D/" --decode_with_NN -m 'model,L=3(7),layer=5x512,epochs=1000,p=' --test -e


// ./simulate -s TORUS --pmin 0 -N TEST --pmax 0.04 --Np 10 --Nq 1 -n 10000 --Lmin 3 --Lmax 17 -v 1 --thread --fname '/ML3D/results/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out' -d /scratch/users/ladmon/ -m 'model,L=3(7),layer=5x512,epochs=1000,p='
// ./simulate -s TORUS --pmin 0 -N TEST --pmax 0.04 --Np 10 --Nq 1 -n 10000 --Lmin 3 --Lmax 17 -v 1 --thread --fname '/ML3D/results/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out' -d /scratch/users/ladmon/ -m 'model,L=3(7),layer=5x512,epochs=1000,p='

//PLANE test

//LOSS