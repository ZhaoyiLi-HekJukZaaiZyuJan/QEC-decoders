# include "vertex.hpp"
# include "cluster.hpp"
# include "functions.hpp"

#include <random>

//for chunk graphics
std::vector<string> colorvector= { "\033[1;31m◯\033[0m","\033[1;32m◯\033[0m","\033[1;33m◯\033[0m","\033[1;34m◯\033[0m","\033[1;35m◯\033[0m","\033[1;36m◯\033[0m","\033[1;37m◯\033[0m"};

cluster::cluster(const coord& S, const surfacetype& this_surf){
	this->S = S;
	this->this_surf = this_surf;
	vector<int> surf(2*S.x*S.y*S.z,1);//decode
	this->surf = surf;//decode
	
	vector<int> c_error_pos(3*S.x*S.y*S.z,1);
	this->c_error_pos = c_error_pos;
	vector<int> c_loss_pos(3*S.x*S.y*S.z,1);
	this->c_loss_pos = c_loss_pos;
	vector<int> stabs(2*S.x*S.y*S.z,1);
	this->stabs = stabs;
}

void cluster::addPauli(const double & p1, const double & p2, const int& seed){
	//Total heuristic Probabilities
	//Initialization of error operator
	fill(c_error_pos.begin(), c_error_pos.end(), 1);
	random_device rd;
	mt19937 engine{rd()};
	if (seed != 0) {
		engine.seed(seed);
	}
	uniform_real_distribution<> dist(0.0, 1.0);
	
	
	//Error Model 1: uncorrelated error distribution for all physical c_error_poss
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		if(dist(engine) < p1) {
		c_error_pos[c] = -1;
		} else {
		c_error_pos[c] = 1;
		}
	}
	//Error Model 2: Correlation Error on the Same Edge
	for (int c = 0; c < S.x*S.y*S.z; c++) {
		coord C(c,S);
		if (C.z == S.z - 1 || C.x == 0) {//remove boundary cubes
			continue;
		}
		for (int face = 0; face < 3; face ++) {
			for (int dir = 0; dir < 2; dir ++) {
				if (dist(engine) < p2) {
					c_error_pos[C.getFaceQubits(S, face, dir, 3)] *= -1;
					c_error_pos[C.getFaceQubits(S, face, dir, 2)] *= -1;
				}
			}
		}
	}
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		coord C(c, S);
		if(this_surf == PLANE && ((C.l == 1 && (C.y == S.y - 1 || C.x == 0)) ||(C.l == 2 && (C.z == S.z - 1 || C.x == 0)))){
			c_error_pos[c] = 1; //	for planar code, remove errors on left, more, and lower boundaries to create edges.
		} if(this_surf == TORUS && C.l == 2 && C.z == S.z - 1){
			c_error_pos[c] = 1;	//	for planar code, remove time like errors on left, more, and lower boundaries to create edges.
		}
	}
}
//differentiate correlated and uncorrelated errors

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

void cluster::addLoss(const double & p1, const double & p2, const int& seed=0){
	fill(c_loss_pos.begin(), c_loss_pos.end(), 1);
	//Total heuristic Probabilities
	random_device rd;
	mt19937 engine{rd()};
	if (seed != 0) {
		engine.seed(seed+1);
	}
	uniform_real_distribution<> dist(0.0, 1.0);
	
	//Initialization of error operator
	//Error Model 1: uncorrelated error distribution for all physical c_error_poss
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		if(dist(engine) < p1) {
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

void cluster::addGateNoise(const double & p, const int& seed=0){
	//Total heuristic Probabilities
	random_device rd;
	mt19937 engine{rd()};
	if (seed != 0) {
		engine.seed(seed);
	}
	uniform_real_distribution<> dist(0.0, 1.0);
	
	//////// Initialization of error operator //////// 
	fill(c_error_pos.begin(), c_error_pos.end(), 1);
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		for (int i = 0; i < 3; i++) { //p_S, p_M, p_P
			if(dist(engine) < p*2/3) {
				c_error_pos[c] *= -1;
			}
		}
		for (int i = 0; i < 2; i++) { //single qubit errors
			if(dist(engine) < p*8/15) {
				c_error_pos[c] *= -1;
			}
		}
	}

	//////// CZ gates //////// 
	for (int c = 0; c < S.x*S.y*S.z; c++) {
		coord C(c,S);
		for (int face = 0; face < 3; face++) {
			for (int dir = 0; dir < 2; dir ++) {//verticle/horizontal
				double p1 = dist(engine);
				if (p1 < p*4/15) {
					c_error_pos[C.getFaceQubits(S, face, dir, 3)] *= -1;//3
				} else if (p*4/15 < p1 && p1 < p*8/15) {
					c_error_pos[C.getFaceQubits(S, face, dir, 2)] *= -1;//2
				} else if (p*8/15 < p1 && p1 < p*12/15){
					c_error_pos[C.getFaceQubits(S, face, dir, 3)] *= -1;//3
					c_error_pos[C.getFaceQubits(S, face, dir, 2)] *= -1;//2
				}
			}
		}
	}
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		coord C(c, S);
		if(this_surf == PLANE && ((C.l == 1 && (C.y == S.y - 1 || C.x == 0)) ||(C.l == 2 && (C.z == S.z - 1 || C.x == 0)))){
			c_error_pos[c] = 1;
		}//	for planar code, remove errors on left, more, and lower boundaries to create edges.
	}
}

void cluster::addFullGateNoise(const double & p, const int& seed=0){
	//Total heuristic Probabilities
	random_device rd;
	mt19937 engine{rd()};
	if (seed != 0) {
		engine.seed(seed);
	}
	uniform_real_distribution<> dist(0.0, 1.0);
	
	//////// Initialization of error operator ////////
	fill(c_error_pos.begin(), c_error_pos.end(), 1);
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		for (int i = 0; i < 3; i++) { //p_S, p_M, p_P
			if(dist(engine) < p*2/3) {
				c_error_pos[c] *= -1;
			}
		}
	}

	//////// CZ gates ////////
	for (int c = 0; c < S.x*S.y*S.z; c++) {
		coord C(c,S);
		if (C.z == S.z - 1 || C.y == S.y - 1) {//remove boundary cubes
			continue;
		}
		for (int face = 0; face < 3; face++) {
			double p1 = dist(engine); //process 1 (black)
			double p2 = dist(engine); //process 2 (pink)
			double p3 = dist(engine); //process 3 (rose)
			double p4 = dist(engine); //process 4 (purple)

			if (p1 < p*4/15) {
				c_error_pos[C.getFaceQubits(S, face, 0, 3)] *= -1;//3
			} else if (p*4/15 < p1 && p1 < p*8/15){
				c_error_pos[C.getFaceQubits(S, face, 0, 0)] *= -1;//0
				c_error_pos[C.getFaceQubits(S, face, 0, 1)] *= -1;//1
				c_error_pos[C.getFaceQubits(S, face, 0, 2)] *= -1;//2
			} else if (p*8/15 < p1 && p1 < p*12/15){
				c_error_pos[C.getFaceQubits(S, face, 0, 0)] *= -1;//0
				c_error_pos[C.getFaceQubits(S, face, 0, 1)] *= -1;//1
				c_error_pos[C.getFaceQubits(S, face, 0, 2)] *= -1;//2
				c_error_pos[C.getFaceQubits(S, face, 0, 3)] *= -1;//3
			}
			
			if (p2 < p*4/15) {
				c_error_pos[C.getFaceQubits(S, face, 0, 2)] *= -1;//2
			} else if (p*4/15 < p2 && p2 < p*8/15) {
				c_error_pos[C.getFaceQubits(S, face, 0, 0)] *= -1;//0
				c_error_pos[C.getFaceQubits(S, face, 0, 1)] *= -1;//1
			} else if (p*8/15 < p2 && p2 < p*12/15) {			
				c_error_pos[C.getFaceQubits(S, face, 0, 0)] *= -1;//0
				c_error_pos[C.getFaceQubits(S, face, 0, 1)] *= -1;//1
				c_error_pos[C.getFaceQubits(S, face, 0, 2)] *= -1;//2	
			}
			
			if (p3 < p*4/15) {
				c_error_pos[C.getFaceQubits(S, face, 0, 1)] *= -1;//1
			} else if (p*4/15 < p3 && p3 < p*8/15){
				c_error_pos[C.getFaceQubits(S, face, 0, 0)] *= -1;//0
			} else if (p*8/15 < p3 && p3 < p*12/15){
				c_error_pos[C.getFaceQubits(S, face, 0, 0)] *= -1;//0
				c_error_pos[C.getFaceQubits(S, face, 0, 1)] *= -1;//1
			}

			if (p4 < p*8/15) { //purple
				c_error_pos[C.getFaceQubits(S, face, 0, 0)] *= -1;//0
			}
		}
	}
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		coord C(c, S);
		if(this_surf == PLANE && ((C.l == 1 && (C.y == S.y - 1 || C.x == 0)) ||(C.l == 2 && (C.z == S.z - 1 || C.x == 0)))){
			c_error_pos[c] = 1;
		}//	for planar code, remove errors on left, more, and lower boundaries to create edges.
	}
}

void cluster::addBiasedGateNoise(const double & p, const double & B, const int& seed=0){
	//Total heuristic Probabilities
	random_device rd;
	mt19937 engine{rd()};
	if (seed != 0) {
		engine.seed(seed);
	}
	uniform_real_distribution<> dist(0.0, 1.0);
	
	//Initialization of error operator
	fill(c_error_pos.begin(), c_error_pos.end(), 1);
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		for (int i = 0; i < 2; i++) {
			if(dist(engine) < p*2/3) {
				c_error_pos[c] *= -1;
			}
		}
	}

	for (int c = 0; c < S.x*S.y*S.z; c++) {
		coord C(c,S);
		if (C.z == S.z - 1 || C.x == 0) {//remove boundary cubes
			continue;
		}
		for (int face = 0; face < 3; face++) {
			for (int dir = 0; dir < 2; dir ++) {//verticle/horizontal
				double pz = dist(engine);
				double px = dist(engine);
				if (pz < p*2/3) {
					c_error_pos[C.getFaceQubits(S, face, dir, 3)] *= -1;//3
				}
				if (px < p*2/3/B) {
					c_error_pos[C.getFaceQubits(S, face, dir, 1)] *= -1;//1 //this noise model is symmetric under switch 2<->1
				}
			}
			double p_black_x = dist(engine);
			double p_black_z = dist(engine);
			double p_purple = dist(engine);
			if (p_purple < p*2/3) {
				c_error_pos[C.getFaceQubits(S, face, 0, 0)] *= -1;//0
			}
			if (p_black_z < p*2/3) {
				c_error_pos[C.getFaceQubits(S, face, 0, 3)] *= -1;//3
			}
			if (p_black_x < p*2/3/B) {
				c_error_pos[C.getFaceQubits(S, face, 0, 3)] *= -1;//3
			}
		}
	}
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		coord C(c, S);
		if(this_surf == PLANE && ((C.l == 1 && (C.y == S.y - 1 || C.x == 0)) ||(C.l == 2 && (C.z == S.z - 1 || C.x == 0)))){
			c_error_pos[c] = 1;
		}//	for planar code, remove errors on left, more, and lower boundaries to create edges.
	}
}

void cluster::addError(){ //for debugging
	//Total heuristic Probabilities	
	//add errors manually here
	coord C1(3,3,3,0);
	coord C2(4,3,3,0);
	coord C3(1,1,1,0);
	coord C4(1,2,1,0);
	for (int k = 0; k < 3; k ++){
		c_error_pos[C1.getFaceQubits(S,2,0,k)] *= -1;
		c_error_pos[C2.getFaceQubits(S,2,0,k)] *= -1;
		c_error_pos[C3.getFaceQubits(S,1,0,k)] *= -1;
		c_error_pos[C4.getFaceQubits(S,1,0,k)] *= -1;
	}

	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		coord C(c, S);
		if(this_surf == PLANE && ((C.l == 1 && (C.y == S.y - 1 || C.x == 0)) ||(C.l == 2 && (C.z == S.z - 1 || C.x == 0)))){
			c_error_pos[c] = 1;
		}//	for planar code, remove errors on left, more, and lower boundaries to create edges.
	}
}

//No seed used for seed = 0
void cluster::addNoise(const double & p, const double & q, const noisemodel N, const lossmodel L, const int& seed){
	// Pauli part
	if (N == GATE) {
		addGateNoise(p, seed);
	}
	if (N == GATE_full) {
		addFullGateNoise(p, seed);
	}
	if (N == GATE_biased) { //p:error probability, //q: bias
		addBiasedGateNoise(p, q, seed);
	}
	if (N == INDEP_211) {
		addPauli(p*p,0,seed);
		addLoss(2*p*(1-p),0,seed);
	} 
	if (N != GATE && N != GATE_full && N != GATE_biased && N != INDEP_211) {
		addPauli(NOISEMODELMAP[N](p,0).first, NOISEMODELMAP[N](p,0).second,  seed);
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
			if(C.l == 0 && C.x == 0){//on left face
				int cV = coord(divmod(C.x+1, S.x), C.y, C.z, 0).hash(S);
				int i = findParent(cV);
				boundary_info[i] = 1;//label suppress generation of boundary nodes
				left_boundary_chunks.push_back(i);
			} else if(C.l == 0 && C.x == S.x-1){// on right face
				int cV = coord(C.x, C.y, C.z, 0).hash(S);
				int i = findParent(cV);
				boundary_info[i] = 1;//label suppress generation of boundary nodes
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

void cluster::getStabs(){
	stabs = {};
	for (int c = 0; c < S.x*S.y*S.z; c++){
		stabs.push_back(1);//initialization
		//measurement of vertex operator
		vertex aVertex(c, S);
		for (int pos = 0; pos < 6; pos ++) {
			if (c_error_pos[aVertex.partial[pos]] == -1) {
				stabs[c] *= -1;
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
	printMatrix(print_out);
}

void cluster::printPrimal(surfacetype s = TORUS){
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
				}
			}
		}
		cout << "layer:" << k << endl;
		printMatrix(print_out);
	}
	
}


int cluster::getLeftDistance(const coord& S, const vector<int>& chunk){
	int min_dist = INT_MAX;
	for (int i = 0; i < chunk.size(); i++) {
		int dist = coord(chunk[i],S).x;
		min_dist = (dist < min_dist) ? dist : min_dist;
	}
	for (int i = 0; i < left_boundary_chunks.size(); i++) {
		int dist = getTaxicabDistanceChunks(S, super_chunks[left_boundary_chunks[i]], chunk, PLANE);
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
		int dist = getTaxicabDistanceChunks(S, super_chunks[right_boundary_chunks[i]], chunk, PLANE);
		min_dist = (dist < min_dist) ? dist : min_dist;
	}
	return min_dist;
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

void cluster::surfaceCorrect(PerfectMatching* pm, const int& vertices_num, const vector<int>& vertexPosition, const vector<int>& boundary_nodes, const int& verbose, const surfacetype& s){
	//Find the error operator of each pair.
	vector<int> error_op_Z_pos; //correction operator not reduced
	vector<int> matchPosition; //matching vertex of vertexPosition respectively
	for (int i = 0; i < vertices_num; i++) {
		int ac = vertexPosition[i];
		coord relative_pos; // relative position between two vertices in a pair
		//get relative position (vector) of pair to determine the path in between.
		vertex aVertex(ac, S), bVertex;
		if (pm->GetMatch(i) < vertices_num){
			int bc = vertexPosition[pm->GetMatch(i)];// position of vertexa's match, vertexb
			matchPosition.push_back(bc);
			if (count(matchPosition.begin(), matchPosition.end(), ac)) continue; //Prevent recounting of vertices
			bVertex = vertex(bc, S);
			relative_pos = getTaxicabDisplacement(S, ac, bc, s);
		}else{ // matched to boundary
			matchPosition.push_back(0);
			relative_pos = {boundary_nodes[i], 0, 0, 0};
		}
		int n;
			if (relative_pos.x > 0) {//a to the left of b
				n = aVertex.partial[1];//use right c_error_pos
				for(int i =0; i < abs(relative_pos.x); i++){
					error_op_Z_pos.push_back(coord(divmod(coord(n,S).x + i, S.x), coord(n,S).y, 0, 0).hash(S));
				}
			} else if(relative_pos.x < 0) {//a to the right of b
				n = aVertex.partial[0];//use left c_error_pos
				for(int i =0; i < abs(relative_pos.x); i++){
					error_op_Z_pos.push_back(coord(divmod(coord(n,S).x - i, S.x), coord(n,S).y, 0, 0).hash(S));
				}
			}
			if (relative_pos.y > 0) {//a above b
				n = bVertex.partial[2];//use upper c_error_pos
				for(int i =0; i < abs(relative_pos.y); i++){
					error_op_Z_pos.push_back(coord(coord(n,S).x, divmod(coord(n,S).y - i, S.y), 0, 1).hash(S));
				}
			} else if (relative_pos.y < 0) {//a below b
				n = bVertex.partial[3];//use lower
				for(int i =0; i < abs(relative_pos.y); i++){
					error_op_Z_pos.push_back(coord(coord(n,S).x, divmod(coord(n,S).y + i, S.y), 0, 1).hash(S));
				}
			}
			cout << "this vertex:" << coord(ac,S) << endl; 
			cout << "match vertex:" << coord(vertexPosition[pm->GetMatch(i)],S) << endl; 
			cout << "this vertex op:";
			for (int i = 0 ; i <error_op_Z_pos.size(); i++) {
				cout << coord(error_op_Z_pos[i],S) << " ";
			}
			cout << endl;
	}

	for (int i = 0; i < error_op_Z_pos.size(); i++) {
		surf[error_op_Z_pos[i]] *= -1; //act local Z operator
	}
	if(verbose == 2){
		cout << "vertex/match/bn: "<< endl;
		//printout test
		for (int i = 0 ; i < vertices_num; i++) {
			cout << coord(vertexPosition[i],S)<< "/";
			cout << coord(matchPosition[i],S)<<"/";
			if (s == PLANE){
				cout << boundary_nodes[i] << "/";
			}	
			cout << getTaxicabDistance(S,coord(vertexPosition[i],S),coord(matchPosition[i],S),s) << endl;
		}
		cout << "errorOps: " << endl;
		for (int i = 0 ; i <error_op_Z_pos.size(); i++) {
			cout << coord(error_op_Z_pos[i],S) << " ";
		}
		cout<<endl;
	}
}


void cluster::surfaceCorrectLoss(PerfectMatching* pm, const int& vertices_num, const vector<int>& wrong_chunks, const vector<int>& boundary_info, const int& verbose, const surfacetype& s){
		//Find the error operator of each pair.
		vector<int> error_op_Z_pos; //correction operator not reduced
		vector<int> matchPosition; //matching vertex of wrong_chunks respectively
		for (int i = 0; i < vertices_num; i++) {
			int ac = super_chunks[wrong_chunks[i]][0];
			coord relative_pos; // relative position between two vertices in a pair
			//get relative position (vector) of pair to determine the path in between.
			vertex aVertex(ac, S), bVertex;
			if (pm->GetMatch(i) < vertices_num){
				int bc = super_chunks[wrong_chunks[pm->GetMatch(i)]][0];// position of vertexa's match, vertexb
				matchPosition.push_back(bc);
				if (count(matchPosition.begin(), matchPosition.end(), ac)) continue; //Prevent recounting of vertices
				bVertex = vertex(bc, S);
				relative_pos = getTaxicabDisplacement(S, ac, bc, s);
			}else{ // matched to boundary
				matchPosition.push_back(0);
				int x =  coord(aVertex.c,S).x;
				if (x <= S.x-x) {
					relative_pos = {-x, 0, 0, 0};
				} else {
					relative_pos = {S.x-x, 0, 0, 0};
				}
			}
			int n;
			if (relative_pos.x > 0) {//a to the left of b
				n = aVertex.partial[1];//use right c_error_pos
				for(int i =0; i < abs(relative_pos.x); i++){
					error_op_Z_pos.push_back(coord(divmod(coord(n,S).x + i, S.x), coord(n,S).y, 0, 0).hash(S));
				}
			} else if(relative_pos.x < 0) {//a to the right of b
				n = aVertex.partial[0];//use left c_error_pos
				for(int i =0; i < abs(relative_pos.x); i++){
					error_op_Z_pos.push_back(coord(divmod(coord(n,S).x - i, S.x), coord(n,S).y, 0, 0).hash(S));
				}
			}
			if (relative_pos.y > 0) {//a above b
				n = bVertex.partial[2];//use upper c_error_pos
				for(int i =0; i < abs(relative_pos.y); i++){
					error_op_Z_pos.push_back(coord(coord(n,S).x, divmod(coord(n,S).y - i, S.y), 0, 1).hash(S));
				}
			} else if (relative_pos.y < 0) {//a below b
				n = bVertex.partial[3];//use lower
				for(int i =0; i < abs(relative_pos.y); i++){
					error_op_Z_pos.push_back(coord(coord(n,S).x, divmod(coord(n,S).y + i, S.y), 0, 1).hash(S));
				}
			}			
	
		}
	
		for (int i = 0; i < error_op_Z_pos.size(); i++) {
			surf[error_op_Z_pos[i]] *= -1; //act local Z operator
		}
		if(verbose == 2){
			cout << "vertex/match/bn_info/: "<< endl;
			//printout test
			for (int i = 0 ; i < vertices_num; i++) {
				cout << coord(super_chunks[wrong_chunks[i]][0],S)<< "/";
				cout << coord(matchPosition[i],S)<<"/";
				cout << boundary_info[i] <<"/";
			}
			cout << "errorOps: " << endl;
			for (int i = 0 ; i <error_op_Z_pos.size(); i++) {
				cout << coord(error_op_Z_pos[i],S) << " ";
			}
			cout << endl;
		}
}

int cluster::decodeWithMWPM(int verbosity, bool make_corrections, surfacetype surf){
	//PAIR MATCHING
	getStabs();
    vector<int> vertexPosition;///actual (non boundary) vertices to be matched
	
    for (int c = 0; c < S.x*S.y*S.z; c++) {
        if (stabs[c] == -1) {
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

	if (this_surf == PLANE) {// add in the boundary nodes
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
		surfaceCorrect(pm, vertices_num, vertexPosition, boundary_nodes, verbosity, surf);///debug
	}
	
	//alternative check
	int parity = 1;

	for (int i = 0; i < vertices_num; i++) {
		if (this_surf == PLANE){
			if (pm->GetMatch(i) == i + vertices_num && boundary_nodes[i] > 0) { //matched to right boundary
				parity *= -1;
			}
		} else if (this_surf == TORUS){		
			int aC = vertexPosition[i];
			int bC = vertexPosition[pm->GetMatch(i)];// position of vertexa's match, vertexb

			int x1 = aC % S.x;
			int x2 = bC % S.x;
			if ((2*abs(x1 - x2) > S.x && x1 < x2)){ //crossed boundary
				parity *= -1;
			}
		}
	}

	for (int i = 0; i < S.y; i++){//errors on right boundary
		for (int j = 0; j < S.z; j++){
			parity *= c_error_pos[coord(S.x-1, i, j, 0).hash(S)];
		}
	}
	
	delete pm;
	return parity;
}


int cluster::decodeWithMWPMLoss(int verbosity, bool make_corrections, surfacetype surf){
	//PAIR MATCHING
	getStabs();
    vector<int> wrong_chunks;//actual (non boundary) chunks to be matched
	//e.g. [1,3,5,7]
	//get wrong_chunks
	for (int i = 0; i < super_chunks.size(); i++) {
		if (boundary_info[i] != 0) {//don't consider boundary chunks
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
			int dist = getTaxicabDistanceChunks(S, super_chunks[wrong_chunks[i]], super_chunks[wrong_chunks[j]], surf);
	//			cout << wrong_chunks[i] << ":" << wrong_chunks[j]  << dist << endl; /// debug
			if (dist < dmin1 + dmin2) {//optimization
				pm->AddEdge(i,j,dist);
				pm->AddEdge(i + vertices_num, j + vertices_num, 0);//boundary interconnection
				
			}
		}
	}
	
	//solve the graph using MWPM decoder
	pm->Solve();

		//surface correction
	if(make_corrections){
		surfaceCorrectLoss(pm, vertices_num, wrong_chunks, boundary_info, verbosity, PLANE);
	}
	
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