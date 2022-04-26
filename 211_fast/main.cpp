//########################## [[2,1,1]] Toric Code ##########################
//########################## Faster Main function ##########################

# include <fstream>
# include <functional>
# include <algorithm>
# include <utility>
# include <iomanip>
# include "libs/cxxopts/cxxopts.hpp"
# include <assert.h>
# include "functions.cpp"
# include "thread_for.hpp"
# include <time.h>
#include <chrono>
# include "libs/blossom5-v2.05.src/PerfectMatching.h"
# include "libs/blossom5-v2.05.src/GEOM/GeomPerfectMatching.h"

using namespace std;

clock_t start_t = clock();
auto start = chrono::steady_clock::now();

class cluster {
	//422 structures
	vector<vector<int>> z_error_pos_211;
	//high level structures
	vector<int> z_error_pos;
	vector<int> x_vec;
	//loss, superchecks, and superchunks
	vector<int> loss_pos;
	
	vector<vector<int>> super_chunks;
	vector<int> left_boundary_chunks;
	vector<int> right_boundary_chunks;
	vector<int> boundary_info;
	
	coord S;
	surfacetype this_surf = PLANE;
	
	
	
public:
	cluster(const coord&, const surfacetype&);
	
	void addNoise(const double&, const noisemodel, const int&);
	void addGateNoise(const double&, const int&);
	
	void decode211(int);
	void print();
	
	void addLoss(const double&, const lossmodel, const int&);
	void getSuperChunks();
	void printSuperChunks();
	int findParent(int);
	int getLeftDistance(const coord&, const vector<int>&);
	int getRightDistance(const coord&, const vector<int>&);
	
	void getx_vec();
	void getXVec();
	int decodeWithMWPM(int, int);
	
//	vector<int> surf;///decode
//	void getSurf();///decode
//	void printSurf();///decode
//	void surfaceCorrect(PerfectMatching*, const int&, const vector<int>&, const vector<int>&, const int&);///decode
	
	int checkMeasurementOutcomeX1();
	int checkMeasurementOutcomeX2();
};

cluster::cluster(const coord& S, const surfacetype& this_surf){
	this->S = S;
	this->this_surf = this_surf;
	
//	vector<int> surf(2*S.x*S.y*S.z,1);
//	this->surf = surf;
	
	vector<int> z_error_pos(3*S.x*S.y*S.z,1);
	this->z_error_pos = z_error_pos;
	
	vector<int> loss_pos(3*S.x*S.y*S.z,1);
	this->loss_pos = loss_pos;
	
	vector<vector<int>> z_error_pos_211(3*S.x*S.y*S.z, vector<int>(2, 1));
	this->z_error_pos_211 = z_error_pos_211;
}
void cluster::addNoise(const double & p, const noisemodel N, const int& seed=0){
	//Total heuristic Probabilities
	random_device rd;
	mt19937 engine{rd()};
	if (seed != 0) {
		engine.seed(seed);
	}
	uniform_real_distribution<> dist(0.0, 1.0);
	//Initialization of error operator
	//Error Model 1: uncorrelated error distribution for all physical z_error_poss
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		for (int k = 0; k < 2; k++) {
			if(dist(engine) < NOISEMODELMAP[N](p).first) {
				z_error_pos_211[c][k] = -1;
			} else {
				z_error_pos_211[c][k] = 1;
			}
		}
	}
	//	//Error Model 2: Correlation Error on the Same Edge (check correctness)
	//	for (int c = 0; c < S.x*S.y*S.z; c++) {
	//		vertex aVertex(c, S); //get vertices as when getting x_vec's
	//		for (int pos = 0; pos < 6; pos = pos + 2) {
	//			if(dist(engine) < NOISEMODELMAP[N](p).second) {
	//				z_error_pos[aVertex.partial[pos]] *= -1;
	//				z_error_pos[aVertex.partial[pos + 1]] *= -1;
	//			}
	//		}
	//	}
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		coord C(c, S);
		if(this_surf == PLANE && ((C.l == 1 && (C.y == S.y - 1 || C.x == 0)) ||(C.l == 2 && (C.z == S.z - 1 || C.x == 0)))){
			z_error_pos[c] = 1;
		}//	for planar code, remove errors on left, more, and lower boundaries to create edges.
	}
}

void cluster::addGateNoise(const double & p, const int& seed=0){
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		z_error_pos_211[c][0] = 1;
		z_error_pos_211[c][1] = 1;
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
//			z_error_pos_211[c][1] *= -1;
//		}
//		if(p2 < p*2/3) {
//			z_error_pos_211[c][0] *= -1;
//			z_error_pos_211[c][1] *= -1;
//		}
//		if(pCX < p*4/15) {
//			z_error_pos_211[c][1] *= -1;
//		} else if(p*4/15 < pCX && pCX < p*8/15) {
//			z_error_pos_211[c][0] *= -1;
//		} else if(p*8/15 < pCX && pCX < p*12/15) {
//			z_error_pos_211[c][0] *= -1;
//			z_error_pos_211[c][1] *= -1;
//		}
//	}
	//initialization, measurement and storage errors
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		for (int i = 0; i < 2; i++) {// 211
//			if(dist(engine) < p*2-pow(p,2)*4/3) {
			if(dist(engine) < p*4/3) {
				z_error_pos_211[c][i] = -1;
			} else {
				z_error_pos_211[c][i] = 1;
			}
		}
	}
	
	for (int c = 0; c < S.x*S.y*S.z; c++) {
		coord C(c,S);
		if (C.z == S.z - 1 || C.x == 0) {//remove boundary cubes
			continue;
		}
		for (int i = 0; i < 2; i++) {//211
			for (int face = 0; face < 3; face ++) {
				for (int direction = 0; direction < 2; direction++) {//horizontal/vertical
				
					double p1 = dist(engine); //process 1 (black)
					double p2 = dist(engine); //process 2 (pink)
					double p3 = dist(engine); //process 3 (rose)
					double p4 = dist(engine); //process 4 (purple)
					

					if (p1 < p*8/15) { //black
						z_error_pos_211[C.getFaceQubits(S, face, direction, 3)][i] *= -1;
					}
					if (p4 < p*4/15) { //purple
						z_error_pos_211[C.getFaceQubits(S, face, direction, 0)][i] *= -1;
					} else if (p*4/15 < p4 && p4 < p*8/15) { //purple
						z_error_pos_211[C.getFaceQubits(S, face, direction, 1)][i] *= -1;
						z_error_pos_211[C.getFaceQubits(S, face, direction, 2)][i] *= -1;
						z_error_pos_211[C.getFaceQubits(S, face, direction, 3)][i] *= -1;
					} else if (p*8/15 < p4 && p4 < p*12/15) { //purple
						z_error_pos_211[C.getFaceQubits(S, face, direction, 0)][i] *= -1;
						z_error_pos_211[C.getFaceQubits(S, face, direction, 1)][i] *= -1;
						z_error_pos_211[C.getFaceQubits(S, face, direction, 2)][i] *= -1;
						z_error_pos_211[C.getFaceQubits(S, face, direction, 3)][i] *= -1;
					}
					
					
					if (p2 < p*4/15) { //pink
						z_error_pos_211[C.getFaceQubits(S, face, direction, 2)][i] *= -1;
					} else if (p*4/15 < p2 && p2 < p*8/15) {
						z_error_pos_211[C.getFaceQubits(S, face, direction, 3)][i] *= -1;
					} else if (p*8/15 < p2 && p2 < p*12/15) {
						z_error_pos_211[C.getFaceQubits(S, face, direction, 2)][i] *= -1;
						z_error_pos_211[C.getFaceQubits(S, face, direction, 3)][i] *= -1;
					}

					if (p3 < p*4/15) { //rose
						z_error_pos_211[C.getFaceQubits(S, face, direction, 1)][i] *= -1;
					} else if (p*4/15 < p3 && p3 < p*8/15) {
						z_error_pos_211[C.getFaceQubits(S, face, direction, 1)][i] *= -1;
						z_error_pos_211[C.getFaceQubits(S, face, direction, 2)][i] *= -1;
						z_error_pos_211[C.getFaceQubits(S, face, direction, 3)][i] *= -1;
					} else if (p*8/15 < p3 && p3 < p*12/15) {
						z_error_pos_211[C.getFaceQubits(S, face, direction, 2)][i] *= -1;
						z_error_pos_211[C.getFaceQubits(S, face, direction, 3)][i] *= -1;
					}
					
				}
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


void cluster::decode211(int verbosity = 0){
	//initialization loss and z_error
	
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		//stabalizer measurement (xx)
		int xx = 1; //stabalizer operator
		for (int k = 0; k < 2; k++) {
			xx *= z_error_pos_211[c][k];
		}
//		if (xx == -1) {//correction scheme
//			z_error_pos_211[c][0] *= -1;
//		}
		loss_pos[c] = xx;//discard qubit if bad
		
		//pauli measurement (xi)
		z_error_pos[c] = z_error_pos_211[c][0];//erronous qubit
	}
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		coord C(c, S);
		if(this_surf == PLANE && ((C.l == 1 && (C.y == S.y - 1 || C.x == 0)) ||(C.l == 2 && (C.z == S.z - 1 || C.x == 0)))){
			z_error_pos[c] = 1;
			loss_pos[c] = 1;
		}//	for planar code, remove errors on left, more, and lower boundaries to create edges.
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
		if (loss_pos[c] < 0) {//z_error_pos loss encountered
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

void cluster::getx_vec(){
	x_vec = {};
	for (int c = 0; c < S.x*S.y*S.z; c++){
		x_vec.push_back(1);//initialization
		//measurement of vertex operator
		vertex aVertex(c, S);
		for (int pos = 0; pos < 6; pos ++) {
			if (z_error_pos[aVertex.partial[pos]] == -1) {
				x_vec[c] *= -1;
			}
		}
	}
}

//void cluster::getSurf(){
//	for (int i =0; i < S.x; i++) {
//		for (int j =0; j < S.y; j++) {
//			for (int l=0; l < 2; l++) {
//				int bit = 1;
//				for (int k =0; k < S.z; k++) {
//					bit *= z_error_pos[coord(i, j, k, l).hash(S)];
//				}
//				surf[coord(i, j, 0, l).hash(S)] = bit;
//			}
//		}
//	}
//	cout << "here" <<endl;
//}
//
//void cluster::printSurf(){
//	vector<vector<string>> print_out;
//	for(int i = 0; i < 2*S.x; i++){
//		vector<string> a_line;
//		for(int j = 0; j < 2*S.y; j++){
//			a_line.push_back(" ");
//		}
//		print_out.push_back(a_line);
//	}
//	for (int x = 0; x < S.x; x++) {
//		for (int y = 0; y < S.y; y++) {
//			print_out[2 * y][2 * x + 1] = to_string(surf[coord(x, y, 0, 0).hash(S)])[0];
//			print_out[2 * y +1][2 * x] = to_string(surf[coord(x, y, 0, 1).hash(S)])[0];
//		}
//	}
//	cout << "surface:" << endl;
//	printMatrix(print_out);
//}

void cluster::print(){
	getx_vec();
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
					print_out[2 * y][2 * x + 1] = to_string(z_error_pos[c])[0];
					print_out[2 * y][2 * x] = to_string(x_vec[c])[0];
					print_out[2 * y +1][2 * x] = to_string(z_error_pos[coord(x, y, k/2, 1).hash(S)])[0];
				}
			}
		} else {//odd layer
			for (int x = 0; x < S.x; x++) {
				for (int y = 0; y < S.y; y++) {
					print_out[2 * y][2 * x] = to_string(z_error_pos[coord(x, y, k/2, 2).hash(S)])[0];
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
				a_line.push_back("   ");
			}
			print_out.push_back(a_line);
		}
		if (k % 2 == 0) {//even layer
			for (int x = 0; x < S.x; x++) {
				for (int y = 0; y < S.y; y++) {
					int c = coord(x, y, k/2, 0).hash(S);
					stringstream s1, s2, s3;
					s1 << setw(3) << to_string(findParent(c));
					s2 << setw(3) << to_string(loss_pos[c]);
					s3 << setw(3) << to_string(loss_pos[coord(x, y, k/2, 1).hash(S)]);
					print_out[2 * y][2 * x + 1] = s2.str();
					print_out[2 * y][2 * x] = s1.str();
					print_out[2 * y +1][2 * x] = s3.str();
				}
			}
		} else {//odd layer
			for (int x = 0; x < S.x; x++) {
				for (int y = 0; y < S.y; y++) {
					stringstream s1;
					s1 << setw(3) << to_string(loss_pos[coord(x, y, k/2, 2).hash(S)]);
					print_out[2 * y][2 * x] = s1.str();
				}
			}
		}
		cout << "layer:" << k << endl;
		printMatrix(print_out);
	}
}

//void cluster::surfaceCorrect(PerfectMatching* pm, const int& vertices_num, const vector<int>& wrong_chunks, const vector<int>& boundary_info, const int& verbosity = 0){
//	//Find the error operator of each pair.
//	vector<int> error_op_Z_pos; //correction operator not reduced
//	vector<int> matchPosition; //matching vertex of wrong_chunks respectively
//	for (int i = 0; i < vertices_num; i++) {
//		int ac = super_chunks[wrong_chunks[i]][0];
//		coord relative_pos; // relative position between two vertices in a pair
//		//get relative position (vector) of pair to determine the path in between.
//		vertex aVertex(ac, S), bVertex;
//		if (pm->GetMatch(i) < vertices_num){
//			int bc = super_chunks[wrong_chunks[pm->GetMatch(i)]][0];// position of vertexa's match, vertexb
//			matchPosition.push_back(bc);
//			if (count(matchPosition.begin(), matchPosition.end(), ac)) continue; //Prevent recounting of vertices
//			bVertex = vertex(bc, S);
//			relative_pos = getTaxicabDisplacement(S, ac, bc);
//		}else{ // matched to boundary
//			matchPosition.push_back(0);
//			int x =  coord(aVertex.c,S).x;
//			if (x <= S.x-x) {
//				relative_pos = {-x, 0, 0, 0};
//			} else {
//				relative_pos = {S.x-x, 0, 0, 0};
//			}
//		}
//		int n;
//		if (relative_pos.x > 0) {//a to the left of b
//			n = aVertex.partial[1];//use right z_error_pos
//			for(int i =0; i < abs(relative_pos.x); i++){
//				error_op_Z_pos.push_back(coord(divmod(coord(n,S).x + i, S.x), coord(n,S).y, 0, 0).hash(S));
//			}
//		} else if(relative_pos.x < 0) {//a to the right of b
//			n = aVertex.partial[0];//use left z_error_pos
//			for(int i =0; i < abs(relative_pos.x); i++){
//				error_op_Z_pos.push_back(coord(divmod(coord(n,S).x - i, S.x), coord(n,S).y, 0, 0).hash(S));
//			}
//		}
//		if (relative_pos.y > 0) {//a above b
//			n = bVertex.partial[2];//use upper z_error_pos
//			for(int i =0; i < abs(relative_pos.y); i++){
//				error_op_Z_pos.push_back(coord(coord(n,S).x, divmod(coord(n,S).y - i, S.y), 0, 1).hash(S));
//			}
//		} else if (relative_pos.y < 0) {//a below b
//			n = bVertex.partial[3];//use lower
//			for(int i =0; i < abs(relative_pos.y); i++){
//				error_op_Z_pos.push_back(coord(coord(n,S).x, divmod(coord(n,S).y + i, S.y), 0, 1).hash(S));
//			}
//		}
//
//	}
//
//	for (int i = 0; i < error_op_Z_pos.size(); i++) {
//		surf[error_op_Z_pos[i]] *= -1; //act local Z operator
//	}
//	if(verbosity == 2){
//		cout << "vertex/match/bn_info/: "<< endl;
//		//printout test
//		for (int i = 0 ; i < vertices_num; i++) {
//			for (int cV = 0; cV < super_chunks[wrong_chunks[i]].size(); cV ++) {
//				cout << coord(super_chunks[wrong_chunks[i]][cV],S);
//			}
//			cout << "/";
//			cout << coord(matchPosition[i],S)<<"/";
//			cout << boundary_info[i] <<"/";
//			cout << getTaxicabDistance(S,super_chunks[wrong_chunks[i]][0],matchPosition[i]) << endl;
//		}
//		cout << "errorOps: " << endl;
//		for (int i = 0 ; i <error_op_Z_pos.size(); i++) {
//			cout << coord(error_op_Z_pos[i],S) << " ";
//		}
//		cout << endl;
//	}
//}

int cluster::decodeWithMWPM(int verbosity = 0, int test = 0){
	//PAIR MATCHING
	getx_vec();
	vector<int> wrong_chunks;//actual (non boundary) chunks to be matched
	//e.g. [1,3,5,7]
	//get wrong_chunks
	for (int i = 0; i < super_chunks.size(); i++) {
		if (boundary_info[i] != 0) {// remove boundary chunks
			continue;
		}
		int measurement = 1;
		for (int j = 0; j < super_chunks[i].size(); j++) {
			measurement *= x_vec[super_chunks[i][j]];
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
//	surfaceCorrect(pm, vertices_num, wrong_chunks, boundary_info, verbosity);
	
	int parity = 1;
	for (int i = 0; i < vertices_num; i++) {
		if (pm->GetMatch(i) == i + vertices_num && boundary_info[wrong_chunks[i]] < 0) { //matched to left boundary
			parity *= -1;
		}
	}
	for (int i = 0; i < S.x; i++){//errors on left plane
		for (int j = 0; j < S.z; j++){
			parity *= z_error_pos[coord(0, i, j, 0).hash(S)];
		}
	}
	for (int i = 0; i < left_boundary_chunks.size(); i++) {
		for(int cV = 0; cV < super_chunks[left_boundary_chunks[i]].size(); cV++){
			parity *= x_vec[super_chunks[left_boundary_chunks[i]][cV]];
		}
	}
	
	delete pm;
	return parity;
}

void testDecoding(cluster& test_cluster, const double& p, const double& pl, const int& seed, surfacetype surf, noisemodel N, lossmodel Nl = INDEP_loss, int verbosity = 0){
	start_t = clock();
	//add noise
	if (N == GATE) {
		test_cluster.addGateNoise(p, seed);
	} else {
		test_cluster.addNoise(p, N, seed);
	}
	
	test_cluster.decode211(0);
	
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
//			test_cluster.getSurf();
//			test_cluster.printSurf();
		}
		cout << "t(getSuperChunks):" << double(clock()-start_t)/CLOCKS_PER_SEC << endl;//timing
		start_t = clock();
	}
	int parity = test_cluster.decodeWithMWPM(verbosity, 1);
	
	if (verbosity >= 1){
		if (verbosity == 2) {
//m
		}
		cout << "t(decodeWithMWPM):"<< double(clock()-start_t)/CLOCKS_PER_SEC << endl;//timing
		start_t = clock();
	}
	
	cout << "X1 check:" << parity <<endl;
}

int loopDecoding(const int Lmin, const int Lmax, const int trials, const double pmin, const double pmax, const int Np, const double plmin, const double plmax, const int Npl, const string fname, surfacetype surf, noisemodel N, lossmodel Nl, int verbosity = 0, int thread = 0){
	cout << "hardware_concurrency" << thread::hardware_concurrency() <<endl;
	ofstream outfile;
	outfile.open(fname);
	outfile << "L,p_error,num_success\n";
	int bins = 0 == Np ? 1:Np;
	int binsl = 0 == Npl ? 1:Npl;
	int num_correct;
	
	for (int il = 0; il <= Npl; il ++) {
		double pl = plmin + (plmax-plmin)/binsl * il;
		for (int L = Lmin; L <= Lmax; L = L+2) {
			int M = L;
			cout << endl;
			cluster test_cluster({L,L,L,0}, surf);
			for (int i = 0; i <= Np; i ++) {
				double p = pmin + (pmax-pmin)/bins * i;
				if (!thread) {
					cluster test_cluster({L,L,L,0}, surf);
				}
				
				//run simulation
				num_correct = 0;
				if (thread) {
					parallel_for(trials, [&](int start, int end){
						cluster test_cluster({L,L,L,0}, surf);
						for(int i = start; i < end; ++i){
							if (N == GATE) {
								test_cluster.addGateNoise(p, 0);
							} else {
								test_cluster.addNoise(p, N, 0);
							}
							test_cluster.decode211(0);
							try {
								test_cluster.getSuperChunks();
							} catch (...) {
								continue;
							}
							if (surf == PLANE && test_cluster.decodeWithMWPM(verbosity) == 1) {
								num_correct ++; //correction successful
							}
						}
					});
				} else {
					for(int i = 0; i < trials; ++i){
						if (N == GATE) {
							test_cluster.addGateNoise(p, 0);
						} else {
							test_cluster.addNoise(p, N, 0);
						}
						test_cluster.decode211(0);
						test_cluster.getSuperChunks();
						try {
							test_cluster.getSuperChunks();
						} catch (...) {
							cout << "h " <<endl;
							continue;
						}
						if (surf == PLANE && test_cluster.decodeWithMWPM(verbosity) == 1) {
							num_correct ++; //correction successful
						}
					}
				}
				
				//printout/outfile
				if (verbosity >= 2) {
					cout << "t=" << double(clock()-start_t)/CLOCKS_PER_SEC;
					start_t = clock();
				} if (verbosity >= 1){
					cout << L << "," << p << "," << pl << "," << num_correct << endl;
				} else {
					cout << ".";
					cout.flush();
				}
				outfile << L << "," << p << "," << pl << "," << num_correct << "\n";
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
	cout << "hardware_concurrency" << thread::hardware_concurrency() <<endl;
	string fname;
	surfacetype s;
	noisemodel N;
	lossmodel Nl;
	bool test, use_env, thread;
	int verbosity, return_value = 0;
	int Lmin, Lmax, n, Np, Npl, seed, times;
	float pmin, pmax, plmin, plmax;
	
	//getting options
	cxxopts::Options options(*argv,
							 "Simulator for fault-tolerant measurement-based quantum "
							 "computation on foliated surface code cluster states"
							 );
	options.add_options()
	("fname", "filename", cxxopts::value(fname)->default_value(""))
	("s", "surface type", cxxopts::value(s)->default_value("PLANE"))
	("Lmin", "Minimal size of mesh", cxxopts::value(Lmin)->default_value("3"))
	("Lmax", "Maximal size of mesh", cxxopts::value(Lmax)->default_value("17"))
	("n", "Number of trials", cxxopts::value(n)->default_value("100"))
	
	("N", "noise model", cxxopts::value(N)->default_value("GATE"))
	("Np", "z error p Points", cxxopts::value(Np)->default_value("10"))
	("pmin", "Minimal z error probability", cxxopts::value(pmin)->default_value("0.001"))
	("pmax", "Maximal z error probability", cxxopts::value(pmax)->default_value("0.008"))
	
	("Nl", "loss model", cxxopts::value(Nl)->default_value("INDEP_loss"))
	("Npl", "loss p points", cxxopts::value(Npl)->default_value("1"))
	("plmin", "Minimal loss probability", cxxopts::value(plmin)->default_value("0"))
	("plmax", "Maximal loss probability", cxxopts::value(plmax)->default_value("0"))
	
	("v, verbosity", "verbosity switch", cxxopts::value(verbosity)->default_value("0"))
	("seed", "seed switch", cxxopts::value(seed)->default_value("0"))
	("thread", "thread switch", cxxopts::value(thread)->default_value("0"))
	("use_env", "use environment variables", cxxopts::value(use_env)->default_value("0"))
	("times", "test times switch", cxxopts::value(times)->default_value("1"))
	("test", "test switch", cxxopts::value(test)->default_value("0"));
	options.parse(argc, argv);
	
	cout << "hardware_concurrency" << thread::hardware_concurrency() << endl;
	if (use_env) {
		cout << "SLURM_CPUS_PER_TASK:" << atoi(getenv("SLURM_CPUS_PER_TASK")) << endl;
	}
	//outputing options
	cout << "Lmin:" << Lmin << ";Lmax:" << Lmax << endl;
	cout << "Pmin:" << pmin << ";Pmax:" << pmax << ";nP" << Np << endl;
	cout << "Plmin:" << plmin << ";Plmax:" << plmax << ";nPl" << Npl << endl;
	cout << "n:" << n << ";seed:" << seed << ";thread" << thread << endl;
	cout << "N:" << N << ";verbosity:" << verbosity << "Nl;" << Nl << endl;
	
	
	if (fname == "") {
		fname = "L=" + to_string(Lmin) + ",P=(" + to_string(pmin).substr(3,2) + "," + to_string(pmax).substr(3,2) + "),n=" +to_string(n) + to_string(s) + "," + to_string(N) + ".out";
	}
	if (test) {
		cluster test_cluster({Lmin,Lmin,Lmin,0}, s);
		int i = 0;
		do {
			testDecoding(test_cluster, pmin, plmin, seed, s, N, Nl, verbosity);
			i++;
		}
		while (i < times);
	} else{
		return_value = loopDecoding(Lmin, Lmax, n, pmin, pmax, Np-1, plmin, plmax, Npl-1, fname, s, N, Nl, verbosity, thread);
	}
	return return_value;
}
///./simulate -s PLANE --pmin 0 --pmax 0.008  --Np 10 --Npl 1 -n 1000 --Lmin 3 -v 1

///./simulate -s PLANE --pmin 0.004 --pmax 0.008 --Np 10 -n 100 --Lmin 11 -v 1
///./simulate -N INDEP -s PLANE --pmin 0 -n 1000 --Lmin 17 -v 1 --test 1
///
//timing test
///./simulate -s PLANE --pmin 0.01 --pmax 0.05  --Np 10 --Npl 1 -n 500 --Lmin 3 -v 1

//one step large
	//(simple time)
	///./simulate -s PLANE --pmin 0.01 --pmax 0.05 --Np 10 --Npl 1 -n 500 --test 1 --Lmin 20 -v 2

	//(long time)
	//./simulate -s PLANE --plmin 0.05 --plmax 0.05 --pmin 0.004 --pmax 0.008  --Np 10 --Npl 1 -n 500 --Lmin 30 -v 1 --test 1

//loop run
	//(INDEP)
		///./simulate -N INDEP  --pmin 0.01 --pmax 0.05 --Lmin 3 -v 1 ###(pE ~0.04)
	//(GATE)
		///./simulate -N GATE --pmin 0 --pmax 0.012 --Lmin 3 -v 1 ###(pE ~0.0?)
