//########################## Faster Main function ##########################
//########################## 3D Fast ##########################

# include <iostream>
# include <fstream>
# include <random>
# include <vector>
# include <cstdlib>
# include <string>
# include <functional>
# include <algorithm>
# include <cmath>
# include <utility>
# include "libs/cxxopts/cxxopts.hpp"
# include "thread_for.hpp"
# include <assert.h>
# include <time.h>
# include "libs/blossom5-v2.05.src/PerfectMatching.h"
# include "libs/blossom5-v2.05.src/GEOM/GeomPerfectMatching.h"

using namespace std;

clock_t start_t = clock();
auto start = chrono::steady_clock::now();

int divmod (int a, int b) {return  (a%b+b) % b;} //modular division
enum surfacetype {PLANE};
istream& operator>> (istream& is, surfacetype& aSurf){
	string str;
	is >> str;
	if (str == "PLANE") aSurf = PLANE;
	else;
	return is;
}
string to_string(surfacetype& surf){
	if (surf == PLANE) return "PLANE";
}

enum noisemodel {INDEP, EM2, EM2_full, GATE};
istream& operator>> (istream& is, noisemodel& N){
	string str;
	is >> str;
	if (str == "INDEP") N = INDEP;
	else if (str == "EM2") N = EM2;
	else if (str == "EM2_full") N = EM2_full;
	else if (str == "GATE") N = GATE;
	else;
	return is;
}
string to_string(noisemodel& N){
	if (N == INDEP) return "INDEP";
	else if (N == EM2) return "EM2";
	else if (N == EM2_full) return "EM2_full";
	else if (N == GATE) return "GATE";
}


typedef function<pair<double, double>(double, double)> noisemodelfunc;

map<noisemodel, noisemodelfunc> NOISEMODELMAP {
	{EM2, [](const double& p,const double& q)->pair<double,double>{return {p*(32/15+2),p*4/15};}},
	{INDEP, [](const double& p,const double& q)->pair<double,double>{return {p, q};}},
	{EM2_full, [](const double& p,const double& q)->pair<double,double>{return {p*62/15 - pow(p,2)*1096/75 + pow(p,3)*32224/1125 - pow(p,4)*567296 /16875 + pow(p,5)*1196032/50625 - pow(p,6)*4194304/455625 + pow(p,7)*2097152/1366875, 1/2-sqrt(1/4-p*4/15)};}}
};

		
template <class T>
ostream& operator<<(ostream& os, const vector<T> vec) {
    os << "[";
    for (int i = 0; i < vec.size(); i++) {
        os << vec[i] << ", ";
    }
    os << "]";
    return os;
}

template <class T, class S>
ostream& operator<<(ostream& os, const pair<T,S> aPair) {
	os << "(" << aPair.first << "," << aPair.second << ")";
	return os;
}

template <class T>
void printMatrix(const vector<T> vec) {
	for (int i=0; i<vec.size(); i++) {
		cout << vec[i] << endl;
	}
}

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
	int getFaceQubits(const coord& S, const int& face, const int& k){
		if (face == 0) { // x-y plane
			switch(k) {
				case 0:
					return coord(x, divmod(y+1, S.y), z, 0).hash(S);
				case 1:
					return coord(x, y, z, 0).hash(S);
				case 2:
					return coord(divmod(x + 1, S.x), y, z, 1).hash(S);
				case 3:
					return hash(S);
			}
		} else if (face == 1) { // z-y plane
			switch(k) {
				case 0:
					return coord(x, y, z, 2).hash(S);//4
				case 1:
					return coord(x, divmod(y+1, S.y), z, 2).hash(S);//6
				case 2:
					return coord(x, divmod(y+1, S.y), divmod(z+1, S.z), 0).hash(S);//8
				case 3:
					return coord(x, divmod(y+1, S.y), z, 0).hash(S);//0
			}
		} else { //x-z plane
			switch(k) {
				case 0:
					return coord(x, divmod(y+1, S.y), z, 0).hash(S);//2
				case 1:
					return coord(x, y, divmod(z+1, S.z), 1).hash(S);//10
				case 2:
					return coord(divmod(x+1, S.x), y, z, 2).hash(S);//5
				case 3:
					return coord(x, y, z, 2).hash(S);//4
			}
		}
	}
};

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
	vertex(const int& hash, const coord& S){
		coord c(hash,S);
		this->c = c;
		physical_z_error_poss ={
			coord(divmod(c.x-1, S.x), c.y, c.z, 0).hash(S), coord(c.x, c.y, c.z, 0).hash(S),
			coord(c.x, divmod(c.y-1, S.y), c.z, 1).hash(S), coord(c.x, c.y, c.z, 1).hash(S),
			coord(c.x, c.y, divmod(c.z-1, S.z), 2).hash(S), coord(c.x, c.y, c.z, 2).hash(S),
		};
		//each vertex contains 6 physical_z_error_poss, left (Mx-), right(Mx+), up (Ly-), down (Ly+), less (Nz-), more (Nz+)
	};
	coord c;
	vector<int> physical_z_error_poss;
};

class cluster {
	vector<int> z_error_pos;
	coord S;
	surfacetype this_surf = PLANE;
	
	vector<int> surf;///debug
public:
	cluster(const coord&, const surfacetype&);
	void addNoise(const double&, const double&, const noisemodel, const int&);
	void addGateNoise(const double&, const noisemodel, const int&);
	vector<int> getZVec();
	void print();
	int decodeWithMWPM(int);
	
	void getSurf();
	void printSurf();
	void surfaceCorrect(PerfectMatching*, const int&, const vector<int>&, const vector<int>&, const int&);///debug
	
	int checkMeasurementOutcomeX1();
	int checkMeasurementOutcomeX2();
};

cluster::cluster(const coord& S, const surfacetype& this_surf){
	this->S = S;
	this->this_surf = this_surf;
	
	vector<int> surf(2*S.x*S.y*S.z,1);///debug
	this->surf = surf;///debug
	
	vector<int> z_error_pos(3*S.x*S.y*S.z,1);
	this->z_error_pos = z_error_pos;
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

	//Initialization of error operator
	//Error Model 1: uncorrelated error distribution for all physical z_error_poss
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		if(dist(engine) < NOISEMODELMAP[N](p, q).first) {
		z_error_pos[c] = -1;
		} else {
		z_error_pos[c] = 1;
		}
	}
	//Error Model 2: Correlation Error on the Same Edge (opposite faces) (check correctness)
	for (int c = 0; c < S.x*S.y*S.z; c++) {
		coord C(c,S);
		for (int face = 0; face < 3; face ++) {
			double p1 = dist(engine);
			double p2 = dist(engine);
			if (p1 < NOISEMODELMAP[N](p, q).second) { //vertical
				z_error_pos[C.getFaceQubits(S, face, 3)] *= -1;//3
				z_error_pos[C.getFaceQubits(S, face, 2)] *= -1;//2
			}
			if (p2 < NOISEMODELMAP[N](p, q).second) { //horizontal
				z_error_pos[C.getFaceQubits(S, face, 1)] *= -1;//1
				z_error_pos[C.getFaceQubits(S, face, 0)] *= -1;//0
			}
		}
	}

	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		coord C(c, S);
		if(this_surf == PLANE && ((C.l == 1 && (C.y == S.y - 1 || C.x == 0)) ||(C.l == 2 && (C.z == S.z - 1 || C.x == 0)))){
			z_error_pos[c] = 1;
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
	
	//Initialization of error operator
	//Error Model 1: uncorrelated error distribution for all physical z_error_poss
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		if(dist(engine) < p*16/15) {
			z_error_pos[c] = -1;
		} else {
			z_error_pos[c] = 1;
		}
	}
	
	for (int c = 0; c < S.x*S.y*S.z; c++) {
		coord C(c,S);
		for (int face = 0; face < 3; face ++) {
			double p1 = dist(engine);
			double p2 = dist(engine);
			if (p1 < p*4/15) { //vertical
				z_error_pos[C.getFaceQubits(S, face, 3)] *= -1;//3
			}
			if (p*4/15 < p1 < p*8/15) {
				z_error_pos[C.getFaceQubits(S, face, 2)] *= -1;//2
			} else if (p*8/15 < p1 < p*12/15){
				z_error_pos[C.getFaceQubits(S, face, 3)] *= -1;//3
				z_error_pos[C.getFaceQubits(S, face, 2)] *= -1;//2
			}
			if (p2 < p*4/15) { //horizontal
				z_error_pos[C.getFaceQubits(S, face, 1)] *= -1;//1
			}
			if (p*4/15 < p1 < p*8/15){
				z_error_pos[C.getFaceQubits(S, face, 0)] *= -1;//0
			} else if (p*8/15 < p1 < p*12/15){
				z_error_pos[C.getFaceQubits(S, face, 1)] *= -1;//1
				z_error_pos[C.getFaceQubits(S, face, 0)] *= -1;//0
			}
		}
	}
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		coord C(c, S);
		if(this_surf == PLANE && ((C.l == 1 && (C.y == S.y - 1 || C.x == 0)) ||(C.l == 2 && (C.z == S.z - 1 || C.x == 0)))){
			z_error_pos[c] = 1;
		}//	for planar code, remove errors on left, more, and lower boundaries to create edges.
	}
}

int getTaxicabDistance(const coord& S, const coord& c1, const coord& c2){  // compute taxicab distance between two coords
	return abs(c1.x - c2.x) + abs(c1.y - c2.y) + abs(c1.z - c2.z);
}
//
coord getTaxicabDisplacement(const coord& S, const coord& c1, const coord& c2){  //
	return {c2.x - c1.x, c2.y - c1.y, c2.z - c1.z, 0};
}


vector<int> cluster::getZVec(){
	vector<int> ZVec;
	for (int c = 0; c < S.x*S.y*S.z; c++){
		ZVec.push_back(1);
		//measurement of vertex operator
		vertex aVertex(c, S);
		for (int pos = 0; pos < 6; pos ++) {
			if (z_error_pos[aVertex.physical_z_error_poss[pos]] == -1) {
				ZVec[c] *= -1;
			}
		}
	}
	return ZVec;
}

void cluster::getSurf(){
	for (int i =0; i < S.x; i++) {
		for (int j =0; j < S.y; j++) {
			for (int l=0; l < 2; l++) {
				int bit = 1;
				for (int k =0; k < S.z; k++) {
					bit *= z_error_pos[coord(i, j, k, l).hash(S)];
				}
				surf[coord(i, j, 0, l).hash(S)] = bit;
			}
		}
	}
}

void cluster::printSurf(){
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
			print_out[2 * y][2 * x + 1] = to_string(surf[coord(x, y, 0, 0).hash(S)])[0];
			print_out[2 * y +1][2 * x] = to_string(surf[coord(x, y, 0, 1).hash(S)])[0];
		}
	}
	cout << "surface:" << endl;
	printMatrix(print_out);
}

void cluster::print(){
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
					print_out[2 * y][2 * x] = to_string(getZVec()[c])[0];
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

void cluster::surfaceCorrect(PerfectMatching* pm, const int& vertices_num, const vector<int>& vertexPosition, const vector<int>& boundary_nodes, const int& verbose = 0){
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
			relative_pos = getTaxicabDisplacement(S, coord(ac, S), coord(bc, S));
		}else{ // matched to boundary
			matchPosition.push_back(0);
			relative_pos = {boundary_nodes[i], 0, 0, 0};
		}
		int n;
		if (relative_pos.x > 0) {//a to the left of b
			n = aVertex.physical_z_error_poss[1];//use right z_error_pos
			for(int i =0; i < abs(relative_pos.x); i++){
				error_op_Z_pos.push_back(coord(divmod(coord(n,S).x + i, S.x), coord(n,S).y, 0, 0).hash(S));
			}
		} else if(relative_pos.x < 0) {//a to the right of b
			n = aVertex.physical_z_error_poss[0];//use left z_error_pos
			for(int i =0; i < abs(relative_pos.x); i++){
				error_op_Z_pos.push_back(coord(divmod(coord(n,S).x - i, S.x), coord(n,S).y, 0, 0).hash(S));
			}
		}
		if (relative_pos.y > 0) {//a above b
			n = bVertex.physical_z_error_poss[2];//use upper z_error_pos
			for(int i =0; i < abs(relative_pos.y); i++){
				error_op_Z_pos.push_back(coord(coord(n,S).x, divmod(coord(n,S).y - i, S.y), 0, 1).hash(S));
			}
		} else if (relative_pos.y < 0) {//a below b
			n = bVertex.physical_z_error_poss[3];//use lower z_error_pos
			for(int i =0; i < abs(relative_pos.y); i++){
				error_op_Z_pos.push_back(coord(coord(n,S).x, divmod(coord(n,S).y + i, S.y), 0, 1).hash(S));
			}
		}

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
			cout << boundary_nodes[i] << "/";
			cout << getTaxicabDistance(S,coord(vertexPosition[i],S),coord(matchPosition[i],S)) << endl;
		}
		cout << "errorOps: " << endl;
		for (int i = 0 ; i <error_op_Z_pos.size(); i++) {
			cout << coord(error_op_Z_pos[i],S) << " ";
		}
		cout<<endl;
	}
}

int cluster::decodeWithMWPM(int verbose = 0){
	//PAIR MATCHING
	vector<int> ZVec = getZVec();
    vector<int> vertexPosition;///actual (non boundary) vertices to be matched
	
    for (int c = 0; c < S.x*S.y*S.z; c++) {
        if (ZVec[c] == -1) {
			if(coord(c, S).x == 0){///for PLANAR code only: only take the vertices in the PLANE
				continue;
			}
            vertexPosition.push_back(c);
        }
    }
	
	int vertices_num = vertexPosition.size(), matches_num, edges_num;
	matches_num = 2 * vertices_num;
	edges_num = vertices_num * (vertices_num - 1) + vertices_num;

    PerfectMatching *pm = new PerfectMatching(matches_num, edges_num);
	pm->options.verbose = false;
	
	vector<int> boundary_nodes; //vector to keep track of x distance to boundary nodes

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
	for (int i = 0; i < vertices_num; i++) {
		coord c1(vertexPosition[i], S);
        for (int j = i + 1; j < vertices_num; j++) {// add in interconnections
            coord c2(vertexPosition[j],S);
			int dist = getTaxicabDistance(S, c1, c2);
			if (dist < abs(boundary_nodes[i]) + abs(boundary_nodes[j])) {//optimization
				pm->AddEdge(i,j,dist);
				pm->AddEdge(i+vertices_num,j+vertices_num,0);//boundary interconnection			}
//			pm->AddEdge(i,j,dist);
			}
        }
    }

	//solve the graph using MWPM decoder
	pm->Solve();
	
	//surface correction
	surfaceCorrect(pm, vertices_num, vertexPosition, boundary_nodes, verbose);///debug
	
	//alternative check
	int parity = 1;
	for (int i = 0; i < vertices_num; i++) {
		if (pm->GetMatch(i) == i + vertices_num && boundary_nodes[i] <= 0) {//matched to left boundary
			parity *= -1;
		}
	}
	for (int i = 0; i < S.x; i++){//errors on left boundary
		for (int j = 0; j < S.z; j++){
			parity *= z_error_pos[coord(0, i, j, 0).hash(S)];
		}
	}
	
	delete pm;
	return parity;
}

void testDecoding(cluster& test_cluster, const double& p, const int& seed, surfacetype surf = PLANE, noisemodel N = EM2, int verbose = 0){
	start_t = clock();
	if (N == GATE) {
		test_cluster.addGateNoise(p, N, seed);
	} else {
		test_cluster.addNoise(p, 0, N, seed);
	}
	if (verbose >= 1){
		if (verbose == 2) {
			test_cluster.print();
			test_cluster.getSurf();///debug
			test_cluster.printSurf();///debug
		}
		cout << "t(Generation):" << double(clock()-start_t)/CLOCKS_PER_SEC << endl;//timing
		start_t = clock();
	}
	
	int parity = test_cluster.decodeWithMWPM(verbose);
	
	if (verbose >= 1){
		if (verbose == 2) {
			test_cluster.printSurf();///debug
		}
		cout << "t(decodeWithMWPM):"<< double(clock()-start_t)/CLOCKS_PER_SEC << endl;//timing
		start_t = clock();
	}
	
	cout << "X1 check:" << parity <<endl;
}

int loopDecoding(const int L_min, const int L_max, const int trials, const double P_min, const double P_max, const int Np, const double Q_min, const double Q_max, const int Nq, const string fname, surfacetype surf = PLANE, noisemodel N = EM2, int verbose = 0, int thread = 1){
	ofstream outfile;
	outfile.open(fname);
	outfile << "L,p_error,num_success\n";
	int p_bins = 0 == Np ? 1 : Np;
	int q_bins = 0 == Nq ? 1 : Nq;
	int num_correct;

	for (int L = L_min; L <= L_max; L=L+2) {
		cluster test_cluster({L,L,L,0}, surf);
		int M = L;
		cout << endl;
		for (int iq = 0; iq <= Nq; iq ++) {
			double q = Q_min + (Q_max-Q_min)/q_bins * iq;
			
			for (int ip = 0; ip <= Np; ip ++) {
				double p = P_min + (P_max-P_min)/p_bins * ip;
				//run simulation
				int num_correct = 0;
				if (thread) {
					parallel_for(trials, [&](int start, int end){
						cluster test_cluster({L,L,L,0}, surf);
						for(int i = start; i < end; ++i){
							if (N == GATE) {
								test_cluster.addGateNoise(p, N);
							} else {
								test_cluster.addNoise(p, q, N);
							}
							if (surf == PLANE && test_cluster.decodeWithMWPM(verbose) == 1) {
								num_correct ++; //correction successful
							}
						}
					});
				} else {
					for(int i = 0; i < trials; ++i){
						if (N == GATE) {
							test_cluster.addGateNoise(p, N);
						} else {
							test_cluster.addNoise(p, q, N);
						}
						if (surf == PLANE && test_cluster.decodeWithMWPM(verbose) == 1) {
							num_correct ++; //correction successful
						}
					}
				}
				//printout/outfile
				if (verbose >= 2) {
					cout << "t=" << double(clock()-start_t)/CLOCKS_PER_SEC << endl;
					start_t = clock();
				} else if (verbose >= 1){
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
	string fname;
	surfacetype s;
	noisemodel N;
	bool test, use_env, thread;
	int L_min, L_max, n, Np, Nq, seed, times, verbose, return_value;
	float P_min, P_max, Q_min, Q_max;
	cxxopts::Options options(*argv,
							 "Simulator for fault-tolerant measurement-based quantum "
							 "computation on foliated surface code cluster states"
							 );
	options.add_options()
	("fname", "filename", cxxopts::value(fname)->default_value(""))
	("s", "surface type", cxxopts::value(s)->default_value("PLANE"))
	("Lmin", "Minimal size of mesh", cxxopts::value(L_min)->default_value("3"))
	("Lmax", "Maximal size of mesh", cxxopts::value(L_max)->default_value("17"))
	("n", "Number of trials", cxxopts::value(n)->default_value("10000"))
	
	("N", "noise model", cxxopts::value(N)->default_value("INDEP"))
	("Np", "single z_error_pos error Points", cxxopts::value(Np)->default_value("10"))
	("pmin", "Minimal single z_error_pos error probability", cxxopts::value(P_min)->default_value("0"))
	("pmax", "Maximal single z_error_pos error probability", cxxopts::value(P_max)->default_value("0.04"))
	
	("Nq", "correlated error points", cxxopts::value(Nq)->default_value("10"))
	("qmin", "Minimal correlated error probability", cxxopts::value(Q_min)->default_value("0"))
	("qmax", "Maximal correlated error probability", cxxopts::value(Q_max)->default_value("0.02"))
	
	("v", "verbosity switch", cxxopts::value(verbose)->default_value("0"))
	("seed", "seed switch", cxxopts::value(seed)->default_value("0"))
	("times", "test times switch", cxxopts::value(times))
	("thread", "thread switch", cxxopts::value(thread))
	("use_env", "use environment variables", cxxopts::value(use_env))
	("test", "test switch", cxxopts::value(test));
	
	options.parse(argc, argv);
	cout << "use_env:" << use_env << "thread" << thread << "test" << test << endl;
	cout << "Lmin:" << L_min << ";Lmax:" << L_max << endl;
	cout << "Pmin:" << P_min << ";Pmax:" << P_max << ";nP" << Np << endl;
	cout << "Qmin:" << Q_min << ";Qmax:" << Q_max << ";nQ" << Nq << endl;
	cout << "n:" << n << "seed:" << seed << ";thread" << thread << endl;
	
	//outputing options
	cout << "hardware_concurrency:" << thread::hardware_concurrency() << endl;
	if (use_env) {
		cout << "SLURM_CPUS_PER_TASK:" << atoi(getenv("SLURM_CPUS_PER_TASK")) << endl;
	}
	
	if (fname == "") {
		fname = "L=" + to_string(L_min) + ",P=(" + to_string(P_min).substr(3,2) + "," + to_string(P_max).substr(3,2) + "),n=" +to_string(n) + to_string(s) + "," + to_string(N) + ".out";
	}
	if (test) {
		cluster test_cluster({L_min,L_min,L_min,0}, s);
		int i = 0;
		do {
			testDecoding(test_cluster, P_min, seed, s, N, verbose);
			i++;
		}
		while (i < times);
	} else{
		loopDecoding(L_min, L_max, n, P_min, P_max, Np-1, Q_min, Q_max, Nq-1, fname, s, N, verbose, thread);
	}
	return return_value;
}
/// INDEP noise
/// ///./simulate -s PLANE --pmin 0 --pmax 0.05  --Np 20 -n 1000 --Lmin 3
///
///./simulate -s PLANE --pmin 0.004 --pmax 0.008  --Np 10 -n 500 --Lmin 15
///./simulate --pmin 0.004 --pmax 0.008 --Np 10 -n 100 --Lmin 11 -v 1
///./simulate -N INDEP -s PLANE --pmin 0 -n 1000 --Lmin 17 -v 1 --test 1
///./simulate  -s PLANE --pmin 0.004 -n 1000 --Lmin 17 -v 1 --test 1 --times 1 --seed 3



//large test
///./simulate -s PLANE --pmin 0 --pmax 0.04  --qmin 0 --qmax 0.04 --Np 10  --Nq 10 -n 500 --Lmin 3 -v 1
///./simulate  -s PLANE --pmin 0.004 -n 1000 --Lmin 17 -v 1 --test 1 --times 1 --seed 3
