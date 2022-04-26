//===================================================================//
//===================   ML Assisted decoding (make)   =====================//
//=====================        main.cpp     =========================//
//===================================================================//
#include "main.hpp"
#include "thread_for.hpp"
#include <filesystem>
#include <iostream>

using namespace std;
namespace fs = filesystem;


static int window_size = 5;

//===================================================================//

int divmod (int, int); //modular division

struct subcluster::subvertex_x {
	subvertex_x(){}
	subvertex_x(const int& hash){
		int L = S.x;
		int M = S.y;
		x = hash % L; //x subcoordinate
		y = hash / L; //y subcoordinate
		physical_qubits = {y*L+divmod(x-1,M), hash, M*L+divmod(y-1,L)*L+x, M*L+hash};
		//each subvertex_x contains 4 physical_qubits, left (Mx-), right (Mx+), up, down
	}
	int x;
	int y;
	vector<int> physical_qubits;
};

struct subcluster::subvertex_z {//writing here
	subvertex_z(){}
	subvertex_z(const int& hash){
		int L = S.x;
		int M = S.y;
		x = hash % L; //x subcoordinate
		y = hash / L; //y subcoordinate
		physical_qubits = {M*L+hash, M*L+y*L+divmod(x+1,M), hash, divmod(y+1,L)*L+x};
		//each subvertex_x contains 4 physical_qubits, left (Mx-), right (Mx+), up, down
	}
	int x;
	int y;
	vector<int> physical_qubits;
};

//===================================================================//
subcluster::subcoord subcluster::S = subcoord(0,0,0); //global variable

subcluster::subcoord::subcoord(const int& x, const int& y, const int& l){
	this->x = x;
	this->y = y;
	this->l = l;
}

subcluster::subcoord::subcoord(const int& c){
	int L = S.x;
	int M = S.y;
	*this = {c % (M*L) % L, c % (M*L) / L, c / (M*L)};
}

subcluster::subcoord::subcoord(const int& c, const int& window_size){
	int L = window_size;
	int M = window_size;
	*this = {c % (M*L) % L, c % (M*L) / L, c / (M*L)};
}


int subcluster::subcoord::hash(){
	int L = S.x;
	int M = S.y;
	return l*L*M + y*L + x;
}

void subcluster::subcoord::operator=(const subcluster::subcoord& c) {
	x = c.x;
	y = c.y;
	l = c.l;
}


//print subcoord
ostream& operator<<(ostream& os, const subcluster::subcoord& c) {
	os << "(" << c.x << "," << c.y << "," << c.l << ")";
	return os;
}

//coord comparison
bool operator<(const subcluster::subcoord& c1, const subcluster::subcoord& c2) {
	if (c1.x < c2.x) {
		return true;
	} else if (c1.x > c2.x) {
		return false;
	} else {
		if (c1.y < c2.y){
			return true;
		} else {
			return false;
		}
	}
}

bool operator==(const subcluster::subcoord& c1, const subcluster::subcoord& c2) {
	if (c1.x == c2.x && c1.y == c2.y && c1.l == c2.l) return true;
	else return false;
}


//===================================================================//

int getTaxicabDistance(const subcluster::subcoord& S, const int& c1, const int& c2, const subsurfacetype& surf){  // compute taxicab distance between two cubes, given their cube numbers
	int L = S.x, M = S.y, x1 = subcluster::subcoord(c1).x, x2 = subcluster::subcoord(c2).x, y1 = subcluster::subcoord(c1).y, y2 = subcluster::subcoord(c2).y;
	
	if (surf == TORUS) {
		return min(abs(x1 - x2), M-abs(x1 - x2)) + min(abs(y1 - y2), L-abs(y1 - y2));
	} else if (surf == PLANE){
		return abs(x1 - x2) + abs(y1 - y2);
	} else{
		return 0;
	}
}

subcluster::subcoord getTaxicabDisplacement(const int& c1, const int& c2, const subsurfacetype& surf){  // compute taxicab distance between two points with given coordinates
	int L = subcluster::S.x, M = subcluster::S.y, x1 = subcluster::subcoord(c1).x, x2 = subcluster::subcoord(c2).x, y1 = subcluster::subcoord(c1).y, y2 = subcluster::subcoord(c2).y;
	
	if (surf == TORUS) {
		int x_dist, y_dist, x_relative_pos, y_relative_pos;
		x_dist = min(M-abs(x1 - x2), abs(x1 - x2));
		y_dist = min(L-abs(y1 - y2), abs(y1 - y2));
		
		//replace this part with displacement function
		if (x_dist == 0) { //two vertices on same x
			x_relative_pos = 0;
			
		} else if ((x_dist == abs(x1 - x2) && x1 < x2) ||
				   (x_dist == M - abs(x1 - x2) && x1 > x2)){
			x_relative_pos = 1;
		} else {
			x_relative_pos = -1;
		}
		
		if (y_dist == 0) { //two vertices on same y
			y_relative_pos = 0;
		} else if ((y_dist == abs(y1 - y2) && y1 < y2) ||
				   (y_dist == L - abs(y1 - y2) && y1 > y2)){
			y_relative_pos = 1;
			
		} else {
			y_relative_pos = -1;
		}
		return {x_relative_pos * x_dist, y_relative_pos * y_dist,0};
		
	} else if (surf == PLANE){
		return {x2 - x1, y2 - y1, 0};
	} else{
		return {0,0,0};
	}
}


subcluster::subcluster(const subcluster::subcoord& S, const subsurfacetype& this_surf){
	this->S = S;
	this->this_surf = this_surf;
	
	vector<int> z_error_pos(2*S.x*S.y,1);
	this->z_error_pos = z_error_pos;
	vector<int> x_error_pos(2*S.x*S.y,1);
	this->x_error_pos = x_error_pos;
	vector<int> stabs(2*S.x*S.y,1);
	this->stabs = stabs;
	
}

void subcluster::printQubit(int window_center = -1){
	cout << "qubit:" <<endl;
	vector<vector<string>> print_out;
	for(int i = 0; i < 2*S.x; i++){
		vector<string> a_line;
		for(int j = 0; j < 2*S.y; j++){
			a_line.push_back(" ");
		}
		print_out.push_back(a_line);
	}
	for(int c = 0; c < S.y*S.x; c++){ //print primal lattice
		int x = c % S.y; //x subcoordinate
		int y = c / S.y; //y subcoordinate
		if (z_error_pos[c] > 0) {
			if (x_error_pos[c] > 0) {
				print_out[2 * y][2 * x + 1] = " ";
			} else {
				print_out[2 * y][2 * x + 1] = "\033[1;33mX\033[0m";
			}
		} else {
			if (x_error_pos[c] > 0) {
				print_out[2 * y][2 * x + 1] = "\033[1;96mZ\033[0m";
			} else {
				print_out[2 * y][2 * x + 1] = "\033[1;32mY\033[0m";
			}
		}
	}
	for(int c = S.y*S.x; c < 2*S.y*S.x; c++){ //print dual lattice
		int x = (c - S.y*S.x) % S.y; //x subcoordinate
		int y = (c - S.y*S.x) / S.y; //y subcoordinate
		if (z_error_pos[c] > 0) {
			if (x_error_pos[c] > 0) {
				print_out[2 * y +1][2 * x] = " ";
			} else {
				print_out[2 * y +1][2 * x] = "\033[1;33mX\033[0m";
			}
		} else {
			if (x_error_pos[c] > 0) {
				print_out[2 * y +1][2 * x] = "\033[1;96mZ\033[0m";
			} else {
				print_out[2 * y +1][2 * x] = "\033[1;32mY\033[0m";
			}
		}
	}
	for(int c = 0; c < S.y*S.x; c++){//print x (green) measurements
		int x = c % S.y; //x subcoordinate
		int y = c / S.y; //y subcoordinate
		if (stabs[c] > 0) {
			print_out[2 * y][2 * x] = "◯";
		} else {
			print_out[2 * y][2 * x] = "\033[1;92m⊕\033[0m";
		}
	}
	for(int c = 0; c < S.y*S.x; c++){//print z measurements
		int x = c % S.y; //x subcoordinate
		int y = c / S.y; //y subcoordinate
		if (stabs[c + S.y*S.x] > 0) {
			print_out[2 * y + 1][2 * x + 1] = "□";
		} else {
			print_out[2 * y + 1][2 * x + 1] = "\033[1;95m⊟\033[0m";
		}
	}
	if (window_center != -1){
		subcoord C(window_center);
		if (!C.l) { //0 qubits
			for (int d = 0; d < 2*window_size*window_size; d++) {
				subcoord D(d,window_size);
				if(!D.l){
					print_out[2 * divmod(C.y+D.y-window_size/2,S.y)][2 * divmod(C.x+D.x-window_size/2+1,S.x)] = "W";
				} else {
					print_out[2 * divmod(C.y+D.y-window_size/2,S.y) + 1][2 * divmod(C.x+D.x-window_size/2,S.x) + 1] = "W"; //purple
				}
			}
		} else { //1 qubits
			for (int d = 0; d < 2*window_size*window_size; d++) {
				subcoord D(d,window_size);
				if(!D.l){
					print_out[2 * divmod(C.y+D.y-window_size/2,S.y) + 1][2 * divmod(C.x+D.x-window_size/2,S.x) + 1] = "W";//purple
				} else {
					print_out[2 * divmod(C.y+D.y-window_size/2+1,S.y)][2 * divmod(C.x+D.x-window_size/2,S.x)] = "W";
				}
			}
		}
	}
	
	printMatrix(print_out);
}

void subcluster::addError(){ //for debugging
	//Total heuristic Probabilities	
	//add errors manually here
	subcoord C1(3,3,0);
	z_error_pos[C1.hash()] *= -1;
	subcoord C2(3,1,1);
	z_error_pos[C2.hash()] *= -1;


	if (this_surf == PLANE) {//correction for PLANE removal of boundary errors
		for (int i = 0; i < S.x; i++){
			z_error_pos[S.y*S.x + S.y*i] = 1;
			x_error_pos[S.y*S.x + S.y*i] = 1;
		}
		for (int i = 0; i < S.y; i++){
			z_error_pos[2*S.y*S.x - i - 1] = 1;
			x_error_pos[2*S.y*S.x - i - 1] = 1;
		}
	}
}

void subcluster::addNoise(const double& p, const noisemodel N, const int& seed=0){
	random_device rd;
	mt19937 engine{rd()};
	if (seed != 0) {
		engine.seed(seed);
	}
	uniform_real_distribution<> dist(0.0, 1.0);
	
	//uncorrelated error distribution for all physical qubits
	if(N == INDEP){
		for (int i = 0; i < z_error_pos.size(); i++) {
			if(dist(engine) < p) {
				z_error_pos[i] = -1;
			} else {
				z_error_pos[i] = 1;
			}
		}
		for (int i = 0; i < x_error_pos.size(); i++) {
			if(dist(engine) < p) {
				x_error_pos[i] = -1;
			} else {
				x_error_pos[i] = 1;
			}
		}
	} else if(N == DEPOL){
		for (int i = 0; i < z_error_pos.size(); i++) {
			double r = dist(engine);
			if(r < p*2/3) {
				z_error_pos[i] = -1;
			} else {
				z_error_pos[i] = 1;
			}
			if(r > p*1/3 && r < p) {
				x_error_pos[i] = -1;
			}else {
				x_error_pos[i] = 1;
			}
		}
	}
	
	if (this_surf == PLANE) {//correction for PLANE removal of boundary errors
		for (int i = 0; i < S.x; i++){
			z_error_pos[S.y*S.x + S.y*i] = 1;
			x_error_pos[S.y*S.x + S.y*i] = 1;
		}
		for (int i = 0; i < S.y; i++){
			z_error_pos[2*S.y*S.x - i - 1] = 1;
			x_error_pos[2*S.y*S.x - i - 1] = 1;
		}
	}
}


void subcluster::getx_measurements(){
	for (int c = 0; c < S.x*S.y; c++) {//measurement of subvertex_x operator
		subvertex_x asubvertex = subvertex_x(c);
		stabs[c] = 1;
		for (int pos = 0; pos < 4; pos ++) {
			if (z_error_pos[asubvertex.physical_qubits[pos]] == -1) {
				stabs[c] *= -1;
			}
		}
	}
}

void subcluster::getz_measurements(){
	for (int c = 0; c < S.x*S.y; c++) {//measurement of subvertex_x operator
		subvertex_z asubvertex = subvertex_z(c);
		stabs[c+S.x*S.y] = 1;
		for (int pos = 0; pos < 4; pos ++) {
			if (x_error_pos[asubvertex.physical_qubits[pos]] == -1) {
				stabs[c+S.x*S.y] *= -1;
			}
		}
	}
}
vector<float> subcluster::getWindow(const subcoord& C){
	vector<int> window_int(2*window_size*window_size,0);
	if (!C.l) { //0 qubits
		for (int d = 0; d < 2*window_size*window_size; d++) {
			subcoord D(d,window_size);
			if(!D.l){
				window_int[d] = -(stabs[subcoord(divmod(C.x+D.x-window_size/2+1,S.x),divmod(C.y+D.y-window_size/2,S.y),0).hash()]-1)/2; //green
			} else {
				window_int[d] = (stabs[subcoord(divmod(C.x+D.x-window_size/2,S.x),divmod(C.y+D.y-window_size/2,S.y),1).hash()]-1)/2; //purple
			}
		}
	} else { //1 qubits
		for (int d = 0; d < 2*window_size*window_size; d++) {
			subcoord D(d,window_size);
			if(!D.l){
				window_int[d] = (stabs[subcoord(divmod(C.x+D.x-window_size/2,S.x),divmod(C.y+D.y-window_size/2,S.y),1).hash()]-1)/2; //purple
			} else {
				window_int[d] = -(stabs[subcoord(divmod(C.x+D.x-window_size/2,S.x),divmod(C.y+D.y-window_size/2+1,S.y),0).hash()]-1)/2; //green
			}
		}
	}
	vector<float> window(window_int.begin(), window_int.end());
	return window;
}

void subcluster::decodeWithNN(cppflow::model model, bool binary_output, int verbose = 0, double cutoff = 1){
	for (int c = 0; c < 2*S.x*S.y; c++) {
		subcoord C(c); //for l=0 qubits: left; for l=1 qubits: up(X stab l=0)
		//get adjacent vertices
		int cV1 = subcoord(C.x,C.y,!C.l).hash(), cV2, cV3; //for l=0 qubits: down; for l=1 qubits: right (Z stab l=1)
		if (!C.l) { //l=0 qubits
			cV2 = subcoord(divmod(C.x+1, S.x), C.y, 0).hash(); //right (X stab l=0)
			cV3 = subcoord(C.x, divmod(C.y-1, S.y), 1).hash(); //up (Z stab l=1)
		} else {    //l=1 qubits
			cV2 = subcoord(C.x, divmod(C.y+1, S.y), 0).hash(); //down (Z stab l=1)
			cV3 = subcoord(divmod(C.x-1, S.x), C.y, 1).hash(); //left
		}
		if (stabs[c]<0||stabs[cV1]<0||stabs[cV2]<0||stabs[cV3]<0) {
			vector<float> window = getWindow(C);
			
			auto input = cppflow::tensor(window, {1,window_size*window_size*2});
			auto output = model(input);
			// debug
			// cout << endl;
			// cout << C << endl;
			// cout << window << endl;
			// cout << output.get_data<float>() << endl;
			// printQubit(c);
			// printQubit();
			// cout << stabs << endl;
				
			
			if (binary_output) {
				// vector<int> corrections{2*int(round(output.get_data<float>()[0]))-1, 2*int(round(output.get_data<float>()[1]))-1};
				// z_error_pos[c] *= corrections[0];
				// x_error_pos[c] *= corrections[1];

				vector<int> corrections{2*int(round(1+cutoff*(output.get_data<float>()[0]-1)))-1, 1};
				z_error_pos[c] *= corrections[0];
			} else{
				vector<int> corrections{0,0};
				float max = 0;
				int max_pos = 0;
				for (int i = 0;i<4;i++){
					float this_value = output.get_data<float>()[i];			
					if (this_value > max) {
						// cout << this_value <<endl;
						max = this_value;
						max_pos = i;
					}
				}
				
				switch (max_pos) {
					case 0:
						corrections = {1,1};
						break;
					case 1:
						corrections = {1,-1};
						break;
					case 2:
						corrections = {-1,1};
						break;
					case 3:
						corrections = {-1,-1};
						break;
				}

				z_error_pos[c] *= corrections[0];
				x_error_pos[c] *= corrections[1];
				// getx_measurements();
				// getz_measurements();
			}
		}
	}
}

void subcluster::Predecode(bool binary_output, int verbose = 0){
	for (int c = 0; c < 2*S.x*S.y; c++) {
		subcoord C(c); //for l=0 qubits: left; for l=1 qubits: up(X stab l=0)
		//get adjacent vertices
		int cV1 = subcoord(C.x,C.y,!C.l).hash(), cV2, cV3; //for l=0 qubits: down; for l=1 qubits: right (Z stab l=1)
		if (!C.l) { //l=0 qubits
			cV2 = subcoord(divmod(C.x+1, S.x), C.y, 0).hash(); //right (X stab l=0)
			cV3 = subcoord(C.x, divmod(C.y-1, S.y), 1).hash(); //up (Z stab l=1)
		} else {    //l=1 qubits
			cV2 = subcoord(C.x, divmod(C.y+1, S.y), 0).hash(); //down (Z stab l=1)
			cV3 = subcoord(divmod(C.x-1, S.x), C.y, 1).hash(); //left
		}
		if (stabs[c]<0||stabs[cV1]<0||stabs[cV2]<0||stabs[cV3]<0) {
			vector<float> window = getWindow(C);
			
			auto input = cppflow::tensor(window, {1,window_size*window_size*2});
			auto output = model(input);
			// debug
			// cout << endl;
			// cout << C << endl;
			// cout << window << endl;
			// cout << output.get_data<float>() << endl;
			// printQubit(c);
			// printQubit();
			// cout << stabs << endl;
				
			
			if (binary_output) {
				// vector<int> corrections{2*int(round(output.get_data<float>()[0]))-1, 2*int(round(output.get_data<float>()[1]))-1};
				// z_error_pos[c] *= corrections[0];
				// x_error_pos[c] *= corrections[1];

				vector<int> corrections{2*int(round(1+cutoff*(output.get_data<float>()[0]-1)))-1, 1};
				z_error_pos[c] *= corrections[0];
			} else{
				vector<int> corrections{0,0};
				float max = 0;
				int max_pos = 0;
				for (int i = 0;i<4;i++){
					float this_value = output.get_data<float>()[i];			
					if (this_value > max) {
						// cout << this_value <<endl;
						max = this_value;
						max_pos = i;
					}
				}
				
				switch (max_pos) {
					case 0:
						corrections = {1,1};
						break;
					case 1:
						corrections = {1,-1};
						break;
					case 2:
						corrections = {-1,1};
						break;
					case 3:
						corrections = {-1,-1};
						break;
				}

				z_error_pos[c] *= corrections[0];
				x_error_pos[c] *= corrections[1];
				// getx_measurements();
				// getz_measurements();
			}
		}
	}
}

//direction : 0-z error 1-x error
vector<int> subcluster::decodeWithMWPM(int verbose = 0, bool dir = 0, bool make_corrections = 0){
	//PAIR MATCHING
	
	vector<int> vertexPosition;//actual (non boundary) vertices to be matched
	for (int c = 0 ; c < stabs.size()/2; c++) {
		int c_prime = c;
		if (dir){
			c_prime += stabs.size()/2;
		}
		if (stabs[c_prime] == -1) {
			if(this_surf == PLANE && c_prime  % S.y == 0){//only take the ones in the PLANE
				continue;
			}
			vertexPosition.push_back(c_prime);
		}
	}
	
	int vertices_num = vertexPosition.size();
	int matches_num;
	int edges_num;
	if (this_surf == TORUS) {
		matches_num = vertices_num;
		edges_num = vertices_num * (vertices_num - 1)/2;
	} else if (this_surf == PLANE){
		matches_num = 2 * vertices_num;
		edges_num = vertices_num * (vertices_num - 1) + vertices_num;
	}
	
	PerfectMatching *pm = new PerfectMatching(matches_num, edges_num);
	struct PerfectMatching::Options options;
	options.verbose = false;
	pm->options = options;
	
	vector<int> boundary_nodes; //vector to keep track of boundary nodes
	
	for (int i = 0; i < vertices_num; i++) {
		int c1 = vertexPosition[i];
		if (this_surf == PLANE) {// (x方向不同 需改)add in the boundary nodes
			int x1 = c1 %(S.x * S.y) % S.y;
			int y1 = c1 %(S.x * S.y) / S.y;
			if (x1 * 2 <= S.y) { // use left boundary
				boundary_nodes.push_back(-x1);
				pm->AddEdge(i, i + vertices_num, x1);
			} else{ //use right boundary
				boundary_nodes.push_back(S.y - x1);
				pm->AddEdge(i, i + vertices_num, S.y - x1);
			}
		}
		for (int j = i + 1; j < vertices_num; j++) {
			int c2 = vertexPosition[j];
			pm->AddEdge(i,j,getTaxicabDistance(S, c1, c2, this_surf));
		}
	}
	if (this_surf == PLANE) {//add in the interconnection of the boundary nodes
		for (int i = vertices_num; i < 2 * vertices_num; i++) {
			for (int j = i + 1; j < 2 * vertices_num; j++) {
				pm->AddEdge(i,j,0);
			}
		}
	}
	//solve the graph using MWPM decoder
	pm->Solve();
	
	vector<int>parity{1,1};
	
	if (this_surf == PLANE) {
		for (int i = 0; i < vertices_num; i++) {
			if (pm->GetMatch(i) == i + vertices_num && boundary_nodes[i] <= 0) { //matched to left boundary
				parity[1] *= -1;
			}
		}
		for (int i = 0; i < S.y; i++){//errors on left boundary
			parity[1] *= z_error_pos[S.x*i];
		}
	}
	
	if (make_corrections) {
		//======================= CORRECTION ==============================================
		//Find the error operator of each pair.
		vector<int> error_op_Z_pos; //
		vector<int> matchPosition; //matching vertex of vertexPosition respectively
		for (int i = 0; i < vertices_num; i++) {
			subcluster::subcoord relative_pos; //get relative position (vector) of pair to determine the path in between.
			if (dir == 0){
				int c1 = vertexPosition[i]; //position of vertexa
				subvertex_x aVertex(c1);
				subvertex_x bVertex;
				
				
				if (this_surf == TORUS || (this_surf == PLANE && pm->GetMatch(i) < vertices_num)){
					int c2 = vertexPosition[pm->GetMatch(i)];// position of vertexa's match, vertexb
					matchPosition.push_back(c2);
					if (count(matchPosition.begin(), matchPosition.end(), c1)) continue; //Prevent recounting of vertices
					bVertex = subvertex_x(c2);
					
					relative_pos = getTaxicabDisplacement(c1, c2, this_surf);
				}else{
					relative_pos = {boundary_nodes[pm->GetMatch(i) - vertices_num], 0,0};
				}
				int n, x, y;
				if (relative_pos.x > 0) {//a to the left of b
					n = aVertex.physical_qubits[1];//use right qubit
					x = n % S.x;
					y = n / S.x;
					for(int i = 0; i<abs(relative_pos.x); i++){
						error_op_Z_pos.push_back(divmod(x+i, S.x) + y*S.y);
					}
				} else if(relative_pos.x < 0) {//a to the right of b
					n = aVertex.physical_qubits[0];//use left qubit
					x = n % S.x;
					y = n / S.x;
					for(int i =0; i<abs(relative_pos.x); i++){
						error_op_Z_pos.push_back(divmod(x-i, S.x) + y*S.y);
					}
				}		
				if (relative_pos.y > 0) {//a above b
					n = bVertex.physical_qubits[2];//use upper qubit
					x = (n - S.y*S.x) % S.x;
					y = (n - S.y*S.x) / S.x;
					for(int i =0; i< abs(relative_pos.y); i++){
						error_op_Z_pos.push_back(x + divmod(y-i,S.y)*S.y + S.x*S.y);
					}
				} else if (relative_pos.y < 0) {//a below b
					n = bVertex.physical_qubits[3];//use lower qubit
					x = (n - S.y*S.x) % S.x;
					y = (n - S.y*S.x) / S.x;
					for(int i =0; i< abs(relative_pos.y); i++){
						error_op_Z_pos.push_back(x + divmod(y+i,S.y)*S.y + S.x*S.y);
					}
				}
			} else {
				int c1 = vertexPosition[i]; //position of vertexa
				subvertex_z aVertex(c1 - S.x*S.y);
				subvertex_z bVertex;
				
				if (this_surf == TORUS || (this_surf == PLANE && pm->GetMatch(i) < vertices_num)){
					int c2 = vertexPosition[pm->GetMatch(i)];// position of vertexa's match, vertexb
					matchPosition.push_back(c2);
					if (count(matchPosition.begin(), matchPosition.end(), c1)) continue; //Prevent recounting of vertices
					bVertex = subvertex_z(c2 - S.x*S.y);
					
					relative_pos = getTaxicabDisplacement(c1, c2, this_surf);
				}else{
					relative_pos = {boundary_nodes[pm->GetMatch(i) - vertices_num], 0,0};
				}
				int n, x, y;
				if (relative_pos.x > 0) {//a to the left of b
					n = aVertex.physical_qubits[1];//use right qubit
					x = (n - S.x*S.y) % S.x;
					y = (n - S.x*S.y) / S.x;
					for(int i = 0; i<abs(relative_pos.x); i++){
						error_op_Z_pos.push_back(divmod(x+i, S.x) + y*S.y + S.x*S.y);
					}
				} else if(relative_pos.x < 0) {//a to the right of b
					n = aVertex.physical_qubits[0];//use left qubit
					
					x = (n - S.x*S.y) % S.x;
					y = (n - S.x*S.y) / S.x;
					for(int i =0; i<abs(relative_pos.x); i++){
						error_op_Z_pos.push_back(divmod(x-i, S.x) + y*S.y + S.x*S.y);
					}
				}
				if (relative_pos.y > 0) {//a above b
					n = bVertex.physical_qubits[2];//use upper qubit
					x = n % S.x;
					y = n / S.x;
					for(int i =0; i< abs(relative_pos.y); i++){
						error_op_Z_pos.push_back(x + divmod(y-i,S.y)*S.y);
					}
				} else if (relative_pos.y < 0) {//a below b
					n = bVertex.physical_qubits[3];//use lower qubit
					x = n % S.x;
					y = n / S.x;
					for(int i =0; i< abs(relative_pos.y); i++){
						error_op_Z_pos.push_back(x + divmod(y+i,S.y)*S.y);
					}
				}
			}
		}
		//remove repeated Z in correction op
		vector<int> error_op_Z_pos_reduced;
		for (int i = 0; i < error_op_Z_pos.size(); i++) {
			int c = error_op_Z_pos[i];
			if (count(error_op_Z_pos_reduced.begin(), error_op_Z_pos_reduced.end(), c)) continue;
			error_op_Z_pos_reduced.push_back(c);
			if (!dir){
				z_error_pos[c] *= -1; //act local Z operator
			} else{
				x_error_pos[c] *= -1;
			}	
		}
	}
	
	if(this_surf == TORUS){
		int result = 1;
		if (!dir){
			for (int i = 0; i < S.x ; i++) { // Check X1 operator
				result *= z_error_pos[i * S.y];
			}
			parity[0] = result;
			
			result = 1;
			for (int i = S.x*S.y; i < (S.y+1)*S.x ; i++) { // Check X2 operator
				result *= z_error_pos[i];
			}
			parity[1] = result;
		} else{
			for (int i = 0; i < S.x ; i++) { // Check Z1 operator
				result *= x_error_pos[i];
			}
			parity[0] = result;
			
			result = 1;
			for (int i = 0; i < S.x ; i++) { // Check Z2 operator
				result *= x_error_pos[i * S.y + S.x*S.y];
			}
			parity[1] = result;
		}
		
	}
	
	delete pm;
	return parity;
}

void testDecoding(cppflow::model model, const int L, const int M, const double p, const int seed, bool binary_output, bool make_corrections, bool decode_with_NN, subsurfacetype surf, noisemodel N, bool dir, float cutoff, bool verbose){
	subcluster testcluster({L,M,0},surf);
	if (N == MANUAL){
		testcluster.addError();
	} else {
		testcluster.addNoise(p, N, seed);
	}
	
	testcluster.getx_measurements();
	testcluster.getz_measurements();
	testcluster.printQubit();
	
	testcluster.decodeWithNN(model, binary_output, verbose, cutoff);
	testcluster.getx_measurements();
	testcluster.getz_measurements();
	testcluster.printQubit();
	
	if (surf == PLANE && testcluster.decodeWithMWPM(1,dir,make_corrections)[0] == 1) {
		cout << "success" << endl; //correction successful
	} else if (surf == TORUS) {
		vector<int> parity = testcluster.decodeWithMWPM(1,dir,make_corrections);
		if (parity[0] == 1 && parity[1] == 1) {
			cout << "success" << endl;
		} else {
			cout << "fail" << endl;
		}
	} else {
		cout << "fail" << endl;
	}
	testcluster.getx_measurements();
	testcluster.getz_measurements();
	cout << "qubit: " << endl;
	testcluster.printQubit();
}



void loopDecoding(string directory, string model_name, const int L_min, const int L_max, const int trials, const double P_min, const double P_max, const int Np, const string fname, bool binary_output, subsurfacetype surf, noisemodel N, bool dir, bool verbose, bool thread, bool use_env, bool decode_with_NN){
	ofstream outfile;
	outfile.open(fname);
	outfile << "L,p_error,num_success\n";
	int bins = 0 == Np ? 1 : Np;
	for (int L = L_min; L <= L_max; L = L+2) {
		if (verbose < 1) {
			cout << ".";
			cout.flush();
		} else cout << endl;

		for (int i = 1; i <= Np; i ++) {
			float p = P_min + (P_max-P_min)/(bins+1) * i;
			float p_rounded = round(p*10000)/10000;
			string p_name = to_string(p);
			while (p_name.back()=='0'){ //remove insignificant figures
				p_name.pop_back();
			}

			int num_correct = 0;
			subcluster testcluster({L,L,0},surf);

			if (decode_with_NN){
				cppflow::model model(directory + "models/" + model_name + p_name);//absolute directory to my_model
				if (thread) {
					parallel_for(trials, [&](int start, int end){
						subcluster testcluster({L,L,0},surf);
						int num_add = 0;
						for(int i = start; i < end; ++i){
							testcluster.addNoise(p, N);
							testcluster.getx_measurements();
							testcluster.getz_measurements();

							testcluster.decodeWithNN(model, binary_output);
							testcluster.getx_measurements();
							testcluster.getz_measurements();
							
							if (surf == PLANE && testcluster.decodeWithMWPM(verbose, dir, 0)[0] == 1) {
								num_add++; //correction successful
							} else{
								vector<int> parity = testcluster.decodeWithMWPM(verbose, dir, 1);
								if (parity[0]==1&&parity[1]==1) {
									num_add++;
								}
							}
						}
						num_correct += num_add;
					}, use_env);
				} else {
					for (int i = 0; i < trials; i++) {
						testcluster.addNoise(p, N);
						testcluster.getx_measurements();
						testcluster.getz_measurements();

						testcluster.decodeWithNN(model, binary_output);
						testcluster.getx_measurements();
						testcluster.getz_measurements();
						
						if (surf == PLANE && testcluster.decodeWithMWPM(verbose, dir)[0] == 1) {
							num_correct ++; //correction successful
						} else{
								vector<int> parity = testcluster.decodeWithMWPM(verbose, dir, 1);
								if (parity[0]==1&&parity[1]==1) {
									num_correct++;
								}
							}
					}
				}
			} else {
				if (thread) {
					parallel_for(trials, [&](int start, int end){
						subcluster testcluster({L,L,0},surf);
						int num_add = 0;
						for(int i = start; i < end; ++i){
							testcluster.addNoise(p, N);
							testcluster.getx_measurements();
							testcluster.getz_measurements();
							
							if (surf == PLANE && testcluster.decodeWithMWPM(verbose,dir, 0)[0] == 1) {
								num_add++; //correction successful
							} else{
								vector<int> parity = testcluster.decodeWithMWPM(verbose, dir, 1);
								if (parity[0]==1&&parity[1]==1) {
									num_add++;
								}
							}
						}
						num_correct += num_add;
					}, use_env);
				} else {
					for (int i = 0; i < trials; i++) {
						testcluster.addNoise(p, N);
						testcluster.getx_measurements();
						testcluster.getz_measurements();
						
						if (surf == PLANE && testcluster.decodeWithMWPM(verbose, dir)[0] == 1) {
							num_correct ++; //correction successful
						} else{
								vector<int> parity = testcluster.decodeWithMWPM(verbose, dir, 1);
								if (parity[0]==1&&parity[1]==1) {
									num_correct++;
								}
							}
					}
				}
			}
			
			
			outfile << L << "," << p << "," << num_correct << "\n";
			if (verbose == 1) {
				cout << L << "," << p << "," << num_correct << "\n";
			}
		}
	}
}

void generator(const int trials, const int L_min, const double P_min, const double P_max, int seed, const int Np, const string fname, subsurfacetype surf, noisemodel N, bool verbose, bool thread, bool use_env, bool binary_output, bool new_cluster){
	ofstream outfile;
	outfile.open(fname);
	
	random_device rd;
	mt19937 engine{rd()};
	if (seed != 0) {
		engine.seed(seed);
	}
	float p = P_min;
	
	vector<int> list(window_size*window_size*2);
	
	for (int i = 0; i < trials; i++) {
		if (verbose > 0 && i%(trials/10) == 0) {
			cout << ".";
			cout.flush();
		}
		while(true){
			subcluster testcluster({L_min,L_min,0},surf);
			int num_add = 0;
			if (seed != 0) {
				testcluster.addNoise(p, N, seed + 1);
			} else {
				testcluster.addNoise(p, N, 0);
			}
			testcluster.getx_measurements();
			testcluster.getz_measurements();
			//generate method 1
			//randomly pick an erronous stabilizer measurement, random pick one out of 4 neighbors.
			vector<int> test_pos;
			// for (int c = 0; c < 2*subcluster::S.x*subcluster::S.y; c++) {
			// 	if(testcluster.stabs[c] < 0){
			// 		test_pos.push_back(c);
			// 	}
			// }

			//generate method 2
			//only pick x measurements
			for (int c = 0; c < subcluster::S.x*subcluster::S.y; c++) {
				if(testcluster.stabs[c] < 0){
					test_pos.push_back(c);
				}
			}
			int size = test_pos.size();
			if (size == 0){
				continue;
			} else {
				uniform_int_distribution<int> dist(0,size-1);
				int h = test_pos[dist(engine)];
				int c;
				uniform_int_distribution<int> dist2(0,3);
				if (h < subcluster::S.x*subcluster::S.y){
					subcluster::subvertex_x aVertex(h);
					c = aVertex.physical_qubits[dist2(engine)];
				} else {
					subcluster::subvertex_z aVertex(h-subcluster::S.x*subcluster::S.y);
					c = aVertex.physical_qubits[dist2(engine)];
				}
				
				subcluster::subcoord C(c);
				// testcluster.printQubit();
			    // testcluster.printQubit(c);

				vector<float> window = testcluster.getWindow(C);

				if (binary_output) {
					outfile << window << (testcluster.z_error_pos[c]+1)/2 << ", " <<  (testcluster.x_error_pos[c]+1)/2 << endl;
				} else{
					int choice;
					if ((testcluster.z_error_pos[c]+1)/2) {
						if ((testcluster.x_error_pos[c]+1)/2) {
							choice = 0;
						} else {
							choice = 1;
						}
					} else {
						if ((testcluster.x_error_pos[c]+1)/2) {
							choice = 2;
						} else {
							choice = 3;
						}
					}
					outfile << window << choice << endl;
				}
				i++;
			}
			if (i == trials) {
				return;
			}
		}
	}
}

int main(int argc, const char *argv[]) {
	string fname, directory, model_name;
	int L_min, L_max, verbose, n, Np, seed, rpt, dir;
	bool test, thread, use_env, generate_mode, binary_output, make_corrections, new_cluster, decode_with_NN;
	double P_max, P_min, cutoff;
	subsurfacetype s;
	noisemodel N;
	
	cxxopts::Options options(*argv,"2D Simulator");
	options.add_options()
	("f, fname", "Output filename", cxxopts::value(fname))
	("d, directory", "model directory", cxxopts::value(directory)->default_value("/users/VanLadmon/OneDrive - Stanford/PHYSICS/Research/Patrick/ML make/"))
	("m, model", "model parent name", cxxopts::value(model_name)->default_value("model,L=5(7),layer=5x512,epochs=1000,p="))
	("s, surf_type", "Surface type", cxxopts::value(s)->default_value("TORUS"))
	("Lmin", "Minimal size of mesh", cxxopts::value(L_min)->default_value("3"))
	("Lmax", "Maximal size of mesh", cxxopts::value(L_max)->default_value("17"))
	("n", "Number of trials", cxxopts::value(n)->default_value("10000"))
	("dir", "decoding direction", cxxopts::value(dir)->default_value("0"))
	
	("v", "verbosity switch", cxxopts::value(verbose)->default_value("0"))
	("generate", "generation mode switch", cxxopts::value(generate_mode)->default_value("0"))
	("binary", "binary data format switch", cxxopts::value(binary_output)->default_value("0"))
	("make_corrections", "make corrections to result", cxxopts::value(make_corrections)->default_value("0"))
	("decode_with_NN", "turn on NN decoder", cxxopts::value(decode_with_NN)->default_value("0"))
	("cutoff", "NN decoder acceptance cutoff", cxxopts::value(cutoff)->default_value("1"))
	("new", "new cluster", cxxopts::value(new_cluster)->default_value("0"))
	
	("N", "noise model", cxxopts::value(N)->default_value("DEPOL"))
	("Np", "z error p Points", cxxopts::value(Np)->default_value("10"))
	("pmin", "Minimal z error probability", cxxopts::value(P_min)->default_value("0.01"))
	("pmax", "Maximal z error probability", cxxopts::value(P_max)->default_value("0"))
	
	("test", "test switch", cxxopts::value(test)->default_value("0"))
	("seed", "seed switch", cxxopts::value(seed)->default_value("0"))
	("rpt", "test repeat", cxxopts::value(rpt)->default_value("1"))
	("thread", "thread switch", cxxopts::value(thread)->default_value("0"))
	("use_env", "use environment variables", cxxopts::value(use_env)->default_value("0"));
	options.parse(argc, argv);
	
	if (verbose > 0) {
		cout << "hardware_concurrency:" << thread::hardware_concurrency() << endl;
		if (use_env) {
			cout << "SLURM_CPUS_PER_TASK:" << atoi(getenv("SLURM_CPUS_PER_TASK")) << endl;
		}
		if (P_max < P_min){
			P_max = P_min;
		}
		cout << "Lmin:" << L_min << ";Lmax:" << L_max << endl;
		cout << "Pmin:" << P_min << ";Pmax:" << P_max << ";nP" << Np << endl;
		cout << "n:" << n << ";seed:" << seed << ";thread" << thread << ";make_corrections" << make_corrections <<  ";decode_w_NN:" << decode_with_NN << endl;
		cout << "dir:"<< dir << endl;
		cout << "cutoff" << cutoff << endl;
		
		cout << fs::current_path() << endl;
	}
	
	if (test) {
		cppflow::model model(directory + "models/" + model_name);//absolute directory to my_model
		for(int i = 0; i < rpt; i++){
			testDecoding(model, L_min, L_min, P_min, seed + i, binary_output, make_corrections, decode_with_NN, s, N, dir, cutoff, verbose);
		}
	} else if (generate_mode) {
		if (fname == "") {
			fname = "/users/vanladmon/OneDrive - Stanford/PHYSICS/Research/Patrick/ML make/train_data/train_set_L=" + to_string(window_size) + ",P=(" + to_string(P_min).substr(3,2) + "," + to_string(P_max).substr(3,2) + "),n=" +to_string(n) + ".out";
		}
		generator(n, L_min, P_min, P_max, seed, Np-1, fname, s, N, verbose, thread, use_env, binary_output, new_cluster);
	} else {
		if (fname == "") {
			fname = "L=" + to_string(L_min) + ",P=(" + to_string(P_min).substr(3,2) + "," + to_string(P_max).substr(3,2) + "),n=" +to_string(n) + to_string(s) + ".out";
		}
		loopDecoding(directory, model_name, L_min, L_max, n, P_min, P_max, Np-1, fname, binary_output, s, N, dir, verbose, thread, use_env, decode_with_NN);
	}
	return 1;
}

//./simulate -s TORUS --pmin 0 --pmax 0.18  --Np 25 -n 1000 --Lmin 5 --Lmax 5 -v 1 -d ~/ML/ -m "model,L=5(7),layer=3x128,epochs=10000,p=" --decode_with_NN
//./simulate -s TORUS --pmin 0 --pmax 0.12  --Np 25 -n 1000 --Lmin 3 --Lmax 20 -v 1 -d ~/ML/  --fname test.out

//./simulate -s TORUS --pmin 0.036 --Np 10 --Lmin 10 -v 1 --test --make_corrections -d /scratch/users/ladmon/ML/ -m "model_h,L=5(7),layer=3x128,epochs=100000,p=0.036" --binary
//./simulate -s TORUS --pmin 0.02 --pmax 0.02  --Np 20 -n 1 --Lmin 7 -v 1 --generate -d ~/ML

//model,L=5(7),layer=5x512,epochs=1000,p=0.1068