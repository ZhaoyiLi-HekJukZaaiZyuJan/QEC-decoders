#include "submain.hpp"

using namespace std;

//===================================================================//

int divmod (int, int); //modular division

struct subcluster::subvertex {
	subvertex(){}
	subvertex(const int& hash){
		int L = S.x;
		int M = S.y;
		x = hash % L; //x subcoordinate
		y = hash / L; //y subcoordinate
		physical_qubits = {y*L+divmod(x-1,M), hash, M*L+divmod(y-1,L)*L+x, M*L+hash};
		//each subvertex contains 4 physical_qubits, left (Mx-), right (Mx+), up, down
	}
	int x;
	int y;
	vector<int> physical_qubits;
};

//===================================================================//
subcluster::subcoord subcluster::S = subcoord(0,0,0);

subcluster::subcoord::subcoord(const int& x, const int& y, const int& l){
	this->x = x;
	this->y = y;
	this->l = l;
}

subcluster::subcoord::subcoord(const int& c){
	int L = S.x;
	int M = S.y;
	*this = {c % (M*L) % L, c % (M*L) / L, c % (M*L)};
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
	os << "(" << c.x << "," << c.y << ","  << "," << c.l << ")";
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

int getTaxicabDistance(const subcluster::subcoord& S, const int& c1, const int& c2, const subsurfacetype& surf = subTORUS){  // compute taxicab distance between two cubes, given their cube numbers
	int L = S.x, M = S.y, x1 = subcluster::subcoord(c1).x, x2 = subcluster::subcoord(c2).x, y1 = subcluster::subcoord(c1).y, y2 = subcluster::subcoord(c2).y;

	if (surf == subTORUS) {
		return min(abs(x1 - x2), M-abs(x1 - x2)) + min(abs(y1 - y2), L-abs(y1 - y2));
	} else if (surf == subPLANE){
		return abs(x1 - x2) + abs(y1 - y2);
	} else{
		return 0;
	}
}

subcluster::subcoord getTaxicabDisplacement(const int& c1, const int& c2, const subsurfacetype& surf = subTORUS){  // compute taxicab distance between two points with given coordinates
	int L = subcluster::S.x, M = subcluster::S.y, x1 = subcluster::subcoord(c1).x, x2 = subcluster::subcoord(c2).x, y1 = subcluster::subcoord(c1).y, y2 = subcluster::subcoord(c2).y;

	if (surf == subTORUS) {
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

	} else if (surf == subPLANE){
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
	vector<int> x_vec(S.x*S.y,1);
	this->x_vec = x_vec;

}

void subcluster::printQubit(){
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
		print_out[2 * y][2 * x + 1] = to_string(z_error_pos[c])[0];
	}
	for(int c = S.y*S.x; c < 2*S.y*S.x; c++){ //print dual lattice
		int x = (c - S.y*S.x) % S.y; //x subcoordinate
		int y = (c - S.y*S.x) / S.y; //y subcoordinate
		print_out[2 * y +1][2 * x] = to_string(z_error_pos[c])[0];
	}
	for(int c = 0; c < S.y*S.x; c++){
		int x = c % S.y; //x subcoordinate
		int y = c / S.y; //y subcoordinate
		print_out[2 * y][2 * x] = to_string(x_vec[c])[0];
	}
	printMatrix(print_out);
}

void subcluster::addNoise(const double& p){
	random_device rd;
	mt19937 engine{rd()};
	uniform_real_distribution<> dist(0.0, 1.0);

	//uncorrelated error distribution for all physical qubits
	for (int i = 0; i < z_error_pos.size(); i++) {
		if(dist(engine) < p) {
			z_error_pos[i] = -1;
		} else {
			z_error_pos[i] = 1;
		}
	}
	if (this_surf == subPLANE) {//correction for subPLANE removal of boundary errors
		for (int i = 0; i < S.x; i++){
			z_error_pos[S.y*S.x + S.y*i] = 1;
		}
		for (int i = 0; i < S.y; i++){
			z_error_pos[2*S.y*S.x - i - 1] = 1;
		}
	}
}

void subcluster::getx_vec(){
	for (int c = 0; c < x_vec.size(); c++) {//measurement of subvertex operator
		subvertex asubvertex = subvertex(c);
		x_vec[c] = 1;
		for (int pos = 0; pos < 4; pos ++) {
			if (z_error_pos[asubvertex.physical_qubits[pos]] == -1) {
				x_vec[c] *= -1;
			}
		}
	}
}


int subcluster::decodeWithMWPM(int verbose = 0){
	//PAIR MATCHING

    vector<int> subvertexPosition;//actual (non boundary) verteces to be matched

    for (int c = 0 ; c < x_vec.size(); c++) {
        if (x_vec[c] == -1) {
			if(this_surf == subPLANE && c%S.y == 0){//only take the ones in the subPLANE
				continue;
			}
            subvertexPosition.push_back(c);
        }
    }

	int vertices_num = subvertexPosition.size();
	int matches_num;
	int edges_num;
	if (this_surf == subTORUS) {
		matches_num = vertices_num;
		edges_num = vertices_num * (vertices_num - 1)/2;
	} else if (this_surf == subPLANE){
		matches_num = 2 * vertices_num;
		edges_num = vertices_num * (vertices_num - 1) + vertices_num;
	}

    PerfectMatching *pm = new PerfectMatching(matches_num, edges_num);
    struct PerfectMatching::Options options;
    options.verbose = false;
    pm->options = options;

	vector<int> boundary_nodes; //vector to keep track of boundary nodes

    for (int i = 0; i < vertices_num; i++) {
        int c1 = subvertexPosition[i];
		if (this_surf == subPLANE) {// add in the boundary nodes
			int x1 = c1%S.y;
			int y1 = c1/S.y;
			if (x1 * 2 <= S.y) { // use left boundary
				boundary_nodes.push_back(-x1);
				pm->AddEdge(i, i + vertices_num, x1);
			} else{ //use right boundary
				boundary_nodes.push_back(S.y - x1);
				pm->AddEdge(i, i + vertices_num, S.y - x1);
			}
		}
        for (int j = i + 1; j < vertices_num; j++) {
            int c2 = subvertexPosition[j];
			pm->AddEdge(i,j,getTaxicabDistance(S, c1, c2, this_surf));
        }
    }
	if (this_surf == subPLANE) {//add in the interconnection of the boundary nodes
		for (int i = vertices_num; i < 2 * vertices_num; i++) {
			for (int j = i + 1; j < 2 * vertices_num; j++) {
				pm->AddEdge(i,j,0);
			}
		}
	}

	//solve the graph using MWPM decoder
	pm->Solve();

	int parity = 1;
	for (int i = 0; i < vertices_num; i++) {
		if (pm->GetMatch(i) == i + vertices_num && boundary_nodes[i] <= 0) { //matched to left boundary
			parity *= -1;
		}
	}
	for (int i = 0; i < S.y; i++){//errors on left boundary
		parity *= z_error_pos[S.x*i];
	}
	delete pm;
	return parity;
}

void testDecoding(const int L, const int M, const double p, subsurfacetype surf=subTORUS){
	subcluster testcluster({L,L,L},surf);
	testcluster.addNoise(p);
	testcluster.getx_vec();

	cout << "qubit: " << endl;
	testcluster.printQubit();
	

	cout << testcluster.decodeWithMWPM(1);

}

void loopDecoding(const int L_range, const int trials, const string fname, subsurfacetype surf=subTORUS, bool verbose = 0){
	ofstream outfile;
	outfile.open(fname);
	outfile << "L,p_error,num_success\n";
	for (int L = 3; L <= L_range; L=L+1) {
		if (verbose < 1) {
			cout << ".";
			cout.flush();
		} else cout << endl;
		
		for (float P = 0; P <= 15; P++) {
			float p =  P/100;
			int num_correct = 0;
			subcluster testcluster({L,L,L},surf);
			for (int i = 0; i < trials; i++) {

				testcluster.addNoise(p);
				testcluster.getx_vec();

				if (surf == subPLANE && testcluster.decodeWithMWPM(verbose) == 1) {
					num_correct  ++; //correction successful
				}

			}
			outfile << "{" << L << "," << p << "," << num_correct << "}, \n";
			if (verbose == 1) {
				cout << L << " " << p << " " << num_correct << "\n";
			}
		}
	}
}
