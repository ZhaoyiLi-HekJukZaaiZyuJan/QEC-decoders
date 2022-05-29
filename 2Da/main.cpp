#include "main.hpp"
# include "thread_for.hpp"

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
	if (this_surf == PLANE) {//correction for PLANE removal of boundary errors
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
			if(this_surf == PLANE && c%S.y == 0){//only take the ones in the PLANE
				continue;
			}
            subvertexPosition.push_back(c);
        }
    }
	
	int vertices_num = subvertexPosition.size(), matches_num, edges_num;

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
        int c1 = subvertexPosition[i];
		if (this_surf == PLANE) {// add in the boundary nodes
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
	if (this_surf == PLANE) {//add in the interconnection of the boundary nodes
		for (int i = vertices_num; i < 2 * vertices_num; i++) {
			for (int j = i + 1; j < 2 * vertices_num; j++) {
				pm->AddEdge(i,j,0);
			}
		}
	}
	//solve the graph using MWPM decoder
	pm->Solve();

	int parity = 1;
	if (this_surf == PLANE){
		for (int i = 0; i < vertices_num; i++) {
			if (pm->GetMatch(i) == i + vertices_num && boundary_nodes[i] > 0) { //matched to right boundary
				parity *= -1;
			}
		} 
	}


	if (this_surf == TORUS){
		vector<int> matchPosition; //matching vertex of vertexPosition respectively
		for (int i = 0; i < vertices_num; i++) {
			
			int aC = subvertexPosition[i];
			int bC = subvertexPosition[pm->GetMatch(i)];// position of vertexa's match, vertexb

			matchPosition.push_back(bC);
			int x1 = aC % S.y;
			int x2 = bC % S.y;
			if (!((2*abs(x1 - x2) > S.x && x1 < x2))){  //crossed boundary
				parity *= -1;
			}
		}
	}

	for (int i = 0; i < S.y; i++){//errors on right boundary
		parity *= z_error_pos[S.x + S.x*i - 1];
	}
	

	delete pm;
	return parity;
}

void testDecoding(const int L, const int M, const double p, subsurfacetype surf=TORUS){
	subcluster testcluster({L,L,L},surf);
	testcluster.addNoise(p);
	testcluster.getx_vec();

	cout << "qubit: " << endl;
	testcluster.printQubit();
	

	cout << testcluster.decodeWithMWPM(1) << endl;

}

void loopDecoding(const int L_min, const int L_max, const int trials, const double P_min, const double P_max, const int Np, const string fname, subsurfacetype surf, bool verbose, bool thread, bool use_env){
	ofstream outfile;
	outfile.open(fname);
	outfile << "L,p_error,num_success\n";
	int bins = 0 == Np ? 1 : Np;
	for (int L = L_min; L <= L_max; L = L+2) {
		if (verbose < 1) {
			cout << ".";
			cout.flush();
		} else cout << endl;
		
		for (int i = 0; i <= Np; i ++) {
			double p = P_min + (P_max-P_min)/bins * i;
			
			int num_correct = 0;
			subcluster testcluster({L,L,L},surf);
			if (thread) {
				parallel_for(trials, [&](int start, int end){
					subcluster testcluster({L,L,L},surf);
					int num_add = 0;
					for(int i = start; i < end; ++i){
						testcluster.addNoise(p);
						testcluster.getx_vec();
						
						if (testcluster.decodeWithMWPM(verbose) == 1) {
							num_correct  ++; //correction successful
						}
					}
					num_correct += num_add;
				}, use_env);
			} else {
				for (int i = 0; i < trials; i++) {
					testcluster.addNoise(p);
					testcluster.getx_vec();
					
					if (testcluster.decodeWithMWPM(verbose) == 1) {
						num_correct ++; //correction successful
					}
				}
			}
			
			outfile << "{" << L << "," << p << "," << num_correct << "}, \n";
			if (verbose == 1) {
				cout << L << " " << p << " " << num_correct << "\n";
			}
		}
	}
}

int main(int argc, const char *argv[]) {
	string fname, code_model, str;
	int L_min, L_max, verbose, n, Np, seed;
	bool test, thread, use_env;
	double P_max, P_min;
	subsurfacetype s;

	cxxopts::Options options(*argv,"2D Simulator");

	options.add_options()
	("f, fname", "Output filename", cxxopts::value(fname))
	("s, surf_type", "Surface type", cxxopts::value(s)->default_value("PLANE"))
	("Lmin", "Minimal size of mesh", cxxopts::value(L_min)->default_value("3"))
	("Lmax", "Maximal size of mesh", cxxopts::value(L_max)->default_value("17"))
	("n", "Number of trials", cxxopts::value(n)->default_value("10000"))
	
	("v", "verbosity switch", cxxopts::value(verbose)->default_value("0"))

	("Np", "z error p Points", cxxopts::value(Np)->default_value("10"))
	("pmin", "Minimal z error probability", cxxopts::value(P_min)->default_value("0.001"))
	("pmax", "Maximal z error probability", cxxopts::value(P_max)->default_value("0.008"))
	
	("c", "code model", cxxopts::value(code_model)->default_value("2D"))
	("test", "test switch", cxxopts::value(test)->default_value("0"))
	("seed", "seed switch", cxxopts::value(seed)->default_value("0"))
	("thread", "thread switch", cxxopts::value(thread)->default_value("0"))
	("use_env", "use environment variables", cxxopts::value(use_env)->default_value("0"));
	options.parse(argc, argv);

	cout << "hardware_concurrency:" << thread::hardware_concurrency() << endl;
	if (use_env) {
		cout << "SLURM_CPUS_PER_TASK:" << atoi(getenv("SLURM_CPUS_PER_TASK")) << endl;
	}
	cout << "Lmin:" << L_min << ";Lmax:" << L_max << endl;
	cout << "Pmin:" << P_min << ";Pmax:" << P_max << ";nP" << Np << endl;
	cout << "n:" << n << ";seed:" << seed << ";thread" << thread << endl;

	if (fname == "") {
		fname = "L=" + to_string(L_min) + ",P=(" + to_string(P_min).substr(3,2) + "," + to_string(P_max).substr(3,2) + "),n=" +to_string(n) + to_string(s) + ".out";
	}
	
	if (test) {
		testDecoding(L_min, L_min, P_min, s);
	} else {
		loopDecoding(L_min, L_max, n, P_min, P_max, Np - 1, fname, s ,verbose, thread, use_env);
	}

}


//./simulate -s PLANE --pmin 0 --pmax 0.1 --Np 20 -n 10000 --Lmin 3 -v 1 --thread -v 1
