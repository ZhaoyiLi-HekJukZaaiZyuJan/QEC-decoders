//===================================================================//
//===================   ML Assisted decoding    =====================//
//=====================        main.cpp     =========================//
//===================================================================//
#include "main.hpp"
#include "thread_for.hpp"
#include <filesystem>
#include <iostream>

using namespace std;
namespace fs = filesystem;


static int window_size = 7;

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
		if (z_error_pos[c] > 0) {
			if (x_error_pos[c] > 0) {
				print_out[2 * y][2 * x + 1] = " ";
			} else {
				print_out[2 * y][2 * x + 1] = "Z";
			}
		} else {
			if (x_error_pos[c] > 0) {
				print_out[2 * y][2 * x + 1] = "X";
			} else {
				print_out[2 * y][2 * x + 1] = "Y";
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
				print_out[2 * y +1][2 * x] = "Z";
			}
		} else {
			if (x_error_pos[c] > 0) {
				print_out[2 * y +1][2 * x] = "X";
			} else {
				print_out[2 * y +1][2 * x] = "Y";
			}
		}
	}
	for(int c = 0; c < S.y*S.x; c++){//print x measurements
		int x = c % S.y; //x subcoordinate
		int y = c / S.y; //y subcoordinate
		if (stabs[c] > 0) {
			print_out[2 * y][2 * x] = "◯";
		} else {
			print_out[2 * y][2 * x] = "⊕";
		}
	}
	for(int c = 0; c < S.y*S.x; c++){//print z measurements
		int x = c % S.y; //x subcoordinate
		int y = c / S.y; //y subcoordinate
		if (stabs[c + S.y*S.x] > 0) {
			print_out[2 * y + 1][2 * x + 1] = "□";
		} else {
			print_out[2 * y + 1][2 * x + 1] = "⊟";
		}
	}
	printMatrix(print_out);
}

void subcluster::addNoise(const double& p, const int& seed=0){
	random_device rd;
	mt19937 engine{rd()};
	if (seed != 0) {
		engine.seed(seed);
	}
	uniform_real_distribution<> dist(0.0, 1.0);

	//uncorrelated error distribution for all physical qubits
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
			subcoord D(d);
			if(!D.l){
				window_int[d] = -(stabs[subcoord(divmod(C.x+D.x-window_size/2+1,S.x),divmod(C.y+D.y-window_size/2,S.y),0).hash()]-1)/2; //green
			} else {
				window_int[d] = (stabs[subcoord(divmod(C.x+D.x-window_size/2,S.x),divmod(C.y+D.y-window_size/2,S.y),1).hash()]-1)/2; //purple
			}
		}
	} else { //1 qubits
		for (int d = 0; d < 2*window_size*window_size; d++) {
			subcoord D(d);
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

void subcluster::decodeWithNN(cppflow::model model, int verbose = 0){
	for (int c = 0; c < 2*S.x*S.y; c++) {
		subcoord C(c);
		//get adjacent vertices
		int cV1 = subcoord(C.x,C.y,!C.l).hash(), cV2, cV3;
		if (!C.l) { //0 qubits
			cV2 = subcoord(divmod(C.x+1, S.x), C.y, 0).hash();
			cV3 = subcoord(C.x, divmod(C.y-1, S.y), 1).hash();
		} else {
			cV2 = subcoord(C.x, divmod(C.y+1, S.y), 0).hash();
			cV3 = subcoord(divmod(C.x-1, S.x), C.y, 1).hash();
		}
		if (stabs[c]<0||stabs[cV1]<0||stabs[cV2]<0||stabs[cV3]<0) {
			vector<float> window = getWindow(C);
			
			auto input = cppflow::tensor(window, {1,49*2});
			auto output = model(input);
			vector<int> corrections{2*int(round(output.get_data<float>()[0]))-1, 2*int(round(output.get_data<float>()[1]))-1};
			
			cout << window << endl;
			cout << z_error_pos[c] <<endl;
			cout << 2*int(round(output.get_data<float>()[0]))-1 << endl;
			
			z_error_pos[c] *= corrections[0];
			x_error_pos[c] *= corrections[1];
		}
	}
}
	

int subcluster::decodeWithMWPM(int verbose = 0){
	//PAIR MATCHING

    vector<int> subvertex_xPosition;//actual (non boundary) verteces to be matched
    for (int c = 0 ; c < stabs.size(); c++) {
        if (stabs[c] == -1) {
			if(this_surf == PLANE && c%S.y == 0){//only take the ones in the PLANE
				continue;
			}
            subvertex_xPosition.push_back(c);
        }
    }

	int vertices_num = subvertex_xPosition.size();
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
        int c1 = subvertex_xPosition[i];
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
            int c2 = subvertex_xPosition[j];
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

void testDecoding(cppflow::model model, const int L, const int M, const double p, const int seed, subsurfacetype surf=TORUS){
	subcluster testcluster({L,M,0},surf);
	testcluster.addNoise(p,seed);
	testcluster.getx_measurements();
	testcluster.getz_measurements();
	cout << "qubit: " << endl;
	testcluster.printQubit();
	
	testcluster.decodeWithNN(model);
	
	testcluster.getx_measurements();
	testcluster.getz_measurements();
	cout << "qubit: " << endl;
	testcluster.printQubit();

	cout << testcluster.decodeWithMWPM(1) << endl;
}

void loopDecoding(cppflow::model model, const int L_min, const int L_max, const int trials, const double P_min, const double P_max, const int Np, const string fname, subsurfacetype surf, bool verbose, bool thread, bool use_env){
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
			subcluster testcluster({L,L,0},surf);
			if (thread) {
				parallel_for(trials, [&](int start, int end){
					subcluster testcluster({L,L,0},surf);
					int num_add = 0;
					for(int i = start; i < end; ++i){
						testcluster.addNoise(p);
						testcluster.getx_measurements();
						
						if (surf == PLANE && testcluster.decodeWithMWPM(verbose) == 1) {
							num_correct  ++; //correction successful
						}
					}
					num_correct += num_add;
				}, use_env);
			} else {
				for (int i = 0; i < trials; i++) {
					testcluster.addNoise(p);
					testcluster.getx_measurements();
					
					if (surf == PLANE && testcluster.decodeWithMWPM(verbose) == 1) {
						num_correct  ++; //correction successful
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

void generator(const int trials, const double P_min, const double P_max, int seed, const int Np, const string fname, subsurfacetype surf, bool verbose, bool thread, bool use_env){
	ofstream outfile;
	outfile.open(fname);
	
	random_device rd;
	mt19937 engine{rd()};
	if (seed != 0) {
		engine.seed(seed);
	}
	vector<int> list(window_size*window_size*2);
	iota(list.begin(), list.end(), 1);
	discrete_distribution<> dist(list.begin(), list.end());
	int bins = 0 == Np ? 1 : Np;
	
	for (int i = 0; i <= Np; i ++) {
		if (i%(Np/10) == 0) {
			cout << ".";
			cout.flush();
		}
		
		double p = P_min + (P_max-P_min)/bins * i;
		
		int num_correct = 0;
		int L = window_size;
		subcluster testcluster({window_size,window_size,0},surf);
		
		for (int i = 0; i < trials; i++) {
			subcluster testcluster({window_size,window_size,0},surf);
			int num_add = 0;
			if (seed != 0) {
				testcluster.addNoise(p, seed + 1);
			} else {
				testcluster.addNoise(p, 0);
			}
			testcluster.getx_measurements();
			testcluster.getz_measurements();
			
			int c = dist(engine);
			vector<float> window = testcluster.getWindow(subcluster::subcoord(c));
			outfile << window << (testcluster.z_error_pos[c]+1)/2 << ", " <<  (testcluster.x_error_pos[c]+1)/2 << endl;
		}
	}
	cout << endl;
}

int main(int argc, const char *argv[]) {
	string fname;
	int L_min, L_max, verbose, n, Np, seed;
	bool test, thread, use_env, generate_mode;
	double P_max, P_min, test_probability;
	subsurfacetype s;
	cppflow::model model("/users/VanLadmon/OneDrive - Stanford/PHYSICS/Research/Patrick/ML/models/model2");//absolute directory to my_model

	cxxopts::Options options(*argv,"2D Simulator");
	options.add_options()
	("f, fname", "Output filename", cxxopts::value(fname))
	("s, surf_type", "Surface type", cxxopts::value(s)->default_value("PLANE"))
	("Lmin", "Minimal size of mesh", cxxopts::value(L_min)->default_value("3"))
	("Lmax", "Maximal size of mesh", cxxopts::value(L_max)->default_value("17"))
	("n", "Number of trials", cxxopts::value(n)->default_value("10000"))

	("v", "verbosity switch", cxxopts::value(verbose)->default_value("0"))
	("generate", "generation mode switch", cxxopts::value(generate_mode)->default_value("0"))
	("test", "test switch", cxxopts::value(test)->default_value("0"))

	("Np", "z error p Points", cxxopts::value(Np)->default_value("10"))
	("pmin", "Minimal z error probability", cxxopts::value(P_min)->default_value("0.01"))
	("pmax", "Maximal z error probability", cxxopts::value(P_max)->default_value("0"))

	("seed", "seed switch", cxxopts::value(seed)->default_value("0"))
	("thread", "thread switch", cxxopts::value(thread)->default_value("0"))
	("use_env", "use environment variables", cxxopts::value(use_env)->default_value("0"));
	options.parse(argc, argv);

	cout << "hardware_concurrency:" << thread::hardware_concurrency() << endl;
	if (use_env) {
		cout << "SLURM_CPUS_PER_TASK:" << atoi(getenv("SLURM_CPUS_PER_TASK")) << endl;
	}
	if (P_max < P_min){
		P_max = P_min;
	}
	cout << "Lmin:" << L_min << ";Lmax:" << L_max << endl;
	cout << "Pmin:" << P_min << ";Pmax:" << P_max << ";nP" << Np << endl;
	cout << "n:" << n << ";seed:" << seed << ";thread" << thread << endl;
	
	cout << fs::current_path() << endl;
	
	if (test) {
		testDecoding(model, L_min, L_min, P_min, seed, s);
	} else if (generate_mode) {
		if (fname == "") {
			fname = "/users/vanladmon/OneDrive - Stanford/PHYSICS/Research/Patrick/ML/train_data/train_set_L=" + to_string(window_size) + ",P=(" + to_string(P_min).substr(3,2) + "," + to_string(P_max).substr(3,2) + "),n=" +to_string(n) + ".out";
		}
		generator(n, P_min, P_max, seed, Np-1, fname, s ,verbose, thread, use_env);
	} else {
		if (fname == "") {
			fname = "L=" + to_string(L_min) + ",P=(" + to_string(P_min).substr(3,2) + "," + to_string(P_max).substr(3,2) + "),n=" +to_string(n) + to_string(s) + ".out";
		}
		loopDecoding(model, L_min, L_max, n, P_min, P_max, Np-1, fname, s ,verbose, thread, use_env);
	}
	return 1;
}


//-s PLANE --pmin 0 --pmax 0.008  --Np 20 -n 1000 --Lmin 3 -v 1 --thread -v 1
//./simulate -s TORUS --pmin 0.02 --pmax 0.02  --Np 20 --Lmin 10 --thread -v 1 --test
//./simulate -s TORUS --pmin 0.02 --pmax 0.02  --Np 20 -n 512 --Lmin 7 -v 1 --generate
