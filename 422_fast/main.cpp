//########################## [[4,2,2]] Toric Code ##########################
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
	vector<vector<int>> c_error_pos_422;
	public:
	void addPauli422(const double&, const double&, const int&seed = 0);
	void addNoise(const double&, const double&, const noisemodel, const lossmodel, const int& seed = 0);
	void addGateNoise(const double&, const int& seed = 0);
	Cluster(const coord&, const surfacetype&);
	void decode422(int);
};

Cluster::Cluster(const coord& S, const surfacetype& this_surf) : cluster(S, this_surf){ //this is needed for subclass initializers
	vector<vector<int>> c_error_pos_422(3*S.x*S.y*S.z, vector<int>(4, 1));
	this->c_error_pos_422 = c_error_pos_422;
}

void Cluster::addPauli422(const double & p, const double & q, const int& seed){
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
		for (int k = 0; k < 4; k++) {
			if(dist(engine) < p) {
				c_error_pos_422[c][k] = -1;
			} else {
				c_error_pos_422[c][k] = 1;
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
	//	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
	//		coord C(c, S);
	//		if(this_surf == PLANE && ((C.l == 1 && (C.y == S.y - 1 || C.x == 0)) ||(C.l == 2 && (C.z == S.z - 1 || C.x == 0)))){
	//			z_error_pos[c] = 1;
	//		}//	for planar code, remove errors on left, more, and lower boundaries to create edges.
	//	}
}

void Cluster::addGateNoise(const double & p, const int& seed){
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
		for (int k = 0; k < 4; k++) {
			for (int i = 0; i < 2; i++) { //p_S, p_M		
				if(dist(engine) < p*2/3) {
					c_error_pos_422[c][k] = -1;
				} else {
					c_error_pos_422[c][k] = 1;
				}
			}
			
		}
	}
	for (int c = 0; c < S.x*S.y*S.z; c++) {
		coord C(c,S);
		for (int i = 0; i < 4; i++) {//422
			for (int face = 0; face < 3; face ++) {
				double p1 = dist(engine); //process 1 (black)
				double p2 = dist(engine); //process 2 (pink)
				double p3 = dist(engine); //process 3 (rose)
				double p4 = dist(engine); //process 4 (purple)
				
				if (p1 < p*4/15) {
					c_error_pos_422[C.getFaceQubits(S, face, 0, 3)][i] *= -1;//3
				} else if (p*4/15 < p1 && p1 < p*8/15){
					c_error_pos_422[C.getFaceQubits(S, face, 0, 0)][i] *= -1;//0
					c_error_pos_422[C.getFaceQubits(S, face, 0, 1)][i] *= -1;//1
					c_error_pos_422[C.getFaceQubits(S, face, 0, 2)][i] *= -1;//2
				} else if (p*8/15 < p1 && p1 < p*12/15){
					c_error_pos_422[C.getFaceQubits(S, face, 0, 0)][i] *= -1;//0
					c_error_pos_422[C.getFaceQubits(S, face, 0, 1)][i] *= -1;//1
					c_error_pos_422[C.getFaceQubits(S, face, 0, 2)][i] *= -1;//2
					c_error_pos_422[C.getFaceQubits(S, face, 0, 3)][i] *= -1;//3
				}
				
				if (p2 < p*4/15) {
					c_error_pos_422[C.getFaceQubits(S, face, 0, 2)][i] *= -1;//2
				} else if (p*4/15 < p2 && p2 < p*8/15) {
					c_error_pos_422[C.getFaceQubits(S, face, 0, 0)][i] *= -1;//0
					c_error_pos_422[C.getFaceQubits(S, face, 0, 1)][i] *= -1;//1
				} else if (p*8/15 < p2 && p2 < p*12/15) {			
					c_error_pos_422[C.getFaceQubits(S, face, 0, 0)][i] *= -1;//0
					c_error_pos_422[C.getFaceQubits(S, face, 0, 1)][i] *= -1;//1
					c_error_pos_422[C.getFaceQubits(S, face, 0, 2)][i] *= -1;//2	
				}
				
				if (p3 < p*4/15) {
					c_error_pos_422[C.getFaceQubits(S, face, 0, 1)][i] *= -1;//1
				} else if (p*4/15 < p3 && p3 < p*8/15){
					c_error_pos_422[C.getFaceQubits(S, face, 0, 0)][i] *= -1;//0
				} else if (p*8/15 < p3 && p3 < p*12/15){
					c_error_pos_422[C.getFaceQubits(S, face, 0, 0)][i] *= -1;//0
					c_error_pos_422[C.getFaceQubits(S, face, 0, 1)][i] *= -1;//1
				}

				if (p4 < p*8/15) { //purple
					c_error_pos_422[C.getFaceQubits(S, face, 0, 0)][i] *= -1;//0
				}
			}
		}
			
		// //initialization errors X
		// for (int i = 0; i < 2; i++) {//422
		// 	for (int face = 0; face < 3; face ++) {
		// 		double pA1 = dist(engine); //process 1 (black)
		// 		double pA2 = dist(engine); //process 2 (pink)
				
		// 		if (pA1 < p*4/3*(1-p*2/3)){ // Z H errors
		// 			for(int k=0 ; k <4; k++){
		// 				c_error_pos_422[C.getFaceQubits(S, face, 0, k)][i+2] *= -1;//0123
		// 			}
		// 		}
				
		// 		if (pA2 < p*4/15){ // CZ Errors
		// 			for(int k=0 ; k <4; k++){
		// 				c_error_pos_422[C.getFaceQubits(S, face, 0, k)][i+2] *= -1;//0123
		// 			}
		// 		} else if (p*4/15 < pA2 < p*8/15){
		// 			for(int k=0 ; k <4; k++){
		// 				c_error_pos_422[C.getFaceQubits(S, face, 0, k)][i] *= -1;//0123
		// 			}
		// 		}
		// 	}
		// }
			
		// //initialization errors X
		// for (int k =0; k<4; k++) {
		// 	for (int i = 0; i < 3; i = i + 2) {
		// 		for (int face = 0; face < 3; face ++) {
		// 			double pA1 = dist(engine); //process 1 (black)
		// 			double pA2 = dist(engine); //process 2 (pink)
		// 			double pA3 = dist(engine); //process 2 (pink)
					
		// 			if (pA1 < p*2/3){ // Z H errors
		// 				c_error_pos_422[C.getFaceQubits(S, face, 0, k)][i+1] *= -1;//k
		// 			}
					
		// 			if (pA2 < p*2/3){ // Z H errors
		// 				c_error_pos_422[C.getFaceQubits(S, face, 0, k)][i] *= -1;//k
		// 			}
					
		// 			if (pA3 < p*4/15){ // CZ Errors
		// 				c_error_pos_422[C.getFaceQubits(S, face, 0, k)][i+1] *= -1;//k
		// 			} else if (p*4/15 < pA2 < p*8/15){
		// 				c_error_pos_422[C.getFaceQubits(S, face, 0, k)][i] *= -1;//k
		// 			}
		// 		}
		// 	}
		// }
	}
}

//No seed used for seed = 0
void Cluster::addNoise(const double & p, const double & q, const noisemodel N, const lossmodel L, const int& seed){
	if (N == GATE) { //p:error probability, //q: bias
		addGateNoise(p, seed);
	}
	if (N != GATE && N != GATE_biased) {
		addPauli422(NOISEMODELMAP[N](p,0).first, NOISEMODELMAP[N](p,0).second,  seed);
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

void Cluster::decode422(int verbose = 0){
	//initialization loss and z_error
	
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		//stabalizer measurement
		int xxxx = 1; //stabalizer operator
		
		for (int k = 0; k < 4; k++) {
			xxxx *= c_error_pos_422[c][k];
		}
		c_loss_pos[c] = xxxx;
		
		//pauli measurement
		int xxii = 1; //\bar{X} logical pauli operator
		//get z_error_pos
		for (int k = 0; k < 2; k++) {
			xxii *= c_error_pos_422[c][k];
		}
		c_error_pos[c] = xxii;
	}
	for (int c = 0; c < 3*S.x*S.y*S.z; c++) {
		coord C(c, S);
		if(this_surf == PLANE && ((C.l == 1 && (C.y == S.y - 1 || C.x == 0)) || (C.l == 2 && (C.z == S.z - 1 || C.x == 0)))){
			c_error_pos[c] = 1;
			c_loss_pos[c] = 1;
		}//	for planar code, remove errors on left, more, and lower boundaries to create edges.
	}
}

void testDecoding(Cluster& test_cluster, const double& p, const double& q, const int& seed, surfacetype surf, noisemodel N, lossmodel L, int verbosity){
	start_t = clock();
	//add noise
	if (N == GATE) {
		test_cluster.addGateNoise(p, seed);
	} else {
		test_cluster.addNoise(p, q, N, L, seed);
	}
	
	test_cluster.decode422();
	
	cout << "t(Generation):"<< double(clock()-start_t)/CLOCKS_PER_SEC << endl;//timing
	start_t = clock();
	try {
		test_cluster.getSuperChunks();
	} catch (...) {
		cout << "failure" << endl;
		return;
	}
	
	if (verbosity >= 1){
		test_cluster.printSuperChunks();
		test_cluster.getSurf();
		test_cluster.printSurf();
		cout << "t(getSuperChunks):" << double(clock()-start_t)/CLOCKS_PER_SEC << endl;//timing
		start_t = clock();
	}
	int parity = test_cluster.decodeWithMWPMLoss(verbosity, 1, surf);
	
	if (verbosity >= 1){
		test_cluster.printSurf();
		cout << "t(decodeWithMWPM):"<< double(clock()-start_t)/CLOCKS_PER_SEC << endl;//timing
		start_t = clock();
	}
	
	cout << "X1 check:" << parity <<endl;
}

int loopDecoding(const int lmin, const int lmax, const int trials, const double pmin, const double pmax, const int Np, const double qmin, const double qmax, const int Nq, const string fname, surfacetype surf, noisemodel N, lossmodel L, int verbosity, int thread, int out){
	cout << "hardware_concurrency" << thread::hardware_concurrency() <<endl;
	ofstream outfile;
	if(out){
		cout <<"here" << endl;
		outfile.open(fname);
		outfile << "L,p_error,num_success\n";
	}
	int binp = 0 == Np ? 1:Np;
	int binq = 0 == Nq ? 1:Nq;
	int num_correct;
	
	for (int j = 0; j <= Nq; j ++) {
		double q = qmin + (qmax-qmin)/binq * j;
		for (int l = lmin; l <= lmax; l = l+2) {
			int M = L;
			cout << endl;
			Cluster test_cluster({l,l,l,0}, surf);
			for (int i = 0; i <= Np; i ++) {
				double p = pmin + (pmax-pmin)/binp * i;
				if (!thread) {
					Cluster test_cluster({l,l,l,0}, surf);
				}
				
				//run simulation
				num_correct = 0;
				if (thread) {
					parallel_for(trials, [&](int start, int end){
						Cluster test_cluster({l,l,l,0}, surf);
						for(int i = start; i < end; ++i){
							if (N == GATE) {
								test_cluster.addGateNoise(p, 0);
							} else {
								test_cluster.addNoise(p, q, N, L, 0);
							}
							test_cluster.decode422();
							try {
								test_cluster.getSuperChunks();
							} catch (...) {
								continue;
							}
							if (surf == PLANE && test_cluster.decodeWithMWPM(verbosity, 0, surf) == 1) {
								num_correct ++; //correction successful
							}
						}
					});
				} else {
					for(int i = 0; i < trials; ++i){
						if (N == GATE) {
							test_cluster.addGateNoise(p, 0);
						} else {
							test_cluster.addNoise(p, q, N, L, 0);
						}
						test_cluster.decode422();
						try {
							test_cluster.getSuperChunks();
						} catch (...) {
							continue;
						}
						if (surf == PLANE && test_cluster.decodeWithMWPM(verbosity, 0, surf) == 1) {
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
				if(out){
					outfile << l << "," << p << "," << q << "," << num_correct << "\n";
				}		
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
	lossmodel L;
	bool test, use_env, thread, make_corrections, out;
	int verbosity, return_value = 0;
	int lmin, lmax, n, Np, Nq, seed, times;
	float pmin, pmax, qmin, qmax;
	
	//getting options
	cxxopts::Options options(*argv,
							 "Simulator for fault-tolerant measurement-based quantum "
							 "computation on foliated surface code cluster states"
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
	
	//outputing options
	cout << "hardware_concurrency:" << thread::hardware_concurrency() << endl;
	if (use_env) {
		cout << "SLURM_CPUS_PER_TASK:" << atoi(getenv("SLURM_CPUS_PER_TASK")) << endl;
	}
	cout << "lmin:" << lmin << ";lmax:" << lmax << endl;
	cout << "pmin:" << pmin << ";pmax:" << pmax << ";nP:" << Np << endl;
	cout << "qmin:" << qmin << ";qmax:" << qmax << ";nPl:" << Nq << endl;
	cout << "n:" << n << "seed:" << seed << ";thread:" << thread << endl;
	cout << "noise model:" << N << ";loss model:" << L <<endl;
	
	
	if (fname == "") {
		fname = "l=" + to_string(lmin) + ",P=(" + to_string(pmin).substr(3,2) + "," + to_string(pmax).substr(3,2) + "),n=" +to_string(n) + to_string(s) + "," + to_string(N) + ".out";
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
		return_value = loopDecoding(lmin, lmax, n, pmin, pmax, Np-1, qmin, qmax, Nq-1, fname, s, N, L, verbosity, thread, out);
	}
	return return_value;
}
///################ [[4,2,2]] ################

///######## 1D run ########

//### loop run:
	//(INDEP)  ###(p_th ~ 0.045 heuristics)
		///./simulate -N INDEP -n 10000 --pmin 0 --pmax 0.06 --lmin 3 --Np 30 -v 1
	//(GATE) ###(p_th ~ 0.0043, pl ~ 0.1)
		///./simulate -N GATE -n 10000 --pmin 0 --pmax 0.008 --lmin 3 -v 1
	//(GATE_biased (1)) run \beta = 1000 (*p_ref = *)
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
