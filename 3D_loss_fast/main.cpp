//########################## 3D fast loss ##########################
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


void testDecoding(cluster& test_cluster, const double& p, const double& q, const int& seed, surfacetype surf, noisemodel N, lossmodel L, int verbose){
	start_t = clock();
	//add noise
	test_cluster.addNoise(p, q, N, L, seed);
	cout << "t(Generation):"<< double(clock()-start_t)/CLOCKS_PER_SEC << endl;//timing
	start_t = clock();
	if (L != OFF_loss){
		try {
			test_cluster.getSuperChunks();
		} catch (...) {
			cout << "failure" << endl;
			return;
		}
	}
	
	if (verbose >= 1){
		if (verbose == 2) {
			if (L != OFF_loss){
				test_cluster.printSuperChunks();
			}
			test_cluster.printPrimal(surf);
			test_cluster.getSurf();//decode
			test_cluster.printSurf();//decode
		}
		cout << "t(getSuperChunks):" << double(clock()-start_t)/CLOCKS_PER_SEC << endl;//timing
		start_t = clock();
	}
	int parity = test_cluster.decodeWithMWPM(verbose, 1, surf);

	if (verbose >= 1){
		if (verbose == 2) {
			test_cluster.printPrimal(surf);
			test_cluster.printSurf();//decode
		}
		cout << "t(decodeWithMWPMLoss):"<< double(clock()-start_t)/CLOCKS_PER_SEC << endl;//timing
		start_t = clock();
	}
	cout << "X1 check:" << parity <<endl;
}


int loopDecoding(const int lmin, const int lmax, const int trials, const double pmin, const double pmax, const int Np, const double qmin, const double qmax, const int Nq, const string fname, surfacetype surf, noisemodel N, lossmodel L, int verbose, bool thread, bool use_env, bool make_corrections){
	ofstream outfile;
	outfile.open(fname);
	outfile << "L,p_error,num_success\n";
	int bins = 0 == Np ? 1:Np;
	int binsl = 0 == Nq ? 1:Nq;
	int num_correct;
	
	for (int j = 0; j <= Nq; j++) {
		double q = qmin + (qmax-qmin)/binsl * j;
		for (int l = lmin; l <= lmax; l = l+2) {
			cout << endl;
			cluster test_cluster({l,l,l,0}, surf);
			for (int i = 0; i <= Np; i ++) {
				if (!thread) {
					cluster test_cluster({l,l,l,0}, surf);
				}
				double p = pmin + (pmax-pmin)/bins * i;
				//run simulation
				
				num_correct = 0;

				if (q == 0){  //use regular decoder
					if (thread) {
						parallel_for(trials, [&](int start, int end){
							cluster test_cluster({l,l,l,0}, surf);
							for(int i = start; i < end; ++i){
								test_cluster.addNoise(p, q, N, L);	
								if (test_cluster.decodeWithMWPM(verbose,0,surf) == 1) {
									num_correct ++; //correction successful
								}
							}
						}, use_env);
					} else {
						for(int i = 0; i < trials; ++i){
							test_cluster.addNoise(p, q, N, L);
							if (test_cluster.decodeWithMWPM(verbose,0,surf) == 1) {
								num_correct ++; //correction successful
							}
						}
					}
				} else { //use loss decoder
					if (thread) {
						parallel_for(trials, [&](int start, int end){
							cluster test_cluster({l,l,l,0}, surf);
							for(int i = start; i < end; ++i){
								test_cluster.addNoise(p, q, N, L);	
								try {
									test_cluster.getSuperChunks();
								} catch (...) {
									continue;
								}
								if (test_cluster.decodeWithMWPMLoss(verbose,0,surf) == 1) {
									num_correct ++; //correction successful
								}
							}
						}, use_env);
					} else {
						for(int i = 0; i < trials; ++i){
							test_cluster.addNoise(p, q, N, L);
							try {
								test_cluster.getSuperChunks();
							} catch (...) {
								continue;
							}
							if (test_cluster.decodeWithMWPMLoss(verbose,0,surf) == 1) {
								num_correct ++; //correction successful
							}
						}
					}
				}
				
				
				//printout/outfile
				if (verbose >= 2) {
					cout << "t=" << double(clock()-start_t)/CLOCKS_PER_SEC;
					start_t = clock();
				} if (verbose >= 1){
					cout << l << "," << p << "," << q << "," << num_correct << endl;
				} else {
					cout << ".";
					cout.flush();
				}
				outfile << l << "," << p << "," << q << "," << num_correct << "\n";
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
	cout << "s:" << s << endl;
	cout << "noise model:" << N << ";loss model:" << L <<endl;
	
	
	if (fname == "") {
		fname = "l=" + to_string(lmin) + ",p=(" + to_string(pmin).substr(3,2) + "," + to_string(pmax).substr(3,2) + "),n=" +to_string(n) + to_string(s) + "," + to_string(N) + ".out";
	}
	
	if (test) {
		cluster test_cluster({lmin,lmin,lmin,0}, s);
		int i = 0;
		do {
			testDecoding(test_cluster, pmin, qmin, seed, s, N, L, verbosity);
			i++;
		}
		while (i < times);
	} else{
		loopDecoding(lmin, lmax, n, pmin, pmax,  Np - 1, qmin, qmax, Nq - 1, fname, s, N, L, verbosity, thread, use_env, make_corrections);
	}
}

///######## big 2D run ########
///./simulate -s PLANE --qmin 0 --qmax 0.25 --pmin 0 --pmax 0.008  --Np 20 --Nq 20 -n 1000 --lmin 3 -v 1 #big loss test EM2
///./simulate -s PLANE --qmin 0 --qmax 0.25 --pmin 0 --pmax 0.03  --Np 20 --Nq 20 -n 100 --lmin 3 -v 1 -N INDEP  #big loss test INDEP
///./simulate -s PLANE --pmin 0 --pmax 0.012 --qmin 0 --qmax 0 --Np 10 --Nq 1  -n 1000 --lmin 3 -v 1 -N GATE --L toLoss #pauli error only

///######## 1D run ########
//### 211 run: (*p_th = 8.4%, most recent test: 51162075, note the slow convergence*) 
//./simulate -s PLANE --qmin 0 --qmax 0 --pmin 0.05 --pmax 0.1  --Np 20 --Nq 1 -n 1000 --lmin 3 -v 1 -N INDEP_211
//### INDEP run: (*p_th = 3%*)
///./simulate -s PLANE --qmin 0 --qmax 0 --pmin 0.02 --pmax 0.05  --Np 20 --Nq 1 -n 1000 --lmin 3 --lmax 17 -v 1 -N INDEP
//### GATE:      (*p_th = 0.585%, 6.5% verified: 53018107%*)
//### GATE_full: (*p_th = 0.591%, 7.6% verified: 14348%*) C.z == S.z - 1 || C.x == 0
//				 (*p_th = 0.57%, 3.3% verified: 54759532*) C.z == S.z - 1 || C.y == S.y - 1
//### EM2:       (*PLANE: p_th = 0.572%, 5.9% verified: 53112898%*) 
//				 (*TORUS:  *)
//./simulate -s TORUS --qmin 0 --qmax 0 --pmin 0 --pmax 0.009  --Np 20 --Nq 1 -n 10000 --lmin 3 --lmax 17 -v 1 -N EM2
//### EM2_full: 
//./simulate -s PLANE --qmin 0 --qmax 0 --pmin 0 --pmax 0.009  --Np 20 --Nq 1 -n 10000 --lmin 3 --lmax 17 -v 1 -N GATE

//### GATE_biased run \beta = 1000 (*p_ref = 0.74% ArXiv 1308.4776 *)
//./simulate -s PLANE --qmin 1000 --pmin 0.005 --pmax 0.01  --Np 30 --Nq 1 -n 10000 --lmin 3 --lmax 21 -v 1 -N GATE_biased



//### LOSS run: (*p_th = 25%*)
//### Loss EM2 run:
///./simulate -s PLANE --qmin 0 --qmax 0 --pmin 0 --pmax 0.008  --Np 20 --Nq 1 -n 10000 --lmin 3 -v 1 -N GATE


///./simulate -s PLANE --qmin 0 --qmax 0.3 --pmin 0 --pmax 0  --Np 1 --Nq 20 -n 1000 --lmin 3 --lmax 7 -v 1 -N INDEP --L toLoss

///######## tests ########
//timing test
///./simulate -s PLANE --qmin 0.05 --qmax 0.05 --pmin 0.004 --pmax 0.008  --Np 10 --Nq 1 -n 500 --lmin 3 -v 1

//loss test
///./simulate -N INDEP -s PLANE --pmin 0 --qmin 0 --lmin 7 -v 2 --test
