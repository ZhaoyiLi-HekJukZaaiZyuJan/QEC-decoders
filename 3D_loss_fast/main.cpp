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


void testDecoding(cluster& test_cluster, const double& p, const double& pl, const int& seed, surfacetype surf, noisemodel N, lossmodel Nl, int verbose){
	start_t = clock();
	//add noise
	test_cluster.addNoise(p, pl, N, Nl, seed);
	cout << "t(Generation):"<< double(clock()-start_t)/CLOCKS_PER_SEC << endl;//timing
	start_t = clock();
	try {
		test_cluster.getSuperChunks();
	} catch (...) {
		cout << "failure" << endl;
		return;
	}
	
	if (verbose >= 1){
		if (verbose == 2) {
			test_cluster.printSuperChunks();
			test_cluster.printPrimal(surf);
			test_cluster.getSurf();//decode
			test_cluster.printSurf();//decode
		}
		cout << "t(getSuperChunks):" << double(clock()-start_t)/CLOCKS_PER_SEC << endl;//timing
		start_t = clock();
	}
	int parity = test_cluster.decodeWithMWPMLoss(verbose, 1, 1, surf);

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


int loopDecoding(const int L_min, const int L_max, const int trials, const double P_min, const double P_max, const int Np, const double Pl_min, const double Pl_max, const int Npl, const string fname, surfacetype surf, noisemodel N, lossmodel Nl, int verbose, bool thread, bool use_env, bool make_corrections){
	ofstream outfile;
	outfile.open(fname);
	outfile << "L,p_error,num_success\n";
	int bins = 0 == Np ? 1:Np;
	int binsl = 0 == Npl ? 1:Npl;
	int num_correct;
	
	for (int il = 0; il <= Npl; il ++) {
		double pl = Pl_min + (Pl_max-Pl_min)/binsl * il;
		for (int L = L_min; L <= L_max; L = L+2) {
			int M = L;
			cout << endl;
			cluster test_cluster({L,L,L,0}, surf);
			for (int i = 0; i <= Np; i ++) {
				if (!thread) {
					cluster test_cluster({L,L,L,0}, surf);
				}
				double p = P_min + (P_max-P_min)/bins * i;
				//run simulation
				
				num_correct = 0;
				if (thread) {
					parallel_for(trials, [&](int start, int end){
						cluster test_cluster({L,L,L,0}, surf);
						for(int i = start; i < end; ++i){
							test_cluster.addNoise(p, pl, N, Nl);	
							try {
								test_cluster.getSuperChunks();
							} catch (...) {
								continue;
							}
							if (surf == PLANE && test_cluster.decodeWithMWPMLoss(verbose,0,0,surf) == 1) {
								num_correct ++; //correction successful
							}
						}
					}, use_env);
				} else {
					for(int i = 0; i < trials; ++i){
						test_cluster.addNoise(p, pl, N, Nl);
						try {
							test_cluster.getSuperChunks();
						} catch (...) {
							continue;
						}
						if (surf == PLANE && test_cluster.decodeWithMWPMLoss(verbose,0,0,surf) == 1) {
							num_correct ++; //correction successful
						}
					}
				}
				
				//printout/outfile
				if (verbose >= 2) {
					cout << "t=" << double(clock()-start_t)/CLOCKS_PER_SEC;
					start_t = clock();
				} if (verbose >= 1){
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
	string fname;
	surfacetype s;
	noisemodel N;
	lossmodel Nl;
	bool test, use_env, thread, make_corrections;
	int verbose;
	int L_min, L_max, n, Np, Npl, seed, times;
	float P_min, P_max, Pl_min, Pl_max;
	cxxopts::Options options(*argv,
							 "Simulator for fault-tolerant measurement-based quantum "
							 "computation on foliated surface code cluster states"
							 );
	options.add_options()
	("fname", "filename", cxxopts::value(fname)->default_value(""))
	("s", "surface type", cxxopts::value(s)->default_value("PLANE"))
	("lmin", "Minimal size of mesh", cxxopts::value(L_min)->default_value("3"))
	("lmax", "Maximal size of mesh", cxxopts::value(L_max)->default_value("17"))
	("n", "Number of trials", cxxopts::value(n)->default_value("10000"))
	
	("N", "noise model", cxxopts::value(N)->default_value("EM2"))
	("Np", "z error p Points", cxxopts::value(Np)->default_value("10"))
	("pmin", "Minimal z error probability", cxxopts::value(P_min)->default_value("0.000"))
	("pmax", "Maximal z error probability", cxxopts::value(P_max)->default_value("0.008"))
	
	("Nl", "loss model", cxxopts::value(Nl)->default_value("OFF_loss"))
	("Npl", "loss p points", cxxopts::value(Npl)->default_value("10"))
	("qmin", "Minimal loss probability/bias", cxxopts::value(Pl_min)->default_value("0"))
	("qmax", "Maximal loss probability/bias", cxxopts::value(Pl_max)->default_value("0"))
	
	("v", "verbosity switch", cxxopts::value(verbose)->default_value("0"))
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
	cout << "lmin:" << L_min << ";lmax:" << L_max << endl;
	cout << "pmin:" << P_min << ";pmax:" << P_max << ";nP:" << Np << endl;
	cout << "qmin:" << Pl_min << ";qmax:" << Pl_max << ";nPl:" << Npl << endl;
	cout << "n:" << n << "seed:" << seed << ";thread:" << thread << endl;
	cout << "noise model:" << N << ";loss model:" << Nl <<endl;
	
	if (fname == "") {
		fname = "l=" + to_string(L_min) + ",p=(" + to_string(P_min).substr(3,2) + "," + to_string(P_max).substr(3,2) + "),n=" +to_string(n) + To_string(s) + "," + To_string(N) + ".out";
	}
	
	if (test) {
		cluster test_cluster({L_min,L_min,L_min,0}, s);
		int i = 0;
		do {
			testDecoding(test_cluster, P_min, Pl_min, seed, s, N, Nl, verbose);
			i++;
		}
		while (i < times);
	} else{
		loopDecoding(L_min, L_max, n, P_min, P_max,  Np - 1, Pl_min, Pl_max, Npl - 1, fname, s, N, Nl, verbose, thread, use_env, make_corrections);
	}
}

///######## big 2D run ########
///./simulate -s PLANE --qmin 0 --qmax 0.25 --pmin 0 --pmax 0.008  --Np 20 --Npl 20 -n 1000 --lmin 3 -v 1 #big loss test EM2
///./simulate -s PLANE --qmin 0 --qmax 0.25 --pmin 0 --pmax 0.03  --Np 20 --Npl 20 -n 100 --lmin 3 -v 1 -N INDEP  #big loss test INDEP
///./simulate -s PLANE --pmin 0 --pmax 0.012 --qmin 0 --qmax 0 --Np 10 --Npl 1  -n 1000 --lmin 3 -v 1 -N GATE --Nl toLoss #pauli error only

///######## 1D run ########
//### 211 run: (*p_th = 8.4%, most recent test: 51162075, note the slow convergence*) 
//./simulate -s PLANE --qmin 0 --qmax 0 --pmin 0.05 --pmax 0.1  --Np 20 --Npl 1 -n 1000 --lmin 3 -v 1 -N INDEP_211
//### INDEP run: (*p_th = 3%*)
///./simulate -s PLANE --qmin 0 --qmax 0 --pmin 0.02 --pmax 0.05  --Np 20 --Npl 1 -n 1000 --lmin 3 --lmax 17 -v 1 -N INDEP
//### GATE/EM2/EM2_full run: (*p_th = 0.58%*)
//./simulate -s PLANE --qmin 0 --qmax 0 --pmin 0 --pmax 0.009  --Np 20 --Npl 1 -n 10000 --lmin 3 --lmax 17 -v 1 -N GATE
//./simulate -s PLANE --qmin 0 --qmax 0 --pmin 0 --pmax 0.009  --Np 20 --Npl 1 -n 10000 --lmin 3 --lmax 17 -v 1 -N EM2
//### GATE_biased run \beta = 1000 (*p_ref = 0.74% ArXiv 1308.4776 *)
//./simulate -s PLANE --qmin 1000 --pmin 0.005 --pmax 0.01  --Np 30 --Npl 1 -n 10000 --lmin 3 --lmax 21 -v 1 -N GATE_biased



//### LOSS run: (*p_th = 25%*)

//### Loss EM2 run:
///./simulate -s PLANE --qmin 0 --qmax 0 --pmin 0 --pmax 0.008  --Np 20 --Npl 1 -n 10000 --lmin 3 -v 1 -N GATE


///./simulate -s PLANE --qmin 0 --qmax 0.3 --pmin 0 --pmax 0  --Np 1 --Npl 20 -n 1000 --lmin 3 --lmax 7 -v 1 -N INDEP --Nl toLoss

///######## tests ########
//timing test
///./simulate -s PLANE --qmin 0.05 --qmax 0.05 --pmin 0.004 --pmax 0.008  --Np 10 --Npl 1 -n 500 --lmin 3 -v 1

//loss test
///./simulate -N INDEP -s PLANE --pmin 0 --qmin 0 --lmin 7 -v 2 --test
