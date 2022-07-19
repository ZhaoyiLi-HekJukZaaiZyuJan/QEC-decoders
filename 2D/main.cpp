//########################## Surface Code ##########################
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

void testDecoding(const int L, const int M, const double p, subsurfacetype surf=subTORUS){
	subcluster testcluster({L,L,L},surf);
	testcluster.addNoise(p);

	cout << "qubit: " << endl;
	testcluster.printQubit();
	
	cout << "here" <<endl;
	cout << testcluster.decodeWithMWPM(1) << endl;

}

void innerSweepTestDecoding(const int& level, const int& l, const int& i, int& success, int& trials, subcluster testcluster){
	if (level == 0){
		testcluster.getx_vec();
		trials ++;
		int ifsuccess = testcluster.decodeWithMWPM(0);
		success += ifsuccess;
		if (ifsuccess != 0){
			// testcluster.printQubit();
			// cout << endl;
		}
		return;
	} 
	for(int j = i + 1; j < testcluster.S.x*testcluster.S.y*2; j++) {
		testcluster.addNoiseManually(j);
		innerSweepTestDecoding(level - 1, l, j, success, trials, testcluster);
		testcluster.addNoiseManually(j);
	}
}

void SweepTestDecoding(const int L, const int M, const int& l, subsurfacetype surf=subTORUS){
	subcluster testcluster({L,L,L},surf);
	int success = 0;
	int trials = 0;
	innerSweepTestDecoding(l, l, -1, success, trials, testcluster);
	cout << "success:" << success << endl;
	cout << "fail:" << trials-success << endl;
	cout << "trials:" << trials << endl;
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
				cout << L << "," << p << ",0," << num_correct << "\n";
			}
		}
	}
}

int main(int argc, const char *argv[]) {
	string fname, code_model, str;
	int L_min, L_max, level, verbose, n, Np, seed;
	bool test, sweep, thread, use_env;
	double P_max, P_min;
	subsurfacetype s;

	cxxopts::Options options(*argv,"2D Simulator");

	options.add_options()
	("f, fname", "Output filename", cxxopts::value(fname))
	("s, surf_type", "Surface type", cxxopts::value(s)->default_value("PLANE"))
	("Lmin", "Minimal size of mesh", cxxopts::value(L_min)->default_value("3"))
	("Lmax", "Maximal size of mesh", cxxopts::value(L_max)->default_value("17"))
	("l", "level", cxxopts::value(level)->default_value("0"))
	("n", "Number of trials", cxxopts::value(n)->default_value("10000"))
	
	("v", "verbosity switch", cxxopts::value(verbose)->default_value("0"))

	("Np", "z error p Points", cxxopts::value(Np)->default_value("10"))
	("pmin", "Minimal z error probability", cxxopts::value(P_min)->default_value("0.001"))
	("pmax", "Maximal z error probability", cxxopts::value(P_max)->default_value("0.008"))
	
	("c", "code model", cxxopts::value(code_model)->default_value("2D"))
	("test", "test switch", cxxopts::value(test)->default_value("0"))
	("sweep", "sweep switch", cxxopts::value(sweep)->default_value("0"))
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
	} else if (sweep){
		for (int level = 1; level <= L_min*L_min*2; level ++ ){
			cout << level << endl;
			SweepTestDecoding(L_min, L_min, level, s);
		}
		SweepTestDecoding(L_min, L_min, level, s);
	} else {
		loopDecoding(L_min, L_max, n, P_min, P_max, Np - 1, fname, s ,verbose, thread, use_env);	
	}

}


//./simulate -s PLANE --pmin 0 --pmax 0.1 --Np 20 -n 10000 --Lmin 3 -v 1 -v 1
