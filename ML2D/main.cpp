//===================================================================//
//===================   ML Assisted decoding (make)   =====================//
//=====================        main.cpp     =========================//
//===================================================================//

//system libraries
#include <fstream>
#include <functional>
#include <algorithm>
#include <utility>
#include <iomanip>
#include <assert.h>
#include <time.h>

#include <filesystem>
#include <iostream>

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
namespace fs = filesystem;

static int window_size = 5;//this is defined twice. try to reduce

void testDecoding(cppflow::model model, const int L, const int M, const double p, const int seed, bool binary_output, bool make_corrections, bool decode_with_NN, subsurfacetype surf, noisemodel N, bool dir, float cutoff, bool verbose){
	subcluster testcluster({L,M,0},surf);
	if (N == MANUAL){
		testcluster.addError();
	} else {
		testcluster.addNoiseWithX(p, N, seed);
	}
	
	testcluster.getx_measurements();
	testcluster.getz_measurements();
	testcluster.printQubit();
	
	testcluster.decodeWithNN(model, binary_output, verbose, cutoff);
	testcluster.getx_measurements();
	testcluster.getz_measurements();
	testcluster.printQubit();
	
	if (surf == subPLANE && testcluster.decodeWithMWPMLoss(1,dir,make_corrections)[0] == 1) {
		cout << "success" << endl; //correction successful
	} else if (surf == subTORUS) {
		vector<int> parity = testcluster.decodeWithMWPMLoss(1,dir,make_corrections);
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
							testcluster.addNoiseWithX(p, N);
							testcluster.getx_measurements();
							testcluster.getz_measurements();

							testcluster.decodeWithNN(model, binary_output);
							testcluster.getx_measurements();
							testcluster.getz_measurements();
							
							if (surf == subPLANE && testcluster.decodeWithMWPMLoss(verbose, dir, 0)[0] == 1) {
								num_add++; //correction successful
							} else{
								vector<int> parity = testcluster.decodeWithMWPMLoss(verbose, dir, 1);
								if (parity[0]==1&&parity[1]==1) {
									num_add++;
								}
							}
						}
						num_correct += num_add;
					}, use_env);
				} else {
					for (int i = 0; i < trials; i++) {
						testcluster.addNoiseWithX(p, N);
						testcluster.getx_measurements();
						testcluster.getz_measurements();

						testcluster.decodeWithNN(model, binary_output);
						testcluster.getx_measurements();
						testcluster.getz_measurements();
						
						if (surf == subPLANE && testcluster.decodeWithMWPMLoss(verbose, dir)[0] == 1) {
							num_correct ++; //correction successful
						} else{
								vector<int> parity = testcluster.decodeWithMWPMLoss(verbose, dir, 1);
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
							testcluster.addNoiseWithX(p, N);
							testcluster.getx_measurements();
							testcluster.getz_measurements();
							
							if (surf == subPLANE && testcluster.decodeWithMWPMLoss(verbose,dir, 0)[0] == 1) {
								num_add++; //correction successful
							} else{
								vector<int> parity = testcluster.decodeWithMWPMLoss(verbose, dir, 1);
								if (parity[0]==1&&parity[1]==1) {
									num_add++;
								}
							}
						}
						num_correct += num_add;
					}, use_env);
				} else {
					for (int i = 0; i < trials; i++) {
						testcluster.addNoiseWithX(p, N);
						testcluster.getx_measurements();
						testcluster.getz_measurements();
						
						if (surf == subPLANE && testcluster.decodeWithMWPMLoss(verbose, dir)[0] == 1) {
							num_correct ++; //correction successful
						} else{
								vector<int> parity = testcluster.decodeWithMWPMLoss(verbose, dir, 1);
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
				testcluster.addNoiseWithX(p, N, seed + 1);
			} else {
				testcluster.addNoiseWithX(p, N, 0);
			}
			testcluster.getx_measurements();
			testcluster.getz_measurements();
			//generate method 1
			//randomly pick an erronous stabilizer measurement, random pick one out of 4 neighbors.
			vector<int> test_pos;
			// for (int c = 0; c < 2*testcluster.S.x*testcluster.S.y; c++) {
			// 	if(testcluster.stabs[c] < 0){
			// 		test_pos.push_back(c);
			// 	}
			// }

			//generate method 2
			//only pick x measurements
			for (int c = 0; c < testcluster.S.x*testcluster.S.y; c++) {
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
				if (h < testcluster.S.x*testcluster.S.y){
					subvertex aVertex(h,testcluster.S);
					c = aVertex.physical_qubits[dist2(engine)];
				} else {
					testcluster.Subvertex_z aVertex(h-testcluster.S.x*testcluster.S.y);
					c = aVertex.physical_qubits[dist2(engine)];
				}
				
				testcluster.Subcoord C(c);
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
	("s, surf_type", "Surface type", cxxopts::value(s)->default_value("subTORUS"))
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

//./simulate -s subTORUS --pmin 0 --pmax 0.18  --Np 25 -n 1000 --Lmin 5 --Lmax 5 -v 1 -d ~/ML/ -m "model,L=5(7),layer=3x128,epochs=10000,p=" --decode_with_NN
//./simulate -s subTORUS --pmin 0 --pmax 0.12  --Np 25 -n 1000 --Lmin 3 --Lmax 20 -v 1 -d ~/ML/  --fname test.out

//./simulate -s subTORUS --pmin 0.036 --Np 10 --Lmin 10 -v 1 --test --make_corrections -d /scratch/users/ladmon/ML/ -m "model_h,L=5(7),layer=3x128,epochs=100000,p=0.036" --binary
//./simulate -s subTORUS --pmin 0.02 --pmax 0.02  --Np 20 -n 1 --Lmin 7 -v 1 --generate -d ~/ML

//model,L=5(7),layer=5x512,epochs=1000,p=0.1068