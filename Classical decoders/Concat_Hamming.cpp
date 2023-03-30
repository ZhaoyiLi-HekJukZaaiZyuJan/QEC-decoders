#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace std;

const int K = 4;  // length of the original message bits
const int R1 = 1; // number of parity bits for Hamming code 1
const int R2 = 2; // number of parity bits for Hamming code 2
const int R = R1 + R2; // total number of parity bits for the concatenated code
const int N = K + R; // total length of the encoded message

// generate random binary message of length K
void generate_message(bool message[K]) {
    for (int i = 0; i < K; i++) {
        message[i] = rand() % 2;
    }
}

// calculate Hamming code 1 parity bits
void hamming_code_1(bool message[K], bool parity[R1]) {
    parity[0] = message[0] ^ message[1] ^ message[3];
}

// calculate Hamming code 2 parity bits
void hamming_code_2(bool message[K], bool parity[R2]) {
    parity[0] = message[0] ^ message[1] ^ message[3];
    parity[1] = message[0] ^ message[2] ^ message[3];
}

// calculate concatenated Hamming code parity bits
void concatenated_hamming_code(bool message[K], bool parity[R]) {
    bool parity1[R1];
    hamming_code_1(message, parity1);
    bool encoded[K + R1];
    for (int i = 0; i < K; i++) {
        encoded[i] = message[i];
    }
    for (int i = 0; i < R1; i++) {
        encoded[K + i] = parity1[i];
    }
    bool parity2[R2];
    hamming_code_2(encoded, parity2);
    for (int i = 0; i < R2; i++) {
        parity[R1 + i] = parity2[i];
    }
}

// simulate transmission with probability of error p
bool simulate_transmission(bool encoded[N], float p) {
    for (int i = 0; i < N; i++) {
        if ((float)rand() / RAND_MAX < p) {
            encoded[i] = !encoded[i];
        }
    }
    bool received[N];
    for (int i = 0; i < N; i++) {
        received[i] = encoded[i];
    }
    return received;
}

// check for errors in received message
bool check_for_errors(bool received[N], bool parity[R]) {
    bool error = false;
    bool parity_check[R];
    for (int i = 0; i < R1; i++) {
        parity_check[i] = received[K + i];
    }
    for (int i = 0; i < R2; i++) {
        parity_check[R1 + i] = received[K + R1 + i];
    }
    for (int i = 0; i < R; i++) {
        if (parity_check[i] != parity[i]) {
            error = true;
            break;
        }
    }
    return error;
}

int main() {
    srand(time(NULL)); // seed random number generator with current time
    int num_trials = 100000; // number of Monte Carlo trials
    int num_errors = 0; // number of errors detected
    float p = 0.1; // probability of error
   
