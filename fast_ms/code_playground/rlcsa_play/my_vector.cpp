#include <vector>
#include "rlcsa/bits/bitvector.h"
#include "rlcsa/bits/rlevector.h"

using namespace std;
using namespace CSA;

//typedef pair<unsigned int, unsigned int> pair_type;
//

string p2s(pair_type p) {
    return "(" + to_string(p.first) + ", " + to_string(p.second) + ")";
}

size_t random_index(const size_t max_idx) {
    return static_cast<size_t> (max_idx * static_cast<unsigned long> (std::rand()) / (RAND_MAX + 1UL));
}


int main(int argc, char** argv){
    RLEVector::Encoder encoder(32);
    //           1
    // 01234567890123456789
    // 11110000101111000010
    size_t v = 0;
    for (int i = 0; i < 10; i++){
        size_t a = random_index(5);
        size_t b = random_index(6);
        encoder.addRun(v+a, b);
        v += (a + b);
    }
    encoder.flush();
    cout << v << endl;

    RLEVector rle(encoder, v);
    RLEVector::Iterator* it = new RLEVector::Iterator(rle);
    int nones = static_cast<int>(it->rank(rle.getSize()));
    int nruns = it->countRuns();
    int size = rle.getSize();

    cout << "array of length: " << size << endl;
    cout << "runs: " << nruns << endl;
    cout << "ones: " << nones << endl;
    cout << "contents" << endl;
    for(int i=0; i < size; i++)
        cout << static_cast<int>(it->isSet(i)) << " ";
    cout << endl;

    cout << "*************" << endl;
    cout << "sel(0): " << it->select(0) << endl;
    for(int i=0; i < nruns; i++){
        pair_type run_start = it->selectNextRun(size);
        cout << p2s(run_start) << ": ";
        for (int j = 0; j < run_start.second; j++)
            cout << ".";
        cout << endl;
    }
    cout << "*************" << endl;
    return 0;
}
