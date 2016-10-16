/*
 * fabio_djamal_ms.cpp
 *
 *  Created on: Oct 13, 2016
 *      Author: denas
 */

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/csa_wt.hpp>

using namespace std;
using namespace sdsl;

struct interval{
    uint8_t lb = 0, ub = 0;
};

struct interval bstep(csa_wt<> index, uint8_t lb, uint8_t ub, uint8_t c){
    struct interval i;
    
    i.lb = index.C[index.char2comp[c]] + index.wavelet_tree.rank(lb, c);
    i.ub = index.C[index.char2comp[c]] + index.wavelet_tree.rank(ub + 1, c) - 1;
    return i;
}


int main(int argc, char **argv){
    //char t[] = {'a', 'b', 'b', 'a', 'b', 'a'};
    int_vector<8> t = {'a', 'b', 'b', 'a', 'b', 'a'};
    
    
    //uint8_t w_idx, w_len;
    //uint8_t k = t.size() - 1;
    //uint8_t c;

	csa_wt<> sa_of_s;
	construct_im(sa_of_s, "abbaba", 1);

    
    int_vector<8> runs; // the runs vector

    for(int k = 0; k < sa_of_s.sigma; k++){
        cout << k << sa_of_s.comp2char[k] << endl;
        cout << "C[" << k << "] = " << sa_of_s.C[k] << endl;
        cout << endl;
    }
    
    struct interval I = bstep(sa_of_s, 0, 5, 'b');
    cout << "{" << I.lb << ", " << I.ub << "}" << endl;
    I = bstep(sa_of_s, I.lb, I.ub, 'b');
    cout << "{" << I.lb << ", " << I.ub << "}" << endl;
    

//	cout << " i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << endl;
//	csXprintf(cout, "%2I %2S %3s %3P %2p %3B   %:3T", sa_of_s);
}

