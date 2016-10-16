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
    int_vector<8> t = {'a', 'b', 'b', 'a', 'b', 'a'};
    int_vector<8> runs = {0, 0, 0, 0, 0, 0}; // the runs vector

    uint8_t w_idx, w_len;
    uint8_t k, c;
    struct interval I;

	csa_wt<> sa_of_s;
	construct_im(sa_of_s, "ab", 1);

    k = t.size();
    w_idx = k - 1;
    w_len = 1;
    I.lb = I.ub = 1;
    while(--k > 0){
        c = t[k-1];
        cout << "k = " << (int) k << endl;
        cout << "w = (" << (int) w_idx << ", " << (int) w_len << ")" << endl;
        cout << "c = " << c << endl;
        
        I = bstep(sa_of_s, I.lb, I.ub, t[k - 1]);
        cout << "{" << (int) I.lb << ", " << (int) I.ub << "}" << endl;

        if(I.lb > I.ub){
            runs[k] = 0;
            I = bstep(sa_of_s, 0, sa_of_s.size() - 1, c);
            w_len = 1;
            cout << "-" << endl;
        } else {
            w_len += I.ub - I.lb + 1;
            runs[k] = 1;
            cout << "+" << endl;
        }
        w_idx--;
        cout << endl;
    }
    
    //struct interval I = bstep(sa_of_s, 0, 5, 'b');
    //cout << "{" << (int) I.lb << ", " << (int) I.ub << "}" << endl;
    //I = bstep(sa_of_s, I.lb, I.ub, 'b');
    //cout << "{" << (int) I.lb << ", " << (int) I.ub << "}" << endl;


//	cout << " i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << endl;
//	csXprintf(cout, "%2I %2S %3s %3P %2p %3B   %:3T", sa_of_s);
}

