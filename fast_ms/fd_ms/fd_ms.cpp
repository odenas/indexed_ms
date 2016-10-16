/*
 * fabio_djamal_ms.cpp
 *
 *  Created on: Oct 13, 2016
 *      Author: denas
 */

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_trees.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/bp_support.hpp>
#include <sdsl/csa_wt.hpp>

using namespace std;
using namespace sdsl;

struct interval{
    uint8_t lb = 0, ub = 0;
};

struct omega {
    uint8_t idx = 0, len = 0;
};

struct interval bstep(csa_wt<> sa, uint8_t lb, uint8_t ub, uint8_t c){
    struct interval i;
    
    i.lb = sa.C[sa.char2comp[c]] + sa.wavelet_tree.rank(lb, c);
    i.ub = sa.C[sa.char2comp[c]] + sa.wavelet_tree.rank(ub + 1, c) - 1;
    return i;
}

void output_node(const typename cst_sct3<>::node_type &v, const cst_sct3<> &cst){
    cout << cst.depth(v) << "-[" << cst.lb(v) << "," << cst.rb(v) << "]" << endl;
}
void output_omega(const struct omega w, int_vector<8> t){
    cout << "w = ";
    for(int i=0; i<w.len; i++)
        cout << t[w.idx + i];
    cout << endl;
}
void output_interval(const struct interval I){
    cout << "{";
    if(I.lb > I.ub)
        cout << "-, -";
    else
        cout << (int) I.lb << ", " << (int) I.ub;
    cout << "}" << endl;
}

struct interval bstep_along(csa_wt<> sa, int_vector<8> t, struct omega w){
    struct interval i;
    i.ub = sa.size() - 1;
    uint8_t k = w.len;
    
    while(k-- > 0)
        i = bstep(sa, i.lb, i.ub, t[k]);
    return i;
}


int main(int argc, char **argv){
    int_vector<8> t = {'a', 'b', 'b', 'a', 'b', 'a'};
    int_vector<8> runs = {0, 0, 0, 0, 0, 0}; // the runs vector

    uint8_t k, c;
    struct interval I;
    struct omega w, a;
    cst_sct3<>::node_type v;

    cst_sct3<> st_of_s;
    construct_im(st_of_s, "ab", 1);

    k = t.size();
    // start with the single char interval
    c = t[k - 1]; w.idx = k - 1; w.len = 1;
    // initialize to the interval of t[k - 1]
    I.lb = st_of_s.csa.C[st_of_s.csa.char2comp[c]];
    I.ub = st_of_s.csa.C[st_of_s.csa.char2comp[c] + 1] - 1;
    v = st_of_s.child(st_of_s.root(), c);

    while(--k > 0){
        c = t[k-1];
        cout << "[" << (int) k << "]" << endl;
        output_omega(w, t);
        cout << "c = " << c << endl;

        I = bstep(st_of_s.csa, I.lb, I.ub, t[k - 1]);
        output_interval(I);

        if(I.lb > I.ub){
            runs[k] = 0;
            // update I to the parent of the proper locus of w
            cout << "|parent(v)| = " << (int) st_of_s.depth(st_of_s.parent(v)) << endl;
            w.len = st_of_s.depth(st_of_s.parent(v));
            I = bstep_along(st_of_s.csa, t, w);
            I = bstep(st_of_s.csa, I.lb, I.ub, t[k - 1]);
            cout << "-" << endl;
        } else {
            v = st_of_s.wl(v, c); // update v
            runs[k] = 1;
            cout << "+" << endl;
        }
        w.idx--; w.len++;
        cout << endl;
    }
    
    for(int i = 1; i < runs.size(); i++){
        cout << "runs[" << i << "] = " << (int) runs[i] << endl;
    }

//	cout << " i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << endl;
//	csXprintf(cout, "%2I %2S %3s %3P %2p %3B   %:3T", sa_of_s);
}

