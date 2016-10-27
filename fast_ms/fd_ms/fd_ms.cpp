/*
 * fabio_djamal_ms.cpp
 *
 *  Created on: Oct 13, 2016
 *      Author: denas
 */

#include <iostream>
#include <string>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_trees.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/bp_support.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/bit_vectors.hpp>


using namespace std;
using namespace sdsl;


typedef struct interval{
    int lb = 0, ub = 0;
} Interval;

typedef struct omega {
    int idx = 0, len = 0;
} Omega;

typedef struct w_state {
    Interval I;
    Omega w;
    cst_sct3<>::node_type v; // stree node
} Wstate;


Interval bstep(csa_wt<> sa, uint8_t lb, uint8_t ub, uint8_t c){
    Interval i;
    
    i.lb = (int)(sa.C[sa.char2comp[c]] + sa.wavelet_tree.rank(lb, c));
    i.ub = (int)(sa.C[sa.char2comp[c]] + sa.wavelet_tree.rank(ub + 1, c) - 1);
    return i;
}

void output_node(const typename cst_sct3<>::node_type &v, const cst_sct3<> &cst){
    cout << cst.depth(v) << "-[" << cst.lb(v) << "," << cst.rb(v) << "]" << endl;
}
void output_omega(const Omega w, int_vector<8> t){
    cout << "w = ";
    for(int i=0; i<w.len; i++)
        cout << t[w.idx + i];
    cout << endl;
}
void output_interval(const Interval I){
    cout << "{";
    if(I.lb > I.ub)
        cout << "-, -";
    else
        cout << (int) I.lb << ", " << (int) I.ub;
    cout << "}" << endl;
}


void output_partial_vec(const bit_vector v, const int idx, const char name[], bool verbose){
    if (!verbose)
        return ;
    cout << name << ": ";
    for(int i=0; i<v.size(); i++)
        cout << v[i];
    cout << endl << "    ";
    for(int i=0; i<v.size(); i++)
        cout << (i == idx ? "*" : " ");
    cout << endl;
}

Interval bstep_along(csa_wt<> sa, string t, Omega w){
    Interval i;
    i.ub = (int) sa.size() - 1;
    uint8_t k = w.len;
    
    while(k-- > 0)
        i = bstep(sa, i.lb, i.ub, t[w.idx + k]);
    return i;
}

Interval bstep_along_rev(csa_wt<> sa, string t, Omega w){
    Interval i;
    i.ub = (int) (sa.size() - 1);
    //uint8_t k = 0;
    
    for(uint8_t k = w.idx; k < w.idx + w.len; k++)
        i = bstep(sa, i.lb, i.ub, t[k]);
    return i;
}


uint MS(bit_vector ms, uint k){
    if(k == -1)
        return (uint) 1;
    return ((uint) bit_vector::select_1_type (&ms)(k + 1)) - (2 * k);
}

Wstate init_state(int idx, cst_sct3<> *stree, char c){
    Wstate state;

    state.w.idx = idx;
    state.w.len = 1;
    
    state.v = stree->child(stree->root(), c);
    
    uint8_t cc = stree->csa.char2comp[c];
    
    state.I.lb = (int) stree->csa.C[cc];
    state.I.ub = (int) (stree->csa.C[cc + 1] - 1);

    return state;
}

void phase1(bit_vector runs, cst_sct3<> *st_of_s_ptr, string t, bool verbose){
    cst_sct3<> st_of_s = (*st_of_s_ptr);
    uint8_t k = t.size(), c = t[k - 1];
    Wstate state = init_state(k - 1, st_of_s_ptr, c);
    
    while(--k > 0){
        c = t[k-1];

        state.I = bstep(st_of_s.csa, state.I.lb, state.I.ub, t[k - 1]);
        if(state.I.lb > state.I.ub){
            runs[k] = 0;
            // update I to the parent of the proper locus of w until we can extend by 'c'
            do{
                state.v = st_of_s.parent(state.v);
                state.w.len = (int)st_of_s.depth(state.v);
                state.I = bstep_along(st_of_s.csa, t, state.w);
                state.I = bstep(st_of_s.csa, state.I.lb, state.I.ub, c);
            } while(state.I.lb > state.I.ub);
        } else {
            state.v = st_of_s.wl(state.v, c); // update v
            runs[k] = 1;
        }
        state.w.idx--; state.w.len++;
        output_partial_vec(runs, k, "runs", verbose);
    }
}


void phase2(bit_vector runs, bit_vector *ms, cst_sct3<> *st_of_s_ptr, string t, bool verbose){
    cst_sct3<> st_of_s = (*st_of_s_ptr);
    int k = 0, c = t[k], h_star = 0, k_prim = 0, ms_idx = 0;
    Wstate state = init_state(k, st_of_s_ptr, c);
    
    while(k < t.size()){
        output_partial_vec(*ms, ms_idx, "ms", verbose);
        for(h_star = k + state.w.len; state.I.lb <= state.I.ub &&  h_star < t.size(); ){
            c = t[h_star];
            state.I = bstep(st_of_s.csa, state.I.lb, state.I.ub, c);
            if(state.I.lb <= state.I.ub){
                state.w.len++;
                state.v = st_of_s.wl(state.v, c);
                h_star ++;
            }
        }
        // I.lb > I.ub. MS[k] = h_star - k
        for(int i = 0; i < h_star - k - MS(*ms, k - 1) + 1; i++)
            (*ms)[ms_idx++] = 0;
        if(h_star - k - MS(*ms, k - 1) + 1 > 0)
            (*ms)[ms_idx++] = 1;

        k_prim = state.w.idx + state.w.len;
        if(h_star < t.size()){
            do {
                state.v = st_of_s.parent(state.v);
                state.w.idx += (state.w.len - st_of_s.depth(state.v));
                state.w.len = (int)st_of_s.depth(state.v);

                state.I = bstep_along_rev(st_of_s.csa, t, state.w);
                state.I = bstep(st_of_s.csa, state.I.lb, state.I.ub, t[h_star]);
            } while(state.I.lb > state.I.ub);
            k_prim -= state.w.len;
        }

        for(int i=k+1; i<=k_prim - 1; i++)
            (*ms)[ms_idx++] = 1;
        
        //update omega
        state.w.idx = k_prim; state.w.len = h_star - k_prim + 1;
        // update v
        state.v = st_of_s.wl(state.v, c);
        k = k_prim;
    }
    output_partial_vec(*ms, ms_idx, "ms", verbose);
}

bit_vector compute_ms(string T, string S, bool verbose){
    string Srev {S};
    for(int i=0; i<S.length(); i++)
        Srev[S.length() - i - 1] = S[i];

    bit_vector ms(2 * T.length());
    bit_vector runs(T.length());
    cst_sct3<> st_of_s, st_of_s_rev;
    construct_im(st_of_s, S, 1);
    construct_im(st_of_s_rev, Srev, 1);

    if(verbose){
        cout << " i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << endl;
        csXprintf(cout, "%2I %2S %3s %3P %2p %3B   %:3T", st_of_s.csa);
        cout << " i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << endl;
        csXprintf(cout, "%2I %2S %3s %3P %2p %3B   %:3T", st_of_s_rev.csa);
    }

    phase1(runs, &st_of_s, T, verbose);
    if(verbose)
        cout << endl;
    phase2(runs, &ms, &st_of_s_rev, T, verbose);
    if(verbose){
        for(int i=0; i<T.length(); i++)
            cout << "MS[" << i << "] = " << (int) MS(ms, i) << endl;
    }
    return ms;
}


int main(int argc, char **argv){
    bit_vector ms;

    if(argc == 2){ // process file
        string T, S, Srev ;
        //std::ifstream in_file {"/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/a.txt"};
        std::ifstream in_file {argv[1]};
        if(!in_file){
            cout << "could not open file " << argv[1] << endl;
            return 1;
        }

        while (in_file >> T >> S){
            cout << T + " " + S + " ";

            ms = compute_ms(T, S, 0);
            for (int i = 0; i < T.length(); i++)
                cout << MS(ms, i);
            cout << endl;
        }
    } else if (argc == 3) { // process T and S
        string T {argv[1]};
        string S {argv[2]};
        cout << T + " " + S + " ";

        ms = compute_ms(T, S, 0);
        for (int i = 0; i < T.length(); i++)
            cout << MS(ms, i);
        cout << endl;
    } else {
        string T {"bbab"};
        string S {"acabb"};
        ms = compute_ms(T, S, 1);
        cout << T + " " + S + " ";
        for (int i = 0; i < T.length(); i++)
            cout << MS(ms, i);
        cout << endl;

        //cout << "usage: " << argv[0] << "<file of strings>" << endl;
        //cout << "or" << endl;
        //cout << "usage: " << argv[0] << "<T> <S>" << endl;
        //return 1;
    }
    return 0;
}

