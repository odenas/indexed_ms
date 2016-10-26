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


void output_partial_vec(const bit_vector v, const int idx, const char name[]){
    cout << name << ": ";
    for(int i=0; i<v.size(); i++)
        cout << v[i];
    cout << endl << "    ";
    for(int i=0; i<v.size(); i++)
        cout << (i == idx ? "*" : " ");
    cout << endl;
    
}

Interval bstep_along(csa_wt<> sa, int_vector<8> t, Omega w){
    Interval i;
    i.ub = (int) sa.size() - 1;
    uint8_t k = w.len;
    
    while(k-- > 0)
        i = bstep(sa, i.lb, i.ub, t[k]);
    return i;
}

Interval bstep_along_rev(csa_wt<> sa, int_vector<8> t, Omega w){
    Interval i;
    i.ub = (int) (sa.size() - 1);
    //uint8_t k = 0;
    
    for(uint8_t k = 0; k < w.len; k++)
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

void phase1(bit_vector runs, cst_sct3<> *st_of_s_ptr, int_vector<8> *t_ptr){
    cst_sct3<> st_of_s = (*st_of_s_ptr);
    int_vector<8> t = (*t_ptr);
    uint8_t k = t.size(), c = t[k - 1];
    Wstate state = init_state(k - 1, st_of_s_ptr, c);
    
    while(--k > 0){
        c = t[k-1];

        state.I = bstep(st_of_s.csa, state.I.lb, state.I.ub, t[k - 1]);
        if(state.I.lb > state.I.ub){
            runs[k] = 0;
            // update I to the parent of the proper locus of w
            state.w.len = (int) st_of_s.depth(st_of_s.parent(state.v));
            state.I = bstep_along(st_of_s.csa, t, state.w);
            state.I = bstep(st_of_s.csa, state.I.lb, state.I.ub, t[k - 1]);
        } else {
            state.v = st_of_s.wl(state.v, c); // update v
            runs[k] = 1;
        }
        state.w.idx--; state.w.len++;
        output_partial_vec(runs, k, "runs");
    }
}


void phase2(bit_vector runs, bit_vector *ms, cst_sct3<> *st_of_s_ptr, int_vector<8> *t_ptr){
    cst_sct3<> st_of_s = (*st_of_s_ptr);
    int_vector<8> t = (*t_ptr);
    int k = 0, c = t[k], h_star = 0, k_prim = 0, ms_idx = 0;
    Wstate state = init_state(k, st_of_s_ptr, c);
    
    while(k < t.size()){
        output_partial_vec(*ms, ms_idx, "ms");
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
                k_prim -= state.w.len;
            } while(state.I.lb > state.I.ub);
        }

        for(int i=k+1; i<=k_prim - 1; i++)
            (*ms)[ms_idx++] = 1;
        
        //update omega
        state.w.idx = k_prim; state.w.len = h_star - k_prim + 1;
        // update v
        state.v = st_of_s.wl(state.v, c);
        k = k_prim;
    }
    output_partial_vec(*ms, ms_idx, "ms");
}


int main(int argc, char **argv){
    int_vector<8> t = {'a', 'b', 'b', 'a', 'b', 'a'};
    bit_vector ms(2 * t.size());
    bit_vector runs(t.size());
    cst_sct3<> st_of_s, st_of_s_rev;
    construct_im(st_of_s, "ab", 1);
    construct_im(st_of_s_rev, "ba", 1);

    cout << " i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << endl;
    csXprintf(cout, "%2I %2S %3s %3P %2p %3B   %:3T", st_of_s.csa);
    cout << " i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << endl;
    csXprintf(cout, "%2I %2S %3s %3P %2p %3B   %:3T", st_of_s_rev.csa);

    phase1(runs, &st_of_s, &t);
    cout << endl;
    phase2(runs, &ms, &st_of_s_rev, &t);

    for(int i=0; i<t.size(); i++){
        cout << "MS[" << i << "] = " << (int) MS(ms, i) << endl;
    }
}

