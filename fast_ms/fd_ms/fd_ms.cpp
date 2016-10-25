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


struct range{
    uint8_t w_start, w_len; // word coords
    uint8_t iw, jw;         // interval coords
};

struct interval{
    int lb = 0, ub = 0;
};

struct omega {
    int idx = 0, len = 0;
};

struct interval bstep(csa_wt<> sa, uint8_t lb, uint8_t ub, uint8_t c){
    struct interval i;
    
    i.lb = (int)(sa.C[sa.char2comp[c]] + sa.wavelet_tree.rank(lb, c));
    i.ub = (int)(sa.C[sa.char2comp[c]] + sa.wavelet_tree.rank(ub + 1, c) - 1);
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
    i.ub = (int) sa.size() - 1;
    uint8_t k = w.len;
    
    while(k-- > 0)
        i = bstep(sa, i.lb, i.ub, t[k]);
    return i;
}

struct interval bstep_along_rev(csa_wt<> sa, int_vector<8> t, struct omega w){
    struct interval i;
    i.ub = (int) (sa.size() - 1);
    uint8_t k = 0;
    
    while(k++ < w.len)
        i = bstep(sa, i.lb, i.ub, t[k]);
    return i;
}


void init_interval(struct interval *I, cst_sct3<> *st, const uint8_t c){
    uint8_t cc = st->csa.char2comp[c];
    
    I->lb = (int) st->csa.C[cc];
    I->ub = (int) (st->csa.C[cc + 1] - 1);
}

void phase1(uint *runs, cst_sct3<> *st_of_s_ptr, int_vector<8> *t_ptr){
    uint8_t k, c;
    struct interval I;
    struct omega w;
    cst_sct3<>::node_type v;
    cst_sct3<> st_of_s = (*st_of_s_ptr);
    int_vector<8> t = (*t_ptr);
    
    k = t.size();
    // start with the single char interval
    c = t[k - 1]; w.idx = k - 1; w.len = 1;
    // initialize to the interval of t[k - 1]
    init_interval(&I, &st_of_s, c);
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
            w.len = (int) st_of_s.depth(st_of_s.parent(v));
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
}

uint matching_stat(bit_vector ms, uint k){
    if(k == -1)
        return (uint) 1;
    return ((uint) bit_vector::select_1_type (&ms)(k + 1)) - (2 * k);
}

void output_ms(bit_vector ms, int ms_idx){
    cout << "ms: ";
    for(int i=0; i<ms.size(); i++)
        cout << ms[i];
    cout << endl << "    ";
    for(int i=0; i<ms.size(); i++)
        cout << (i == ms_idx ? "*" : " ");
    cout << endl;
}

void phase2(uint *runs, bit_vector *ms, cst_sct3<> *st_of_s_ptr, int_vector<8> *t_ptr){
    int k, c, h_star = 0, k_prim, ms_idx = 0;
    struct interval I;
    struct omega w;
    cst_sct3<>::node_type v;
    cst_sct3<> st_of_s = (*st_of_s_ptr);
    int_vector<8> t = (*t_ptr);
    


    k = 0, c = t[0];
    // initialize to the interval of t[k - 1]
    init_interval(&I, &st_of_s, c);
    w.idx = 0, w.len = 1;
    v = st_of_s.child(st_of_s.root(), c);

    while(k < t.size() - 1){
        output_ms(*ms, ms_idx);
        for(h_star = k + 1; I.lb <= I.ub; ){
            c = t[h_star];
            I = bstep(st_of_s.csa, I.lb, I.ub, c);
            if(I.lb <= I.ub){
                w.len++;
                v = st_of_s.wl(v, c);
                h_star ++;
            }
        }
        // I.lb > I.ub. MS[k] = h_star - k
        for(int i = 0; i < h_star - k - matching_stat(*ms, k - 1) + 1; i++)
            (*ms)[ms_idx++] = 0;
        if(h_star - k - matching_stat(*ms, k - 1) + 1 > 0)
            (*ms)[ms_idx++] = 1;

        k_prim = w.idx + w.len;
        do {
            v = st_of_s.parent(v);
            w.len = (int)st_of_s.depth(v);
            I = bstep_along_rev(st_of_s.csa, t, w);
            I = bstep(st_of_s.csa, I.lb, I.ub, t[h_star]);
            k_prim -= w.len;
        } while(I.lb > I.ub);
        

        for(int i=k+1; i<=k_prim - 1; i++)
            (*ms)[ms_idx++] = 1;
        
        //update omega
        w.idx = k_prim; w.len = h_star - k_prim + 1;
        k = k_prim;
    }
    output_ms(*ms, ms_idx);
    h_star = k + 1;
    for(int i = 0; i < h_star - k - matching_stat(*ms, k - 1) + 1; i++)
        (*ms)[ms_idx++] = 0;
    if(h_star - k - matching_stat(*ms, k - 1) + 1 > 0)
        (*ms)[ms_idx++] = 1;
    output_ms(*ms, ms_idx);
}


int main(int argc, char **argv){
    int_vector<8> t = {'a', 'b', 'b', 'a', 'b', 'a'};
    bit_vector ms(2 * t.size());
    /*
    ms[2] = 1; ms[3] = 1; ms[5] = 1; ms[8] = 1; ms[9] = 1; ms[11] = 1;
    for(int i=0; i<ms.size(); i++)
        cout << "ms["<< i << "] = " << ms[i] << endl;
    for(int i=-1; i<6; i++)
        cout << "MS["<< i << "] = " << (int) matching_stat(ms, i) << endl;
    */
    uint runs[t.size()];
    cst_sct3<> st_of_s, st_of_s_rev;
    construct_im(st_of_s, "aba", 1);
    construct_im(st_of_s_rev, "aba", 1);

    phase1(runs, &st_of_s, &t);
    phase2(runs, &ms, &st_of_s_rev, &t);
    for(int i=0; i<t.size(); i++){
        cout << "MS[" << i << "] = " << (int) matching_stat(ms, i) << endl;
    }

    //for(int i = 1; i < t.size(); i++){
    //    cout << "runs[" << i << "] = " << (int) runs[i] << endl;
    //}

//	cout << " i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << endl;
//	csXprintf(cout, "%2I %2S %3s %3P %2p %3B   %:3T", sa_of_s);
}

