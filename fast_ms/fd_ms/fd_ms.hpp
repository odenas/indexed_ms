//
//  fd_ms.hpp
//  fast_ms
//
//  Created by denas on 11/2/16.
//  Copyright © 2016 denas. All rights reserved.
//

#ifndef fd_ms_h
#define fd_ms_h

#include <iostream>
#include <string>
#include <sdsl/suffix_trees.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "Bwt.hpp"
#include "stree.hpp"


using namespace std;
using namespace sdsl;

typedef unsigned long size_type;


namespace fdms{


class Interval{
private:
    Bwt *bwt;

public:
    size_type lb, ub;

    Interval(Bwt *bwt_, char c){
        bwt = bwt_;
        lb = bwt->C[bwt->char2int[c]];
        ub = bwt->C[bwt->char2int[c] + 1] - 1;
    }

    void set(size_type l, size_type u){
        lb = l;
        ub = u;
    }

    void bstep(char c){
        int cc = bwt->char2int[c];
        lb = bwt->C[cc] + bwt->rank(lb, c);
        ub = bwt->C[cc] + bwt->rank(ub + 1, c) - 1;
    }

    bool is_empty(){
        return lb > ub;
    }

    void dump(){
        std::cout << "{";
        if(is_empty())
            std::cout << "-, -";
        else
            std::cout << lb << ", " << ub;
        std::cout << "}" << endl;
    }
};

class Interval_sdsl{
private:
    cst_sada<> stree;

public:
    size_type lb, ub;

    Interval_sdsl(cst_sada<> *str, char c){
        stree = (*str);
        lb = stree.csa.C[stree.csa.char2comp[c]];
        ub = stree.csa.C[stree.csa.char2comp[c] + 1] - 1;
    }

    void set(size_type l, size_type u){
        lb = l;
        ub = u;
    }

    void bstep(char c){
        int cc = stree.csa.char2comp[c];
        lb = stree.csa.C[cc] + stree.csa.bwt.rank(lb, c);
        ub = stree.csa.C[cc] + stree.csa.bwt.rank(ub + 1, c) - 1;
    }

    bool is_empty(){
        return lb > ub;
    }
};

class Mstat{
private:
    bit_vector ms, runs;
    size_type ms_size;

    void output_partial_vec(const bit_vector v, size_type idx, const char name[], bool verbose){
        if (!verbose)
            return ;
        cout << name << ": ";
        for(size_type i = 0; i < (size_type)v.size(); i++)
            cout << v[i];
        cout << endl;
        for(size_type i = 0; i < (size_type)strlen(name) + 2; i++)
            cout << " ";
        for(size_type i = 0; i < (size_type)v.size(); i++)
            cout << (i == idx ? "*" : " ");
        cout << endl;
    }

    size_type find_k_prim(bit_vector runs, size_type k){
        size_t total_zeros = rank_support_v<0>(&runs)(runs.size());
        size_t zeros = rank_support_v<0>(&runs)(k + 1);

        if(total_zeros == zeros)
            return runs.size();

        bit_vector::select_0_type b_sel(&runs);
        return b_sel(zeros + 1);
    }

    bit_vector build_runs(string t, string s, Bwt bwt, const bool verbose){
        cst_sada<> st_of_s;
        construct_im(st_of_s, s, 1);

        bit_vector runs(ms_size);
        size_type k = ms_size, c = t[k - 1];
        //Interval I{&bwt, static_cast<char>(c)};
        Interval_sdsl I{&st_of_s, static_cast<char>(c)};

        cst_sada<>::node_type v = st_of_s.child(st_of_s.root(), c); // stree node
        while(--k > 0){
            c = t[k-1];
            I.bstep(c);
            if(I.is_empty()){
                runs[k] = 0;
                // update I to the parent of the proper locus of w until we can extend by 'c'
                do{
                    v = st_of_s.parent(v);
                    I.set(st_of_s.lb(v), st_of_s.rb(v));
                    I.bstep(c);
                } while(I.is_empty());
            } else {
                runs[k] = 1;
            }
            v = st_of_s.wl(v, c); // update v
            output_partial_vec(runs, k, "runs", verbose);
        }
        return runs;
    }

    size_type get_ms(bit_vector ms_, size_type k){
        if(k == -1)
            return (size_type) 1;
        return bit_vector::select_1_type (&ms_)(k + 1) - (2 * k);
    }

    bit_vector build_ms(string t, string s_rev, Bwt bwt, const bool verbose){
        bit_vector ms(ms_size * 2);
        cst_sada<> st_of_s;
        construct_im(st_of_s, s_rev, 1);
        size_type k = 0, h_star = k + 1, k_prim, ms_idx = 0;
        uint8_t c = t[k];
        //Interval I{&bwt, static_cast<char>(c)};
        Interval_sdsl I{&st_of_s, static_cast<char>(c)};

        cst_sada<>::node_type v = st_of_s.child(st_of_s.root(), c); // stree node

        while(k < ms_size){
            output_partial_vec(ms, ms_idx, "ms", verbose);
            for(; !I.is_empty() && h_star < ms_size; ){
                c = t[h_star];
                I.bstep(c);
                if(!I.is_empty()){
                    v = st_of_s.wl(v, c);
                    h_star ++;
                }
            }
            for(int i = 0; i < h_star - k - get_ms(ms, k - 1) + 1; i++)
                ms[ms_idx++] = 0;
            if(h_star - k - get_ms(ms, k - 1) + 1 > 0)
                ms[ms_idx++] = 1;

            if(h_star < ms_size){
                do {
                    v = st_of_s.parent(v);
                    I.set(st_of_s.lb(v), st_of_s.rb(v));
                    I.bstep(t[h_star]);
                } while(I.is_empty());
                h_star = h_star + 1;
            }
            // k_prim: index of the first zero to the right of k in runs
            k_prim = find_k_prim(runs, k);

            for(size_type i = k + 1; i <= k_prim - 1; i++)
                ms[ms_idx++] = 1;

            // update v
            v = st_of_s.wl(v, c);
            k = k_prim;
        }
        output_partial_vec(ms, ms_idx, "ms", verbose);
        return ms;
    }

public:
    Mstat(string& T,
          string& Sfwd, Bwt& Bwtfwd, bp_support_sada<>& bfwdbps,
          string& Srev, Bwt& Bwtrev, bp_support_sada<>& brevbps,
          const bool verbose){

        ms_size = T.size();

        Stree st(bfwdbps, Bwtfwd);
        //{"11011010110100011101001000"};
        //  01 34 6 89 1   567 9  2
        //            1         2
        size_type idx[] {0, 1, 3, 4, 6, 8, 9, 11, 15, 16, 17, 19, 22};
        for(size_type i = 0; i < 12; i++)
            cout << "[" << idx[i] << "] = " << st.depth(idx[i]) << endl;

        for(size_type i = 0; i < 12; i++){
            cout << "[" << idx[i] << "]" << endl;
            if(st.is_leaf(idx[i]))
                cout << endl;
            else{
                cout << "'a'" << st.child(idx[i], 'a');
                cout << "'b'" << st.child(idx[i], 'b');
                cout << "'#'" << st.child(idx[i], '#') << endl;
            }
        }

        for(size_type i = 0; i < 12; i++){
            cout << "[" << idx[i] << "]" << endl;
            cout << "'a'" << st.wl(idx[i], 'a');
            cout << "'b'" << st.wl(idx[i], 'b');
            cout << "'#'" << st.wl(idx[i], '#') << endl;
        }

        for(size_type i = 0; i < 12; i++)
            cout << "[" << idx[i] << "] -> " << st.sl(idx[i]) << endl;

        runs = build_runs(T, Sfwd, Bwtfwd, false);
        ms = build_ms(T, Srev, Bwtrev, false);
    }

    size_type operator[](size_type k){ return get_ms(ms, k); }

    void dump(){
        for (int i = 0; i < ms_size; i++)
            cout << (*this)[i];
        cout << endl;
    }
};

}

#endif /* fd_ms_h */
