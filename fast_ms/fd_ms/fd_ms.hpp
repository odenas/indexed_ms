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

class InputSpec{
private:
    bit_vector parse_bitstr(string& s){
        bit_vector b(s.size());

        for(size_type i = 0; i < s.size(); i++)
            b[i] = ((unsigned char)s[i] - 48);
        return b;
    }
public:
    string s_fname, sbp_fname;

    InputSpec(string s_fn, string sbp_fn) : s_fname(s_fn), sbp_fname(sbp_fn){}

    string load_s(){
        string s;
        std::ifstream s_file {s_fname};
        while(s_file >> s)
            ;
        return s;
    }

    bit_vector load_bps(){
        string sbp;
        std::ifstream sbp_file {sbp_fname};
        sbp_file >> sbp;
        return parse_bitstr(sbp);
    }
};


class Interval{
private:
    Bwt& bwt;

public:
    size_type lb, ub;

    Interval(Bwt& bwt_, char c) : bwt{bwt_} {
        lb = bwt.C[bwt.char2int[c]];
        ub = bwt.C[bwt.char2int[c] + 1] - 1;
    }

    void set(size_type l, size_type u){
        lb = l;
        ub = u;
    }

    void bstep(char c){
        int cc = bwt.char2int[c];
        lb = bwt.C[cc] + bwt.rank(lb, c);
        ub = bwt.C[cc] + bwt.rank(ub + 1, c) - 1;
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

    bit_vector build_runs(string& t, Stree& st_of_s, Bwt& bwt, const bool verbose){
        bit_vector runs(ms_size);
        size_type k = ms_size, c = t[k - 1];
        Interval I{bwt, static_cast<char>(c)};

        node_type v = st_of_s.wl(st_of_s.root(), c); // stree node
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

    bit_vector build_ms(string& t, Stree& st_of_s, Bwt& bwt, const bool verbose){
        bit_vector ms(ms_size * 2);
        size_type k = 0, h_star = k + 1, k_prim, ms_idx = 0;
        uint8_t c = t[k];
        Interval I{bwt, static_cast<char>(c)};

        node_type v = st_of_s.wl(st_of_s.root(), c); // stree node

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

        cout << "{Mstat:" << endl;
        Stree st(bfwdbps, Bwtfwd);
        Stree strev(brevbps, Bwtrev);
        cout << "{st__strev: " << sdsl::size_in_bytes(st.m_bp_rank10) + sdsl::size_in_bytes(st.m_bp_rank10) + sdsl::size_in_bytes(strev.m_bp_rank10) + sdsl::size_in_bytes(strev.m_bp_rank10) << "}, " << endl;
        runs = build_runs(T, st, Bwtfwd, verbose);
        ms = build_ms(T, strev, Bwtrev, verbose);
        cout << "}" << endl;
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
