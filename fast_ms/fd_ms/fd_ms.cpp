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

#include "fd_ms.hpp"


using namespace std;
using namespace fdms;


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


class MS{
private:
    sdsl::bit_vector _ms, ms_runs;

    size_type find_k_prim(size_type k){
        size_t total_zeros = rank_support_v<0>(&ms_runs)(ms_runs.size());
        size_t zeros = rank_support_v<0>(&ms_runs)(k + 1);

        if(total_zeros == zeros)
            return ms_runs.size();

        bit_vector::select_0_type b_sel(&ms_runs);
        return b_sel(zeros + 1);
    }

    size_type get_ms(bit_vector ms_, size_type k){
        if(k == -1)
            return (size_type) 1;
        return bit_vector::select_1_type (&ms_)(k + 1) - (2 * k);
    }

    sdsl::bit_vector build_ms(string& t, string& s, bit_vector& bp, const bool verbose){
        fdms::bp_support_sada<> bp_supp(&bp);
        Bwt bwt(s);
        Stree st(bp_supp, bwt);

        size_type k = 0, h_star = k + 1, k_prim, ms_idx = 0, ms_size = t.size() ;
        sdsl::bit_vector ms(ms_size * 2);

        uint8_t c = t[k];
        Interval I{bwt, static_cast<char>(c)};

        node_type v = st.child(st.root(), c); // stree node
        while(k < ms_size){
            output_partial_vec(ms, ms_idx, "ms", verbose);
            for(; !I.is_empty() && h_star < ms_size; ){
                c = t[h_star];
                I.bstep(c);
                if(!I.is_empty()){
                    v = st.wl(v, c);
                    h_star ++;
                }
            }
            for(int i = 0; i < h_star - k - get_ms(ms, k - 1) + 1; i++)
                ms[ms_idx++] = 0;
            if(h_star - k - get_ms(ms, k - 1) + 1 > 0)
                ms[ms_idx++] = 1;

            if(h_star < ms_size){
                do {
                    v = st.parent(v);
                    I.set(st.lb(v), st.rb(v));
                    I.bstep(t[h_star]);
                } while(I.is_empty());
                h_star = h_star + 1;
            }
            // k_prim: index of the first zero to the right of k in runs
            k_prim = find_k_prim(k);

            for(size_type i = k + 1; i <= k_prim - 1; i++)
                ms[ms_idx++] = 1;

            // update v
            v = st.wl(v, c);
            k = k_prim;
        }
        output_partial_vec(ms, ms_idx, "ms", verbose);
        return ms;
    }

public:
    size_type operator[](size_type k){ return get_ms(_ms, k); }

    void dump(){
        for (int i = 0; i < _ms.size() / 2; i++)
            cout << (*this)[i];
        cout << endl;
    }

    MS(string& t, string& s, bit_vector bp, bit_vector runs, const bool verbose) : ms_runs(runs) {
        _ms = build_ms(t, s, bp, verbose);
    }
};

bit_vector build_runs(string& t, string& s, bit_vector& bp, const bool verbose){
    fdms::bp_support_sada<> bp_supp(&bp);
    Bwt bwt(s);
    Stree st(bp_supp, bwt);

    size_type ms_size = t.size();
    bit_vector runs(ms_size);
    size_type k = ms_size, c = t[k - 1];
    Interval I{bwt, static_cast<char>(c)};

    node_type v = st.child(st.root(), c); // stree node
    while(--k > 0){
        c = t[k-1];
        I.bstep(c);
        if(I.is_empty()){
            runs[k] = 0;
            // update I to the parent of the proper locus of w until we can extend by 'c'
            do{
                v = st.parent(v);
                I.set(st.lb(v), st.rb(v));
                I.bstep(c);
            } while(I.is_empty());
        } else {
            runs[k] = 1;
        }
        v = st.wl(v, c); // update v
        output_partial_vec(runs, k, "runs", verbose);
    }
    return runs;
}


MS comp(InputSpec& T, InputSpec& S_fwd, InputSpec& S_rev, const bool verbose){
    string t = T.load_s();
    string s = S_fwd.load_s();
    bit_vector bp = S_fwd.load_bps();
    bit_vector runs = build_runs(t, s, bp, verbose);

    s = S_rev.load_s();
    bp = S_rev.load_bps();
    return MS(t, s, bp, runs, verbose);
}

Mstat compute_ms(string& T,
                 string& Sfwd, string& Sfwdbp,
                 string& Srev, string& Srevbp,
                 const bool verbose){

    bit_vector bfwd(Sfwdbp.size()), brev(Srevbp.size());
    for(size_type i=0; i<Sfwdbp.size(); i++)
        bfwd[i] = ((unsigned char)Sfwdbp[i] - 48);
    for(size_type i=0; i<Srevbp.size(); i++)
        brev[i] = ((unsigned char)Srevbp[i] - 48);

    fdms::bp_support_sada<> Bpsfwd(&bfwd);
    fdms::bp_support_sada<> Bpsrev(&brev);

    if(verbose){
        cst_sct3<> st_of_s, st_of_s_rev;
        construct_im(st_of_s, Sfwd, 1);
        //construct_im(st_of_s_rev, Srev, 1);

        cout << " i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << endl;
        csXprintf(cout, "%2I %2S %3s %3P %2p %3B   %:3T", st_of_s.csa);
        //cout << " i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << endl;
        //csXprintf(cout, "%2I %2S %3s %3P %2p %3B   %:3T", st_of_s_rev.csa);
    }

    Bwt Bwtfwd(Sfwd);
    Bwt Bwtrev(Srev);

    cout << "{compute_ms:" << endl;
    cout << "{bfwd__brev__Bpsfwd__Bpsrev: " << sdsl::size_in_bytes(bfwd) + sdsl::size_in_bytes(brev) + sdsl::size_in_bytes(Bpsfwd) + sdsl::size_in_bytes(Bpsrev) << "}, " << endl;
    Mstat MS(T,
             Sfwd, Bwtfwd, Bpsfwd,
             Srev, Bwtrev, Bpsrev,
             verbose);
    cout << "}" << endl;
    return MS;
}



int main(int argc, char **argv){
    if(argc == 6) {
        InputSpec tspec(argv[1], string(""));
        InputSpec sfwd_spec(argv[2], argv[3]);
        InputSpec srev_spec(argv[4], argv[5]);
        cout << argv[1] << " ";
        comp(tspec, sfwd_spec, srev_spec, 0).dump();
    }
    else {
        if(true){
            InputSpec tspec("/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/input_data/0t.txt", string(""));
            InputSpec sfwd_spec("/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/input_data/0s_fwd.txt",
                                "/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/input_data/0s_fwd_bp.txt");
            InputSpec srev_spec("/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/input_data/0s_rev.txt",
                                "/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/input_data/0s_rev_bp.txt");
            comp(tspec, sfwd_spec, srev_spec, 1);
        } else {
            string T {"baaaaaa"};
            string S {"bbbaababa"};
            string Sbp {"11011010110100011101010011010000"};
            string Srev {"ababaabbb"};
            string Srevbp {"1101101110100100011011010011010000"};

            cout << "{main:" << endl;
            cout << "{S__T__Sbp__Srev__Srevbp: " << S.size() + T.size() + Sbp.size() + Srevbp.size() << "}, " << endl;
            compute_ms(T, S, Sbp, Srev, Srevbp, 1).dump();
            cout << "}" << endl;
        }
    }
    return 0;
}

