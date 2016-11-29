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

#include "CmdArguments.h"
#include "fd_ms.hpp"


using namespace std;
using namespace fdms;
using timer = std::chrono::high_resolution_clock;



void output_partial_vec(const bit_vector& v, size_type idx, const char name[], bool verbose){
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

void dump_ms(bit_vector& ms){
    auto get_ms = [] (bit_vector& __ms, size_type __k) -> size_type {
        if(__k == -1)
            return (size_type) 1;
        return bit_vector::select_1_type (&__ms)(__k + 1) - (2 * __k);
    };

    for (size_type i = 0; i < ms.size() / 2; i++)
        cout << get_ms(ms, i);
    cout << endl;
}

void build_ms(const string& prefix, string& t, InputSpec& S_rev,
              bit_vector& runs, bit_vector& ms,
              const bool space_usage, const bool verbose){
    auto get_ms = [] (select_support_mcl<1,1>& __ms_select1, size_type __k) -> size_type {
        if(__k == -1)
            return (size_type) 1;
        return __ms_select1(__k + 1) - (2 * __k);
    };

    auto find_k_prim = []  (size_type __k, bit_vector& __runs, rank_support_v<0>& __runs_rank0, select_support_mcl<0, 1> __runs_sel0) -> size_type {
        size_t zeros = __runs_rank0(__k + 1);
        //return (__runs_rank0(__runs.size()) == zeros ? __runs.size() : __runs_sel0(zeros + 1));

        if(__runs_rank0(__runs.size()) == zeros)
            return __runs.size();
        return __runs_sel0(zeros + 1);
    };

    string s = S_rev.load_s();
    bit_vector bp = S_rev.load_bps();
    fdms::bp_support_sada<> bp_supp(&bp);
    Bwt bwt(s);
    Stree st(bp_supp, bwt);
    rank_support_v<0> runs_rank0(&runs);
    select_support_mcl<0, 1> runs_select0(&runs);
    size_type size_in_bytes_ms_select1 = 0;


    size_type k = 0, h_star = k + 1, k_prim, ms_idx = 0, ms_size = t.size() ;
    uint8_t c = t[k];
    Interval I{bwt, static_cast<char>(c)};

    node_type v = st.wl(st.root(), c); // stree node
    while(k < ms_size){
        output_partial_vec(ms, ms_idx, "ms", verbose);

        select_support_mcl<1,1> ms_select1(&ms);
        size_in_bytes_ms_select1 = (size_in_bytes_ms_select1 < sdsl::size_in_bytes(ms_select1) ? sdsl::size_in_bytes(ms_select1) : size_in_bytes_ms_select1);

        for(; !I.is_empty() && h_star < ms_size; ){
            c = t[h_star];
            I.bstep(c);
            if(!I.is_empty()){
                v = st.wl(v, c);
                h_star ++;
            }
        }
        for(int i = 0; i < h_star - k - get_ms(ms_select1, k - 1) + 1; i++)
            ms[ms_idx++] = 0;
        if(h_star - k - get_ms(ms_select1, k - 1) + 1 > 0)
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
        k_prim = find_k_prim(k, runs, runs_rank0, runs_select0);

        for(size_type i = k + 1; i <= k_prim - 1; i++)
            ms[ms_idx++] = 1;

        // update v
        v = st.wl(v, c);
        k = k_prim;
    }
    if(space_usage){
        cout << prefix << ", 2, build_ms, space, byte, s, " << s.size() << endl;
        cout << prefix << ", 2, build_ms, space, byte, s_bwt_wtree, " << bwt.size_in_bytes__wtree << endl;
        cout << prefix << ", 2, build_ms, space, byte, s_bwt_bwt, " << bwt.size_in_bytes__bwt << endl;
        cout << prefix << ", 2, build_ms, space, byte, s_bwt_alp, " << bwt.size_in_bytes__C + bwt.size_in_bytes__char2int + bwt.size_in_bytes__Sigma << endl;
        cout << prefix << ", 2, build_ms, space, byte, s_stree_bp, " << sdsl::size_in_bytes(bp) << endl;
        cout << prefix << ", 2, build_ms, space, byte, s_stree_bpsupp, " << sdsl::size_in_bytes(bp_supp) << endl;
        cout << prefix << ", 2, build_ms, space, byte, s_stree_bpsupp_select, " << st.size_in_bytes__select << endl;
        cout << prefix << ", 2, build_ms, space, byte, s_stree_bpsupp_rank, " << st.size_in_bytes__rank << endl;
        cout << prefix << ", 2, build_ms, space, byte, runs_rank, " << sdsl::size_in_bytes(runs_rank0) << endl;
        cout << prefix << ", 2, build_ms, space, byte, runs_select, " << sdsl::size_in_bytes(runs_select0) << endl;
        cout << prefix << ", 2, build_ms, space, byte, ms_select, " << size_in_bytes_ms_select1 << endl;
    }
}

void build_runs(const string& prefix, string& t, InputSpec& S_fwd,
                bit_vector& runs,
                const bool space_usage, const bool verbose){

    string s = S_fwd.load_s();
    Bwt bwt(s);

    bit_vector bp = S_fwd.load_bps();
    fdms::bp_support_sada<> bp_supp(&bp);
    Stree st(bp_supp, bwt);

    size_type ms_size = t.size();
    size_type k = ms_size, c = t[k - 1];
    Interval I{bwt, static_cast<char>(c)};

    node_type v = st.wl(st.root(), c); // stree node
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
    if(space_usage){
        cout << prefix << ", 2, build_runs, space, byte, s, " << s.size() << endl;
        cout << prefix << ", 2, build_runs, space, byte, s_bwt_wtree, " << bwt.size_in_bytes__wtree << endl;
        cout << prefix << ", 2, build_runs, space, byte, s_bwt_bwt, " << bwt.size_in_bytes__bwt << endl;
        cout << prefix << ", 2, build_runs, space, byte, s_bwt_alp, " << bwt.size_in_bytes__C + bwt.size_in_bytes__char2int + bwt.size_in_bytes__Sigma << endl;
        cout << prefix << ", 2, build_runs, space, byte, s_stree_bp, " << sdsl::size_in_bytes(bp) << endl;
        cout << prefix << ", 2, build_runs, space, byte, s_stree_bpsupp, " << sdsl::size_in_bytes(bp_supp) << endl;
        cout << prefix << ", 2, build_runs, space, byte, s_stree_bpsupp_select, " << st.size_in_bytes__select << endl;
        cout << prefix << ", 2, build_runs, space, byte, s_stree_bpsupp_rank, " << st.size_in_bytes__rank << endl;
    }
}


void comp(const string& prefix, InputSpec& T, InputSpec& S_fwd, InputSpec& S_rev,
          const bool space_usage, const bool time_usage, const bool answer, const bool verbose){
    string t = T.load_s();
    bit_vector runs(t.size());
    bit_vector ms(t.size() * 2);

    if(space_usage){
        cout << "prefix, level, func, measuring, unit, item, value" << endl;
        cout << prefix << ", 1, comp, space, byte, t, " << t.size() << endl;
        cout << prefix << ", 1, comp, space, byte, runs, " << sdsl::size_in_bytes(runs) << endl;
        cout << prefix << ", 1, comp, space, byte, ms, " << sdsl::size_in_bytes(ms) << endl;
    }

    auto runs_start = timer::now();
    build_runs(prefix, t, S_fwd, runs, space_usage, verbose);
    auto runs_stop = timer::now();

    auto ms_start = timer::now();
    build_ms(prefix, t, S_rev, runs, ms, space_usage, verbose);
    auto ms_stop = timer::now();

    if(time_usage){
        cout << prefix << ", 1, comp, time, milliseconds, runs, " << std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count() << endl;
        cout << prefix << ", 1, comp, time, milliseconds, ms, " << std::chrono::duration_cast<std::chrono::milliseconds>(ms_stop - ms_start).count() << endl;
    }

    if(answer){
        cout << prefix << " ";
        for(size_type i = 0; i < ms.size() / 2; i++)
            cout << bit_vector::select_1_type (&ms)(i + 1) - (2 * i);
        cout << endl;
    }
}


int main(int argc, char **argv){
    InputParser input(argc, argv);
    if(argc == 1){
        /**
         cst_sct3<> st_of_s, st_of_s_rev;
         construct_im(st_of_s, Sfwd, 1);
         //construct_im(st_of_s_rev, Srev, 1);

         cout << " i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << endl;
         csXprintf(cout, "%2I %2S %3s %3P %2p %3B   %:3T", st_of_s.csa);
         */

        string prefix {"none"};
        InputSpec tspec("/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/input_data/0t.txt", string(""));
        InputSpec sfwd_spec("/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/input_data/0s_fwd.txt",
                            "/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/input_data/0s_fwd_bp.txt");
        InputSpec srev_spec("/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/input_data/0s_rev.txt",
                            "/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/input_data/0s_rev_bp.txt");
        comp(prefix, tspec, sfwd_spec, srev_spec, 0, 0, 1, 1);
    } else {
        const string& base_dir = input.getCmdOption("-d");
        const string& prefix = input.getCmdOption("-p");
        // flags
        const string& space_usage = input.getCmdOption("-s"); // space usage
        const string& time_usage = input.getCmdOption("-t"); // time usage

        const string& answer = input.getCmdOption("-a"); // answer
        const string& verbose = input.getCmdOption("-v"); // verbose


        InputSpec tspec(base_dir + "/" + prefix + "t.txt", string(""));
        InputSpec sfwd_spec(base_dir + "/" + prefix + "s_fwd.txt", base_dir + "/" + prefix + "s_fwd_bp.txt");
        InputSpec srev_spec(base_dir + "/" + prefix + "s_rev.txt", base_dir + "/" + prefix + "s_rev_bp.txt");
        comp(prefix, tspec, sfwd_spec, srev_spec, space_usage == "1", time_usage == "1", answer == "1", verbose == "1");
    }
    return 0;
}

