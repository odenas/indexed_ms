//
//  main.cpp
//  sandbox_test_maxrep
//
//  Created by denas on 5/20/17.
//  Copyright © 2017 denas. All rights reserved.
//


#include <iostream>
#include "utils.hpp"
#include "cmd_utils.hpp"
#include "maxrep_construction.hpp"
#include "stree_sct3.hpp"


#define NWLCALLS 100000
#define NTRIALS  5
#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))
#define IS_MAXIMAL(node) ( ((node).i != (node).j) && (maxrep[(node).i] == 1) && (maxrep[(node).j] == 1) )

using namespace std;
using namespace fdms;

typedef typename StreeOhleb<>::node_type node_type;

size_type count_wl(const StreeOhleb<>& st, const node_type v){
    size_type wl_count = 0;
    for(size_type cidx = 0; cidx < st.csa.sigma; cidx++){
        if(st.is_root(st.single_rank_wl(v, st.csa.comp2char[cidx])))
            wl_count += 1;
    }
    return wl_count;
}

bool in_same_major_block(const node_type v) { return ((v.i>>8) == (v.j>>8)); }

void show(const StreeOhleb<>& st, const sdsl::bit_vector& maxrep, const string fun_name,
          const node_type v, const size_type t_micro){
    
    cout << (in_same_major_block(v) ? 1 : 0) << "," ; // is close?
    cout << fun_name << ",";                // method name
    cout << count_wl(st, v) << "," ;        // how many wl?
    cout << IS_MAXIMAL(v) << ",";            // is_maximal
    cout << t_micro << endl;
}


void test_single_rank(const StreeOhleb<>& st,  string& s, const sdsl::bit_vector& maxrep, const size_type max_k){
    node_type v = st.root(), u = st.root();

    auto main_start = timer::now();
    for (size_type k = max_k - 1; k > 0; k--){
        auto start = timer::now();
        for(size_type i = 0; i < st.csa.sigma; i++)
            u = st.single_rank_wl(v, st.csa.comp2char[i]);
        size_type t = std::chrono::duration_cast<std::chrono::microseconds>(timer::now() - start).count();
        show(st, maxrep, "s", v, t);
        v = st.single_rank_wl(v, s[k]);
        if (k % 500000 == 0)
            report_progress(main_start, k, max_k);
    }
}

void test_double_rank(const StreeOhleb<>& st, const string& s, const sdsl::bit_vector& maxrep, const size_type max_k){
    node_type v = st.root(), u = st.root();
    
    //auto main_start = timer::now();
    for (size_type k = max_k - 1; k > 0; k--){
        auto start = timer::now();
        for(size_type i = 0; i < st.csa.sigma; i++)
            u = st.double_rank_fail_wl(v, st.csa.comp2char[i]);
        size_type t = std::chrono::duration_cast<std::chrono::microseconds>(timer::now() - start).count();
        show(st, maxrep, "f", v, t);
        v = st.double_rank_fail_wl(v, s[k]);
        //if(k % 50000 == 0)
        //    report_progress(main_start, k, max_k);
    }
}

void test_double_rank_nofail(const StreeOhleb<>& st, const string& s, const sdsl::bit_vector& maxrep, const size_type max_k){
    node_type v = st.root(), u = st.root();
    
    for (size_type k = max_k - 1; k > 0; k--){
        auto start = timer::now();
        for(size_type i = 0; i < st.csa.sigma; i++)
            u = st.double_rank_nofail_wl(v, st.csa.comp2char[i]);
        size_type t = std::chrono::duration_cast<std::chrono::microseconds>(timer::now() - start).count();
        show(st, maxrep, "d", v, t);
        v = st.double_rank_nofail_wl(v, s[k]);
    }
}

void test_mrep_rank(const StreeOhleb<>& st, const string& s, const sdsl::bit_vector& maxrep, const size_type max_k){
    node_type v = st.root(), u = st.root();
    
    for (size_type k = max_k - 1; k > 0; k--){
        auto start = timer::now();
        for(size_type i = 0; i < st.csa.sigma; i++)
            u = st.double_rank_fail_wl_mrep(v, st.csa.comp2char[i], IS_MAXIMAL(v));
        size_type t = std::chrono::duration_cast<std::chrono::microseconds>(timer::now() - start).count();
        show(st, maxrep, "m", v, t);
        v = st.double_rank_fail_wl_mrep(v, s[k], IS_MAXIMAL(v));
    }
}

int main(int argc, char **argv) {
    OptParser input(argc, argv);
    InputFlags flags(input);

    //InputSpec s_spec("/Users/denas/projects/matching_statistics/indexed_ms/tests/datasets/testing/rnd_200_256.s");
    InputSpec s_spec("/Users/denas/projects/matching_statistics/indexed_ms/tests/datasets/big_paper2/rep_100000000s_dis_500000t_abcdefghijklmnopqrst_sim1000.s");
    flags.load_stree = true;
    flags.load_maxrep = true;
    //InputSpec s_spec(input.getCmdOption("-s_path"));
    string s = s_spec.load_s();
    StreeOhleb<> st;
    size_type st_time = fdms::load_st(st, s, s_spec.fwd_cst_fname, flags.load_stree);
    cerr << "DONE (" << st_time / 1000 << " seconds)" << endl;


    sdsl::bit_vector maxrep(s.size() + 1);
    size_type mrep_time = fdms::load_maxrep(maxrep, st, s, s_spec.rev_maxrep_fname, flags.load_maxrep);
    cerr << "DONE ( " << mrep_time << " milliseconds)" << endl;

    size_type max_k = st.size() / 1000;
    cout << "close,method,has_wl,is_maximal,time_micro" << endl;
    //test_single_rank(st, s, maxrep, max_k);
    //test_double_rank(st, s, maxrep, max_k);
    //test_double_rank_nofail(st, s, maxrep, max_k);
    test_mrep_rank(st, s, maxrep, max_k);
    
    return 0;
}
