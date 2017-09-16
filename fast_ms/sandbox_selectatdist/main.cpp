//
//  main.cpp
//  sandbox_selectatdist
//
//  Created by denas on 9/15/17.
//  Copyright © 2017 denas. All rights reserved.
//

#include <iostream>
#include "utils.hpp"
#include "cmd_utils.hpp"
#include "stree_sct3.hpp"

using namespace std;
using namespace fdms;


typedef typename StreeOhleb<>::node_type node_type;

void slow(StreeOhleb<>& st, string& s, const size_type max_k, const size_type cnt){
    node_type v = st.root();
    size_type res = 0;
    for (size_type k = max_k - 1; k > 0; k--){
        char c = s[k];
        v = st.double_rank_fail_wl(v, s[k]);
        
        size_type cnt_c = st.m_csa.C[st.m_csa.char2comp[c] + 1] - st.m_csa.C[st.m_csa.char2comp[c]];
        if(st.m_csa.bwt.rank(v.j + 1, c) < cnt_c - cnt + 1){
            auto start = timer::now();
            res = st.m_csa.wavelet_tree.select(st.m_csa.wavelet_tree.rank(v.j + 1, s[k]) + cnt, s[k]);
            size_type t = std::chrono::duration_cast<std::chrono::microseconds>(timer::now() - start).count();
            
            assert (res == st.m_csa.wavelet_tree.select_at_dist(c, v.j, cnt));
            cout << k << "," << s[k] << "," <<  v.j << "," << res << "," << t << ",s" << endl;
        }
    }
}

void fast(StreeOhleb<>& st, string& s, const size_type max_k, const size_type cnt){
    node_type v = st.root();
    size_type res = 0;
    for (size_type k = max_k; k > 0; k--){
        char c = s[k];
        v = st.double_rank_fail_wl(v, s[k]);
        
        size_type cnt_c = st.m_csa.C[st.m_csa.char2comp[c] + 1] - st.m_csa.C[st.m_csa.char2comp[c]];
        if(st.m_csa.bwt.rank(v.j + 1, c) < cnt_c - cnt + 1){
            auto start = timer::now();
            res = st.m_csa.wavelet_tree.select_at_dist(c, v.j, cnt);
            size_type t = std::chrono::duration_cast<std::chrono::microseconds>(timer::now() - start).count();

            assert (res == st.m_csa.wavelet_tree.select_at_dist(c, v.j, cnt));
            cout << k << "," << s[k] << "," <<  v.j << "," << res << "," << t << ",f" << endl;
        }
    }
}

int main(int argc, char** argv) {
    OptParser input(argc, argv);
    InputFlags flags(input);
    
    InputSpec s_spec("/Users/denas/projects/matching_statistics/indexed_ms/tests/datasets/big_paper2/rep_100000000s_dis_500000t_abcd_sim1000.s");
    flags.load_stree = true;
    //InputSpec s_spec(input.getCmdOption("-s_path"));
    string s = s_spec.load_s();
    StreeOhleb<> st;

    size_type st_time = fdms::load_st(st, s, s_spec.fwd_cst_fname, flags.load_stree);
    cerr << "DONE (" << st_time / 1000 << " seconds)" << endl;

    size_type max_k = (s.size() / 100) - 1;
    size_type cnt = 1;
    cout << "k,char,idx,res,time_micro,method" << endl;
    slow(st, s, max_k, 1);
    fast(st, s, max_k, 1);

/*
    sdsl::wt_huff<> wt;
    sdsl::construct(wt, "/Users/denas/projects/matching_statistics/indexed_ms/tests/datasets/big_paper2/rep_100000000s_dis_500000t_abcd_sim1000.s", 1);
    
    string alp{"ab"};
    sr(wt, alp, 100);
 */
    return 0;
}
