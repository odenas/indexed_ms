//
//  main.cpp
//  sandbox_selectatdist
//
//  Created by denas on 9/15/17.
//  Copyright © 2017 denas. All rights reserved.
//

#include <iostream>
#include <fstream>

#include "fd_ms/opt_parser.hpp"
#include "fd_ms/input_spec.hpp"
#include "fd_ms/stree_sct3.hpp"

using namespace std;
using namespace fdms;


typedef StreeOhleb<> cst_t;
typedef typename cst_t::node_type node_type;
typedef typename cst_t::size_type size_type;
using timer = std::chrono::high_resolution_clock;


class InputFlags{
public:
    bool load_cst, load_maxrep;

    InputFlags(){}
    
    InputFlags(const bool load_cst, const bool load_maxrep) : load_cst{load_cst}, load_maxrep{load_maxrep} {}

    InputFlags(const InputFlags& i) : load_cst{i.load_cst}, load_maxrep{i.load_maxrep} {}

    InputFlags(const OptParser args){
        load_cst = (args.getCmdOption("-load_cst") == "1");
        load_maxrep = (args.getCmdOption("-load_maxrep") == "1");
    }

    void help(const string exec_name) const {
        cerr << exec_name << endl << "\t";
        cerr << "[-load_cst 0/1]" << endl << "\t";
        cerr << "[-load_maxrep 0/1]" << endl << "\t";
        cerr << "-s_path <s path>" << endl;
    }
};

void show(char c, size_type idx, size_type res, size_type t, string method){
    cout << c << "," <<  idx << "," << res << "," << t << "," << method << endl;
}

void slow(cst_t& st, string& s, const size_type max_k, const size_type cnt){
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
            show(s[k], v.j, res, t, "s");
        }
    }
}

void fast(cst_t& st, string& s, const size_type max_k, const size_type cnt){
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
            show(s[k], v.j, res, t, "f");
        }
    }
}

int main(int argc, char** argv) {
    OptParser input(argc, argv);
    InputFlags flags(input);
    
    //InputSpec s_spec("/Users/denas/projects/matching_statistics/indexed_ms/tests/datasets/big_paper2/rep_100000000s_dis_500000t_abcd_sim1000.s");
    //flags.load_stree = true;

    InputSpec s_spec(input.getCmdOption("-s_path"));
    string s = s_spec.load_s();
    cst_t st;

    size_type st_time = load_or_build(st, s, s_spec.fwd_cst_fname, flags.load_cst);
    cerr << "DONE (" << st_time / 1000 << " seconds)" << endl;

    size_type max_k = (s.size() / 100) - 1;
    size_type cnt = 1;

    cout << "char,idx,res,time_micro,method" << endl;
    slow(st, s, max_k, cnt);
    fast(st, s, max_k, cnt);

    return 0;
}
