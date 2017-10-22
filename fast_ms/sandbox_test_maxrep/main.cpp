//
//  main.cpp
//  sandbox_test_maxrep
//
//  Created by denas on 5/20/17.
//  Copyright © 2017 denas. All rights reserved.
//


#include <iostream>
#include <fstream>
#include <vector>

#include "input_spec.hpp"
#include "opt_parser.hpp"

#include "stree_sct3.hpp"
#include "cst_iterator.hpp"
#include "edge_list.hpp"

#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))

using namespace std;
using namespace fdms;


//Declare various wl strategies
typedef StreeOhleb<>::node_type (StreeOhleb<>::*wl_method_t1) (const StreeOhleb<>::node_type& v, const StreeOhleb<>::char_type c) const;
typedef StreeOhleb<>::node_type (StreeOhleb<>::*wl_method_t2) (const StreeOhleb<>::node_type& v, const StreeOhleb<>::char_type c, const bool is_max) const;
typedef std::pair<StreeOhleb<>::size_type, StreeOhleb<>::size_type> (sdsl::bwt_of_csa_wt<sdsl::csa_wt<>>::*double_rank_method)(const StreeOhleb<>::size_type i, const StreeOhleb<>::size_type j, const StreeOhleb<>::char_type c)const;
typedef StreeOhleb<>::node_type (StreeOhleb<>::*parent_seq_method) (const StreeOhleb<>::node_type& v, const StreeOhleb<>::char_type c) const;


class InputFlags{
public:
    bool load_cst, load_maxrep;
    
    InputFlags(){}
    
    InputFlags(const bool load_cst, const bool load_maxrep) : load_cst{load_cst} {}
    
    InputFlags(const InputFlags& i) : load_cst{i.load_cst}, load_maxrep{i.load_maxrep} {}
    
    InputFlags(const OptParser args){
        load_cst = (args.getCmdOption("-load_cst") == "1");
        load_maxrep = (args.getCmdOption("-load_maxrep") == "1");
    }
};



size_type t1(const StreeOhleb<>& st, const vector<edge_type> vec, wl_method_t1 wl_f_ptr){
    node_type u = st.root();
    auto start = timer::now();
    for(auto e : vec)
        u = CALL_MEMBER_FN(st, wl_f_ptr)(e.m_node, e.m_c);
    size_type t = std::chrono::duration_cast<std::chrono::microseconds>(timer::now() - start).count();
    return t;
}

size_type t2(const StreeOhleb<>& st, const vector<edge_type> vec, bool maximal){
    node_type u = st.root();
    auto start = timer::now();
    for(auto e : vec)
        u = st.double_rank_fail_wl_mrep(e.m_node, e.m_c, maximal);
    size_type t = std::chrono::duration_cast<std::chrono::microseconds>(timer::now() - start).count();
    return t;
}


void time_all_methods(const StreeOhleb<>& st, EdgeList elst, const bool maximal, const bool successful){
    size_type n=20;
    vector<edge_type> vec = elst.vec(maximal, successful);
    //size_type s=0, f=0, d=0, m=0;
    for (int i=0; i<n; i++){
        std::random_shuffle(vec.begin(), vec.end());
        size_type s = t1(st, vec, &StreeOhleb<>::single_rank_wl);
        size_type f = t1(st, vec, &StreeOhleb<>::double_rank_fail_wl);
        size_type d = t1(st, vec, &StreeOhleb<>::double_rank_nofail_wl);
        size_type m = t2(st, vec, maximal);
        
        cout << i << ",s," << maximal << "," << successful << "," << vec.size() << "," << s << endl;
        cout << i << ",f," << maximal << "," << successful << "," << vec.size() << "," << f << endl;
        cout << i << ",d," << maximal << "," << successful << "," << vec.size() << "," << d << endl;
        cout << i << ",m," << maximal << "," << successful << "," << vec.size() << "," << m << endl;
    }
}


int main(int argc, char **argv) {
    OptParser input(argc, argv);
    InputFlags flags(input);

    //InputSpec s_spec("/Users/denas/projects/matching_statistics/indexed_ms/tests/datasets/small_paper2/rep_1000000s_dis_200000t_abcd_sim1000.s");
    //flags.load_stree = true;
    //flags.load_maxrep = true;
    InputSpec s_spec(input.getCmdOption("-s_path"));
    string s = s_spec.load_s();
    InputSpec::reverse_in_place(s);

    StreeOhleb<> st;
    size_type st_time = load_or_build(st, s, s_spec.rev_cst_fname, flags.load_cst);
    cerr << "DONE (" << st_time / 1000 << " seconds)" << endl;

    //sdsl::bit_vector maxrep(s.size() + 1); sdsl::util::set_to_value(maxrep, 0);
    Maxrep maxrep;
    size_type mrep_time = Maxrep::load_or_build(maxrep, st, s_spec.rev_maxrep_fname, flags.load_maxrep);
    cerr << "DONE ( " << mrep_time << " milliseconds)" << endl;
    
    EdgeList elst = load_edge_list_bin(s_spec.rev_elst_fname);
    elst.shuffle_vectors();
    
    cerr << "timing ";
    cout << "ntrial,method,maximal,has_wl,ncalls,t_micro" << endl;
    time_all_methods(st, elst, true, true);
    cerr << ".";
    time_all_methods(st, elst, true, false);
    cerr << ".";
    time_all_methods(st, elst, false, true);
    cerr << ".";
    time_all_methods(st, elst, false, false);
    cerr << "DONE." << endl;
    return 0;
}
