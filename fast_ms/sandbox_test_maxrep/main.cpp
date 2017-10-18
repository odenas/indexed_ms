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
#include "stree_sct3.hpp"
#include "cst_iterator.hpp"

#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))
#define IS_MAXIMAL(node) ( ((node).i != (node).j) && (maxrep[(node).i] == 1) && (maxrep[(node).j] == 1) )

using namespace std;
using namespace fdms;

typedef StreeOhleb<>::node_type node_type;
typedef pair<node_type, char> edge_type;


/*
 For every traversed node `alpha` try all Weiner links for all characters `c`.
 
 If at least 2 Weiner links are successful, add pairs (alpha,c) to the list of:
 - successful Weiner links for maximal nodes, where c is a char generating a successful Weiner link
 - unsuccessful Weiner links for maximal nodes, where c' is a character that generated a failed Weiner link.
 If just one character `c` gave a successful Weiner link, add the pair
 - (alpha,c) to the list of successful Weiner links for non maximal nodes
 - (alpha,c') to the list unsuccessful Weiner links for non maximal nodes, for all characters c' != c
 */
void fill_vectors(const StreeOhleb<>& st, vector<edge_type>& l1, vector<edge_type>& l2, vector<edge_type>& l3, vector<edge_type>& l4){
    cerr << "filling vectors ...";
    size_type nodes_visited = 0;
    sdsl::bit_vector wl_presence(st.csa.sigma);
    for(cst_dfslr_iterator<StreeOhleb<>> pos(&st); !pos.end(); ++pos, nodes_visited++){
        node_type currnode = *pos;
        
        for(size_type i = 0; i < st.csa.sigma; i++)
            wl_presence[i] = st.has_wl(currnode, st.csa.comp2char[i]);
        size_type wl_count = sdsl::util::cnt_one_bits(wl_presence);
        
        if(wl_count >= 2){ // maximal node
            for(size_type i = 0; i < st.csa.sigma; i++){
                (wl_presence[i] ? l1 : l2).push_back(make_pair(currnode, st.csa.comp2char[i]));
            }
        } else if (wl_count == 1){  // non-maximal node
            for(size_type i = 0; i < st.csa.sigma; i++){
                (wl_presence[i] ? l3 : l4).push_back(make_pair(currnode, st.csa.comp2char[i]));
            }
        } else {
            cerr << "grrr." << endl;
            exit(1);
        }
    }
    cerr << "DONE. " << nodes_visited << " visited." << endl;
}

void shuffle_vectors(vector<edge_type>& l1, vector<edge_type>& l2, vector<edge_type>& l3, vector<edge_type>& l4){
    cerr << "shuffling vectors ";
    std::random_shuffle(l1.begin(), l1.end());
    cerr << ".";
    std::random_shuffle(l2.begin(), l2.end());
    cerr << ".";
    std::random_shuffle(l3.begin(), l3.end());
    cerr << ".";
    std::random_shuffle(l4.begin(), l4.end());
    cerr << ". DONE." << endl;
}

size_type t1(const StreeOhleb<>& st, const vector<edge_type> vec, wl_method_t1 wl_f_ptr){
    node_type u = st.root();
    auto start = timer::now();
    for(auto e : vec)
        u = CALL_MEMBER_FN(st, wl_f_ptr)(e.first, e.second);
    size_type t = std::chrono::duration_cast<std::chrono::microseconds>(timer::now() - start).count();
    return t;
}

size_type t2(const StreeOhleb<>& st, const vector<edge_type> vec, bool maximal){
    node_type u = st.root();
    auto start = timer::now();
    for(auto e : vec)
        u = st.double_rank_fail_wl_mrep(e.first, e.second, maximal);
    size_type t = std::chrono::duration_cast<std::chrono::microseconds>(timer::now() - start).count();
    return t;
}


void time_all_methods(const StreeOhleb<>& st, const vector<edge_type> vec, const bool maximal, const bool successful){
    size_type t = 0;
    
    t = t1(st, vec, &StreeOhleb<>::single_rank_wl);
    cout << "s," << maximal << "," << successful << "," << vec.size() << "," << t << endl;
    
    t = t1(st, vec, &StreeOhleb<>::double_rank_fail_wl);
    cout << "f," << maximal << "," << successful << "," << vec.size() << "," << t << endl;
    
    t = t1(st, vec, &StreeOhleb<>::double_rank_nofail_wl);
    cout << "d," << maximal << "," << successful << "," << vec.size() << "," << t << endl;
    
    t = t2(st, vec, maximal);
    cout << "m," << maximal << "," << successful << "," << vec.size() << "," << t << endl;
}


int main(int argc, char **argv) {
    OptParser input(argc, argv);
    InputFlags flags(input);

    //InputSpec s_spec("/Users/denas/projects/matching_statistics/indexed_ms/tests/datasets/small_paper2/rep_1000000s_dis_200000t_abcd_sim1000.s");
    //flags.load_stree = true;
    //flags.load_maxrep = true;
    InputSpec s_spec(input.getCmdOption("-s_path"));
    string s = s_spec.load_s();
    reverse_in_place(s);
    StreeOhleb<> st;
    size_type st_time = fdms::load_st(st, s, s_spec.rev_cst_fname, flags.load_stree);
    cerr << "DONE (" << st_time / 1000 << " seconds)" << endl;

    sdsl::bit_vector maxrep(s.size() + 1); sdsl::util::set_to_value(maxrep, 0);
    size_type mrep_time = fdms::load_maxrep(maxrep, st, s, s_spec.rev_maxrep_fname, flags.load_maxrep);
    cerr << "DONE ( " << mrep_time << " milliseconds)" << endl;

    vector<edge_type> l1; // successful Wl, maximal
    vector<edge_type> l2; // un-successful Wl, maximal
    vector<edge_type> l3; // successful Wl, non-maximal
    vector<edge_type> l4; // un-successful Wl, non-maximal

    fill_vectors(st, l1, l2, l3, l4);
    shuffle_vectors(l1, l2, l3, l4);
    
    cerr << "timing ";
    cout << "method,maximal,has_wl,ncalls,t_micro" << endl;
    time_all_methods(st, l1, true, true);
    cerr << ".";
    time_all_methods(st, l2, true, false);
    cerr << ".";
    time_all_methods(st, l3, false, true);
    cerr << ".";
    time_all_methods(st, l4, false, false);
    cerr << "DONE." << endl;
    return 0;
}
