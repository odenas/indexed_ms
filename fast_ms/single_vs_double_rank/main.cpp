//
//  main.cpp
//  single_vs_double_rank
//
//  Created by denas on 1/28/17.
//  Copyright © 2017 denas. All rights reserved.
//

#include <iostream>
#include "utils.hpp"
#include "fd_ms.hpp"
#include "stree_sct3.hpp"
#define NWLCALLS 1000000
#define NTRIALS  5
#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))
#define IS_MAXIMAL(node) ( ((node).i != (node).j) && (maxrep[(node).i] == 1) && (maxrep[(node).j] == 1) )


using namespace std;
using namespace fdms;

typedef typename StreeOhleb<>::node_type node_type;


void rank_and_check(const StreeOhleb<> st, const sdsl::bit_vector& maxrep, const size_type ncalls){
    node_type v = st.rightmost_leaf(st.root());
    node_type u = v;
    
    for(size_type nc = 0; nc < ncalls; nc++){
        for(size_type cidx = 0; cidx < st.csa.sigma; cidx++){
            char_type c = st.csa.comp2char[cidx];
            bool is_maximal = IS_MAXIMAL(v);
            
            if(is_maximal){
                u = st.double_rank_fail_wl(v, c);
            } else {
                size_type c_right = st.csa.bwt.rank_and_check(v.j + 1, c);
                if(c_right == 0){
                    u = st.root();
                } else {
                    std::pair<size_type, size_type> lr = std::make_pair(c_right - (v.j - v.i + 1), c_right);
                    u = st._wl_from_interval(lr, c);
                }
            }
            if(!st.is_root(u)){
                v = u;
                break;
            }
        }
    }
}

void single_rank(const StreeOhleb<> st, const sdsl::bit_vector& maxrep, const size_type ncalls){
    node_type v = st.rightmost_leaf(st.root());
    node_type u = v;
    
    for(size_type nc = 0; nc < ncalls; nc++){
        for(size_type cidx = 0; cidx < st.csa.sigma; cidx++){
            char_type c = st.csa.comp2char[cidx];
            bool is_maximal = IS_MAXIMAL(v);
        
            if(is_maximal){
                u = st.double_rank_fail_wl(v, c);
            } else {
                if (st.csa.bwt[v.j] != c){
                    u = st.root();
                } else {
                    size_type c_right = st.csa.bwt.rank(v.j + 1, c);
                    std::pair<size_type, size_type> lr = std::make_pair(c_right - (v.j - v.i + 1), c_right);
                    u = st._wl_from_interval(lr, c);
                }
            }
            if(!st.is_root(u)){
                v = u;
                break;
            }
        }
    }
}

void double_rank(const StreeOhleb<> st, const sdsl::bit_vector& maxrep, const size_type ncalls){
    node_type v = st.rightmost_leaf(st.root());
    node_type u = v;
    for(size_type nc = 0; nc < ncalls; nc++){
        for(size_type cidx = 0; cidx < st.csa.sigma; cidx++){
            char_type c = st.csa.comp2char[cidx];
            u = st.double_rank_fail_wl(v, c);
            if(!st.is_root(u)){
                v = u;
                break;
            }
        }
    }
}


size_type count_wl(const StreeOhleb<>& st, const node_type v){
    size_type wl_count = 0;
    for(size_type cidx = 0; cidx < st.csa.sigma; cidx++){
        if(!st.is_root(st.double_rank_fail_wl(v, st.csa.comp2char[cidx])))
            wl_count += 1;
    }
    return wl_count;
}

bool in_same_major_block(const node_type v) { return ((v.i>>8) == (v.j>>8)); }

node_type find_small_node(const StreeOhleb<>& st, const string s){
    node_type v = st.rightmost_leaf(st.root());
    for(size_type i = 0; i < 4; i++){
        cerr << "(" << v.i << "," << v.j << ") : " << in_same_major_block(v) << endl;
        v = st.parent(v);
    }
    
    assert(v.i != v.j);
    
    size_type wl_cnt = 0;
    cerr << "(" << v.i << "," << v.j << ") : " << in_same_major_block(v) << endl;
    while(!st.is_root(v)){
        wl_cnt = count_wl(st, v);
        assert (!st.is_root(v));
        assert (wl_cnt > 0);
        
        cerr << "(" << v.i << "," << v.j << ") : " << in_same_major_block(v) << " : " << wl_cnt << endl;
        if(wl_cnt > 1 && wl_cnt < st.csa.sigma - 1 && in_same_major_block(v))
            return v;
        
        if (in_same_major_block(v)){
            v = st.sl(v);
        } else {
            for(size_type cidx = 0; cidx < st.csa.sigma; cidx++){
                node_type u = st.double_rank_fail_wl(v, st.csa.comp2char[cidx]);
                if(!st.is_root(u)){
                    v = u;
                    break;
                }
            }
        }
    }
    return st.root(); // not found
}


/* small = in same major block */
node_type find_node(const StreeOhleb<>& st, const string s, const bool small){
    node_type v = st.root();
    size_type i = s.size() - 1;
    size_type wl_cnt = 0;
    
    auto block_condition = [] (node_type v, bool sm) -> bool
    {
        return (sm ? in_same_major_block(v) : !in_same_major_block(v));
    };
    
    i -= 1;
    v = st.single_rank_wl(v, s[i]);
    
    while(--i > 0){
        v = st.single_rank_wl(v, s[i]);
        wl_cnt = count_wl(st, v);
        assert (!st.is_root(v));
        assert (wl_cnt > 0);
        
        if(wl_cnt > 1 && wl_cnt < st.csa.sigma - 1 && block_condition(v, small))
            return v;
        cerr << i << endl;
    }
    return st.root(); // not found
}

void dump_report_slice(StreeOhleb<>& st, string& s, size_type ntrial, size_type close, string method,
                       std::vector<size_type>& time_usage, std::vector<bool>& has_wl){
    // cout << "len_s,nwlcalls,ntrial,close,method,char,charidx,time_ms" << endl;
    for(size_type i = 0; i < time_usage.size(); i++){
        cout << s.size() << "," ;  // s_len
        cout << NWLCALLS << "," ;  // nwlcalls
        cout << ntrial << "," ;    // ntrial
        cout << close << "," ;     // is close?
        cout << method << ",";     // method name
        cout << (char) (i == 0 ? '#' : st.csa.comp2char[i]) << ","; // char
        cout << i << "," ;         // char idx
        cout << has_wl[i] << "," ;         // has wl?
        cout << time_usage[i] << endl;
    }
}



int main(int argc, char **argv) {
    OptParser input(argc, argv);
    InputFlags flags(input);
    
    InputSpec s_spec("/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/datasets/testing/rnd_200_32768.s");
    flags.load_stree = false;
    flags.use_maxrep = true;
    //InputSpec s_spec(input.getCmdOption("-s_path"));
    string s = s_spec.load_s();
    StreeOhleb<> st;
    size_type st_time = fdms::load_st(st, s, s_spec.fwd_cst_fname, flags.load_stree);
    cerr << "DONE (" << st_time / 1000 << " seconds)" << endl;

    
    sdsl::bit_vector maxrep(s.size());
    size_type maxrep_time = load_maxrep(maxrep, st, s, s_spec.rev_maxrep_fname, flags.use_maxrep);
    cerr << "DONE (" << maxrep_time / 1000 << " seconds)" << endl;

    
    auto start_time = timer::now();
    rank_and_check(st, maxrep, NWLCALLS);
    size_type rank_and_check_time = std::chrono::duration_cast<std::chrono::milliseconds>(timer::now() - start_time).count();
    cout << "rank_and_check," << rank_and_check_time << endl;

    start_time = timer::now();
    single_rank(st, maxrep, NWLCALLS);
    size_type single_rank_time = std::chrono::duration_cast<std::chrono::milliseconds>(timer::now() - start_time).count();
    cout << "single_rank," << single_rank_time << endl;

    start_time = timer::now();
    double_rank(st, maxrep, NWLCALLS);
    size_type double_rank_time = std::chrono::duration_cast<std::chrono::milliseconds>(timer::now() - start_time).count();
    cout << "double_rank," << double_rank_time << endl;

    return 0;
}
