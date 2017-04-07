//
//  main.cpp
//  sada_fd_ms
//
//  Created by denas on 4/6/17.
//  Copyright © 2017 denas. All rights reserved.
//

#include <iostream>
#include <string>
#include <tuple>

#include "utils.hpp"
#include "Bwt.hpp"
#include "stree_sada.hpp"

using namespace std;
using namespace fdms;


typedef fdms::bp_support_g<> t_bp_support;
typedef StreeSada<t_bp_support> Cst;
typedef typename Cst::node_type node_type;
typedef tuple<size_type, size_type, node_type> runs_rt;


string s, t;
Counter time_usage, space_usage;


size_type find_k_prim(size_type __k, size_type max__k, sdsl::bit_vector &__runs){
    while(++__k < max__k && __runs[__k] != 0)
        ;
    return __k;
}

Interval bstep_interval(const Bwt& bwt_, Interval cur_i, char c){
    int cc = bwt_.char2int[c];
    return std::make_pair(bwt_.C[cc] + bwt_.rank(cur_i.first, c),
                          bwt_.C[cc] + bwt_.rank(cur_i.second + 1, c) - 1);
}

sdsl::bit_vector bps(const string &s){
    sdsl::cst_sada<> st_of_s;
    sdsl::construct_im(st_of_s, s, 1);
    return st_of_s.bp;
}

runs_rt fill_runs_sada(const string &s, const string &t, sdsl::bit_vector &runs){
    cerr << "building indexes ";
    Bwt bwt(s);
    cerr << ". ";
    sdsl::bit_vector bp = bps(s);
    cerr << ". ";
    t_bp_support bp_supp(&bp);
    cerr << ". ";
    Cst st(bp_supp, bwt);
    cerr << "DONE"  << endl;

    size_type k = t.size();
    uint8_t c = t[k - 1];
    node_type v = st.wl(st.root(), c);

    Interval I(bwt.C[bwt.char2int[c]], bwt.C[bwt.char2int[c] + 1] - 1);
    while(--k > 0){
        c = t[k - 1];
        I = bstep_interval(bwt, I, c);

        if(I.first > I.second){
            runs[k] = 0;
            do{
                v = st.parent(v);
                I = bstep_interval(bwt, make_pair(st.lb(v), st.rb(v)), c);
            } while(I.first > I.second);
            // idx of last 0 in runs - 1 (within this block) and corresponding wl(node)
        } else {
            runs[k] = 1;
        }
        v = st.wl(v, c);
    }
    return make_tuple(0, 0, 0);
}


void build_ms_ohleb(const string &s, const string &t, sdsl::bit_vector &runs, sdsl::bit_vector &ms){
    cerr << "building indexes ";
    Bwt bwt(s);
    cerr << ". ";
    sdsl::bit_vector bp = bps(s);
    cerr << ". ";
    t_bp_support bp_supp(&bp);
    cerr << ". ";
    Cst st(bp_supp, bwt);
    cerr << "DONE"  << endl;

    size_type k = 0, h_star = k + 1, h = h_star, k_prim, ms_idx = 0;
    char_type c = t[k];
    Interval I(bwt.C[bwt.char2int[c]], bwt.C[bwt.char2int[c] + 1] - 1);

    node_type v = st.wl(st.root(), c);
    while(k < t.size()){
        h = h_star;
        for(; I.first <= I.second && h_star < t.size(); ){
            c = t[h_star];
            I = bstep_interval(bwt, make_pair(st.lb(v), st.rb(v)), c);
            if(I.first <= I.second){
                v = st.wl(v, c);
                h_star += 1;
            }
        }
        ms_idx += (h_star - h + 1);
        if(h_star - h + 1)
            ms[ms_idx++] = 1;

        if(h_star < t.size()){
            do{
                v = st.parent(v);
                I = bstep_interval(bwt, make_pair(st.lb(v), st.rb(v)), t[h_star]);
            } while(I.first > I.second);
            h_star += 1;
        }
        v = st.wl(v, c);

        k_prim = find_k_prim(k, t.size(), runs);
        for(size_type i = k + 1; i <= k_prim - 1; i++)
            ms[ms_idx++] = 1;
        k = k_prim;
    }
}


void comp(InputSpec& T, InputSpec& S_fwd, const string& out_path, InputFlags& flags){
    cerr << "loading input ";
    t = T.load_s();
    cerr << ". ";
    s = S_fwd.load_s();
    cerr << ". ";
    cerr << "DONE" << endl;

    sdsl::bit_vector runs(t.size());
    sdsl::bit_vector ms(2 * t.size()); // the ms vector for each thread

    cerr << "build ms ";
    fill_runs_sada(s, t, runs);
    cerr << "DONE" << endl;

    cerr << "reverse s ";
    reverse_in_place(s);
    cerr << "DONE" << endl;

    cerr << "build ms ";
    build_ms_ohleb(s, t, runs, ms);
    cerr << "DONE" << endl;

    if(flags.answer){
        dump_ms(ms);
        cout << endl;
    }
}

int main(int argc, char **argv) {
    OptParser input(argc, argv);
    if(argc == 1){
        const string base_dir = {"/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/datasets/testing/"};
        InputFlags flags(false, // lazy_wl
                         false, // sada cst
                         false, // space
                         false, // time
                         true,  // ans
                         false, // verbose
                         10,    // nr. progress messages for runs construction
                         10,    // nr. progress messages for ms construction
                         false, // load CST
                         1      // nthreads
                         );
        InputSpec tspec(base_dir + "abcde200_64t.txt");
        InputSpec sfwd_spec(base_dir + "abcde200_64s.txt");
        const string out_path = "0";
        comp(tspec, sfwd_spec, out_path, flags);
    } else {
        InputFlags flags(input);
        InputSpec tspec(input.getCmdOption("-t_path"));
        InputSpec sfwd_spec(input.getCmdOption("-s_path"));
        const string out_path = input.getCmdOption("-out_path");
        comp(tspec, sfwd_spec, out_path, flags);
    }
}
