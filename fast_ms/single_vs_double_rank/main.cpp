//
//  main.cpp
//  single_vs_double_rank
//
//  Created by denas on 1/28/17.
//  Copyright © 2017 denas. All rights reserved.
//

#include <iostream>
#include "basic.hpp"
#include "fd_ms.hpp"
#include "stree_sct3.hpp"
#define NTRIALS 500000

using namespace std;
using namespace fdms;
using timer = std::chrono::high_resolution_clock;


monitor::size_dict time_wl_calls(string& s_rev, StreeOhleb<>& st, const size_type ntrials, monitor::size_dict& time_usage){
    typedef typename StreeOhleb<>::node_type node_type;

    size_type k = s_rev.size() - 1;
    node_type v = st.root();
    auto start_time = timer::now();
    for(size_type i=0; i<ntrials; i++)
        v = st.single_rank_wl(v, s_rev[k--]);

    time_usage["single_rank"] = std::chrono::duration_cast<std::chrono::milliseconds>(timer::now() - start_time).count();
    cerr << ntrials << " single rank wl calls took " << time_usage["single_rank"] << " ms" << endl;

    k = s_rev.size() - 1;
    v = st.root();
    start_time = timer::now();
    for(size_type i=0; i<ntrials; i++)
        v = st.wl(v, s_rev[k--]);

    time_usage["double_rank"] = std::chrono::duration_cast<std::chrono::milliseconds>(timer::now() - start_time).count();
    cerr << ntrials << " double rank wl calls took " << time_usage["double_rank"] << " ms" << endl;

    return time_usage;
}


int main(int argc, char **argv) {
    InputParser input(argc, argv);
    InputFlags flags(input);
    InputSpec S_fwd(input.getCmdOption("-s_path"));
    string s = S_fwd.load_s();
    size_type ntrials = (s.size() - 1 > NTRIALS ? NTRIALS : s.size() - 1);

    monitor::size_dict time_usage;
    StreeOhleb<> st;
    auto runs_start = timer::now();
    if(flags.load_stree){
        cerr << "loading the CST T(s) from " << S_fwd.s_fname + ".fwd.stree" << "... ";
        sdsl::load_from_file(st, S_fwd.s_fname + ".fwd.stree");
    } else {
        cerr << "building the CST T(s) of lentgth " << s.size() << "... ";
        sdsl::construct_im(st, s, 1);
    }
    auto runs_stop = timer::now();
    time_usage["dstruct"] = std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count();
    cerr << "DONE (" << time_usage["dstruct"] / 1000 << " seconds)" << endl;


    cout << "s_path,s,ntrials,nwlcalls,item,time_ms" << endl;
    for(size_type i = 1; i < 20; i++){
        monitor::size_dict tu = time_wl_calls(s, st, ntrials, time_usage);
        for(auto item : tu){
            cout << S_fwd.s_fname << "," << s.size() << "," << ntrials << "," << i << "," << item.first << "," << item.second << endl;
        }
    }
    return 0;
}
