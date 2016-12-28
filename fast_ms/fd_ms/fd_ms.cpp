/*
 * fabio_djamal_ms.cpp
 *
 *  Created on: Oct 13, 2016
 *      Author: denas
 */

#include <iostream>
#include <string>

#include <sdsl/suffix_trees.hpp>

#include "basic.hpp"
#include "fd_ms.hpp"
#include "stree_sada.hpp"
#include "stree_sct3.hpp"


//#define STREE_SADA

using namespace std;
using namespace fdms;
using timer = std::chrono::high_resolution_clock;



bvector construct_bp(string& _s){
    sdsl::cst_sada<> temp_st;
    sdsl::construct_im(temp_st, _s, 1);
    bvector _bp(temp_st.bp.size());
    for(int i = 0; i < _bp.size(); i++)
        _bp[i] = temp_st.bp[i];
    return _bp;
}

void reverse_in_place(string& s){
    size_type n = s.size();

    for(int i = 0; i < n / 2; i++){
        char c = s[i];
        s[i] = s[n - 1 - i];
        s[n - 1 - i] = c;
    }
}

monitor::size_dict build_ms_sada(const string& prefix, string& t, string& s_rev, bvector& runs, bvector& ms, const InputFlags& flags){
    typedef fdms::bp_support_g<> t_bp_support;
    monitor::size_dict space_usage;

    Bwt bwt(s_rev);
    bvector bp = construct_bp(s_rev);
    t_bp_support bp_supp(&bp);
    StreeSada<t_bp_support> st(bp_supp, bwt);

    size_type size_in_bytes_alg = build_ms_from_st_and_bwt(st, bwt, t, prefix, runs, ms, flags.lazy);
    if(flags.space_usage){
        space_usage["bwt_wtree"] = bwt.size_in_bytes__wtree;
        space_usage["bwt_bwt"] = bwt.size_in_bytes__bwt;
        space_usage["bwt_alp"] = bwt.size_in_bytes__C + bwt.size_in_bytes__char2int + bwt.size_in_bytes__Sigma;
        space_usage["stree_bp"] = sdsl::size_in_bytes(bp);
        space_usage["stree_bpsupp"] = bp_supp.size_in_bytes_total;
        space_usage["stree_bpsupp_select"] = st.size_in_bytes__select;
        space_usage["stree_bpsupp_rank"] = st.size_in_bytes__rank;
        space_usage["alg"]                 = size_in_bytes_alg;
    }
    return space_usage;
}

monitor::size_dict build_ms_ohleb(const string& prefix, string& t, string& s_rev, bvector& runs, bvector& ms, const InputFlags& flags){
    monitor::size_dict space_usage;
    StreeOhleb<> st;
    sdsl::construct_im(st, s_rev, 1);
    Bwt bwt(s_rev);

    size_type size_in_bytes_alg = build_ms_from_st_and_bwt(st, bwt, t, prefix, runs, ms, flags.lazy);
    if(flags.space_usage){
        space_usage["bwt_wtree"]           = bwt.size_in_bytes__wtree;
        space_usage["bwt_bwt"]             = bwt.size_in_bytes__bwt;
        space_usage["bwt_alp"]             = bwt.size_in_bytes__C + bwt.size_in_bytes__char2int + bwt.size_in_bytes__Sigma;
        space_usage["stree_csa"]           = sdsl::size_in_bytes(st.csa);

        space_usage["stree_bp"]            = sdsl::size_in_bytes(st.bp);
        space_usage["stree_bpsupp"]        = sdsl::size_in_bytes(st.bp_support);
        space_usage["alg"]                 = size_in_bytes_alg;
    }
    return space_usage;
}


monitor::size_dict build_runs_sada(const string& prefix, string& t, string& s, bvector& runs, const InputFlags& flags){
    typedef fdms::bp_support_g<> t_bp_support;
    monitor::size_dict space_usage;

    Bwt bwt(s);
    bvector bp = construct_bp(s);
    t_bp_support bp_supp(&bp);
    StreeSada<t_bp_support> st(bp_supp, bwt);

    size_type size_in_bytes_alg = build_runs_from_st_and_bwt(st, bwt, t, runs);
    if(flags.space_usage){
        space_usage["bwt_wtree"]           = bwt.size_in_bytes__wtree;
        space_usage["bwt_bwt"]             = bwt.size_in_bytes__bwt;
        space_usage["bwt_alp"]             = bwt.size_in_bytes__C + bwt.size_in_bytes__char2int + bwt.size_in_bytes__Sigma;
        space_usage["stree_bp"]            = sdsl::size_in_bytes(bp);
        space_usage["stree_bpsupp"]        = bp_supp.size_in_bytes_total;
        space_usage["stree_bpsupp_select"] = st.size_in_bytes__select;
        space_usage["stree_bpsupp_rank"]   = st.size_in_bytes__rank;
        space_usage["alg"]                 = size_in_bytes_alg;
    }
    return space_usage;
}

monitor::size_dict build_runs_ohleb(const string& prefix, string& t, string& s, bvector& runs, const InputFlags& flags){
    monitor::size_dict space_usage;
    StreeOhleb<> st;
    sdsl::construct_im(st, s, 1);
    Bwt bwt(s);
    size_type size_in_bytes_alg = build_runs_from_st_and_bwt(st, bwt, t, runs);
    if(flags.space_usage){
        space_usage["bwt_wtree"]           = bwt.size_in_bytes__wtree;
        space_usage["bwt_bwt"]             = bwt.size_in_bytes__bwt;
        space_usage["bwt_alp"]             = bwt.size_in_bytes__C + bwt.size_in_bytes__char2int + bwt.size_in_bytes__Sigma;
        space_usage["stree_csa"]           = sdsl::size_in_bytes(st.csa);

        space_usage["stree_bp"]            = sdsl::size_in_bytes(st.bp);
        space_usage["stree_bpsupp"]        = sdsl::size_in_bytes(st.bp_support);
        space_usage["alg"]                 = size_in_bytes_alg;
    }
    return space_usage;
}

void comp(const string& prefix, InputSpec& T, InputSpec& S_fwd, const InputFlags& flags){
    monitor::size_dict runs_space_usage, ms_space_usage, space_usage, time_usage;

    string t = T.load_s();
    string s = S_fwd.load_s();
    bvector runs(t.size());
    bvector ms(t.size() * 2);

    auto runs_start = timer::now();
    //runs_space_usage = build_runs_sada(prefix, t, s, runs, flags);
    runs_space_usage = build_runs_ohleb(prefix, t, s, runs, flags);
    auto runs_stop = timer::now();
    time_usage["runs"] = std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count();

    auto ms_start = timer::now();
    reverse_in_place(s);
    //ms_space_usage = build_ms_sada(prefix, t, s, runs, ms, flags);
    ms_space_usage = build_ms_ohleb(prefix, t, s, runs, ms, flags);
    auto ms_stop = timer::now();
    time_usage["ms"] = std::chrono::duration_cast<std::chrono::milliseconds>(ms_stop - ms_start).count();
    time_usage["total_time"] = time_usage["ms"] + time_usage["runs"];

    space_usage["s"] = s.size();
    space_usage["t"] = t.size();
    space_usage["runs"] = runs.size();
    space_usage["ms"] = ms.size();
    {
        unsigned long rss, vs;
        getmem(&rss, &vs);
        space_usage["resident_mem"] = rss;
        space_usage["virtual_mem"]  = vs;
    }
    for(auto item : ms_space_usage)
        space_usage[item.first] = max(ms_space_usage[item.first], runs_space_usage[item.first]);


    if(flags.space_or_time_usage){
        cout << "prefix,measuring,unit,item,value" << endl;
        for(auto item : space_usage)
            cout << prefix << ",space,byte," << item.first << "," << item.second << endl;
        for(auto item : time_usage)
            cout << prefix << ",time,ms," << item.first << "," << item.second << endl;
    }

    if(flags.answer){
        cout << prefix << " ";
        dump_ms(ms);
    }

}


int main(int argc, char **argv){
    InputParser input(argc, argv);
    if(argc == 1){
        const string base_dir = {"/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/input_data/"};
        string prefix {"ab200_8"};
        InputFlags flags(true, false, false, true, true);
        InputSpec tspec(base_dir + prefix + "t.txt");
        InputSpec sfwd_spec(base_dir + prefix + "s.txt");
        comp(prefix, tspec, sfwd_spec, flags);
    } else {
        const string& base_dir = input.getCmdOption("-d");
        const string& prefix = input.getCmdOption("-p");
        InputFlags flags(input);
        InputSpec tspec(base_dir + "/" + prefix + "t.txt");
        InputSpec sfwd_spec(base_dir + "/" + prefix + "s.txt");
        comp(prefix, tspec, sfwd_spec, flags);
    }
    return 0;
}

