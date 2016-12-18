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

void build_ms_sada(const string& prefix, string& t, string& s_rev, bvector& runs, bvector& ms, const InputFlags& flags){
    typedef fdms::bp_support_g<> t_bp_support;

    Bwt bwt(s_rev);
    bvector bp = construct_bp(s_rev);
    t_bp_support bp_supp(&bp);
    StreeSada<t_bp_support> st(bp_supp, bwt);

    build_ms_from_st_and_bwt(st, bwt, t, prefix, runs, ms, flags.space_usage, flags.verbose);
    if(flags.space_usage){
        cout << prefix << ", 2, build_ms, space, byte, s, " << s_rev.size() << endl;
        cout << prefix << ", 2, build_ms, space, byte, s_bwt_wtree, " << bwt.size_in_bytes__wtree << endl;
        cout << prefix << ", 2, build_ms, space, byte, s_bwt_bwt, " << bwt.size_in_bytes__bwt << endl;
        cout << prefix << ", 2, build_ms, space, byte, s_bwt_alp, " << bwt.size_in_bytes__C + bwt.size_in_bytes__char2int + bwt.size_in_bytes__Sigma << endl;
        cout << prefix << ", 2, build_ms, space, byte, s_stree_bp, " << sdsl::size_in_bytes(bp) << endl;
        cout << prefix << ", 2, build_ms, space, byte, s_stree_bpsupp, " << bp_supp.size_in_bytes_total << endl;
        cout << prefix << ", 2, build_ms, space, byte, s_stree_bpsupp_select, " << st.size_in_bytes__select << endl;
        cout << prefix << ", 2, build_ms, space, byte, s_stree_bpsupp_rank, " << st.size_in_bytes__rank << endl;
    }
}

void build_ms_ohleb(const string& prefix, string& t, string& s_rev, bvector& runs, bvector& ms, const InputFlags& flags){
    StreeOhleb<> st;
    sdsl::construct_im(st, s_rev, 1);
    Bwt bwt(s_rev);

    build_ms_from_st_and_bwt(st, bwt, t, prefix, runs, ms, flags.space_usage, flags.verbose);
    if(flags.space_usage){
        cout << prefix << ", 2, build_ms, space, byte, s, " << s_rev.size() << endl;
        cout << prefix << ", 2, build_ms, space, byte, s_bwt_wtree, " << bwt.size_in_bytes__wtree << endl;
        cout << prefix << ", 2, build_ms, space, byte, s_bwt_bwt, " << bwt.size_in_bytes__bwt << endl;
        cout << prefix << ", 2, build_ms, space, byte, s_bwt_alp, " << bwt.size_in_bytes__C + bwt.size_in_bytes__char2int + bwt.size_in_bytes__Sigma << endl;
        cout << prefix << ", 2, build_runs_sada, space, byte, s_stree_bpsupp_select, NA" << endl;
    }
}



void build_runs_sada(const string& prefix, string& t, string& s, bvector& runs, const InputFlags& flags){
    typedef fdms::bp_support_g<> t_bp_support;

    Bwt bwt(s);
    bvector bp = construct_bp(s);
    t_bp_support bp_supp(&bp);
    StreeSada<t_bp_support> st(bp_supp, bwt);

    build_runs_from_st_and_bwt(st, bwt, t, runs, flags.verbose);
    if(flags.space_usage){
        cout << prefix << ", 2, build_runs_sada, space, byte, s, " << s.size() << endl;
        cout << prefix << ", 2, build_runs_sada, space, byte, s_bwt_wtree, " << bwt.size_in_bytes__wtree << endl;
        cout << prefix << ", 2, build_runs_sada, space, byte, s_bwt_bwt, " << bwt.size_in_bytes__bwt << endl;
        cout << prefix << ", 2, build_runs_sada, space, byte, s_bwt_alp, " << bwt.size_in_bytes__C + bwt.size_in_bytes__char2int + bwt.size_in_bytes__Sigma << endl;
        //cout << prefix << ", 2, build_runs_sada, space, byte, s_stree_bp, " << sdsl::size_in_bytes(bp) << endl;
        cout << prefix << ", 2, build_runs_sada, space, byte, s_stree_bpsupp, " << bp_supp.size_in_bytes_total << endl;
        cout << prefix << ", 2, build_runs_sada, space, byte, s_stree_bpsupp_select, " << st.size_in_bytes__select << endl;
        cout << prefix << ", 2, build_runs_sada, space, byte, s_stree_bpsupp_rank, " << st.size_in_bytes__rank << endl;
    }
}

void build_runs_ohleb(const string& prefix, string& t, string& s, bvector& runs, const InputFlags& flags){
    StreeOhleb<> st;
    sdsl::construct_im(st, s, 1);
    Bwt bwt(s);
    build_runs_from_st_and_bwt(st, bwt, t, runs, flags.verbose);
    if(flags.space_usage){
        cout << prefix << ", 2, build_runs_sada, space, byte, s, " << s.size() << endl;
        cout << prefix << ", 2, build_runs_sada, space, byte, s_bwt_wtree, " << bwt.size_in_bytes__wtree << endl;
        cout << prefix << ", 2, build_runs_sada, space, byte, s_bwt_bwt, " << bwt.size_in_bytes__bwt << endl;
        cout << prefix << ", 2, build_runs_sada, space, byte, s_bwt_alp, " << bwt.size_in_bytes__C + bwt.size_in_bytes__char2int + bwt.size_in_bytes__Sigma << endl;
        cout << prefix << ", 2, build_runs_sada, space, byte, s_stree_bpsupp_select, NA" << endl;
    }
}

void comp(const string& prefix, InputSpec& T, InputSpec& S_fwd, const InputFlags& flags){
    string t = T.load_s();
    string s = S_fwd.load_s();
    bvector runs(t.size());
    bvector ms(t.size() * 2);

    if(flags.space_usage || flags.mach_space_usage || flags.time_usage)
        cout << "prefix, level, func, measuring, unit, item, value" << endl;

    if(flags.space_usage){
        cout << prefix << ", 1, comp, space, byte, t, " << t.size() << endl;
        cout << prefix << ", 1, comp, space, byte, runs, " << sdsl::size_in_bytes(runs) << endl;
        cout << prefix << ", 1, comp, space, byte, ms, " << sdsl::size_in_bytes(ms) << endl;
    }

    auto runs_start = timer::now();
    build_runs_sada(prefix, t, s, runs, flags);
    //build_runs_ohleb(prefix, t, s, runs, flags);
    auto runs_stop = timer::now();

    auto ms_start = timer::now();
    reverse_in_place(s);
    build_ms_sada(prefix, t, s, runs, ms, flags);
    //build_ms_ohleb(prefix, t, s, runs, ms, flags);
    auto ms_stop = timer::now();

    if(flags.time_usage){
        cout << prefix << ", 1, comp, time, milliseconds, runs, " << std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count() << endl;
        cout << prefix << ", 1, comp, time, milliseconds, ms, " << std::chrono::duration_cast<std::chrono::milliseconds>(ms_stop - ms_start).count() << endl;
    }

    if(flags.answer){
        cout << prefix << " ";
        dump_ms(ms);
    }

    if(flags.mach_space_usage){
        unsigned long rss, vs;
        getmem(&rss, &vs);
        cout << prefix << ", 1, comp, resident_memory, byte, total, " << rss << endl;
        cout << prefix << ", 1, comp, virtual_memory, byte, total, " << vs << endl;
    }
}


int main(int argc, char **argv){
    InputParser input(argc, argv);
    if(argc == 1){
        const string base_dir = {"/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/input_data/"};
        string prefix {"ab200_8"};
        InputFlags flags(false, false, false, true, true);
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

