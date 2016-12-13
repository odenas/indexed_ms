/*
 * fabio_djamal_ms.cpp
 *
 *  Created on: Oct 13, 2016
 *      Author: denas
 */

#include <iostream>
#include <string>

#include "basic.hpp"
#include "fd_ms.hpp"
#include "stree_sada.hpp"
#include "stree_sct3.hpp"


#define STREE_SADA

using namespace std;
using namespace fdms;
using timer = std::chrono::high_resolution_clock;


void build_ms_sada(const string& prefix, string& t, InputSpec& S_rev, bvector& runs, bvector& ms, const InputFlags& flags){
    typedef fdms::bp_support_g<> t_bp_support;

    string s = S_rev.load_s();
    bvector bp = S_rev.load_bps();
    t_bp_support bp_supp(&bp);
    Bwt bwt(s);
    StreeSada<t_bp_support> st(bp_supp, bwt);
    build_ms_from_st_and_bwt(st, bwt, t, prefix, runs, ms, flags.space_usage, flags.verbose);
    if(flags.space_usage){
        cout << prefix << ", 2, build_ms, space, byte, s, " << s.size() << endl;
        cout << prefix << ", 2, build_ms, space, byte, s_bwt_wtree, " << bwt.size_in_bytes__wtree << endl;
        cout << prefix << ", 2, build_ms, space, byte, s_bwt_bwt, " << bwt.size_in_bytes__bwt << endl;
        cout << prefix << ", 2, build_ms, space, byte, s_bwt_alp, " << bwt.size_in_bytes__C + bwt.size_in_bytes__char2int + bwt.size_in_bytes__Sigma << endl;
        cout << prefix << ", 2, build_ms, space, byte, s_stree_bp, " << sdsl::size_in_bytes(bp) << endl;
        cout << prefix << ", 2, build_ms, space, byte, s_stree_bpsupp, " << bp_supp.size_in_bytes_total << endl;
        cout << prefix << ", 2, build_ms, space, byte, s_stree_bpsupp_select, " << st.size_in_bytes__select << endl;
        cout << prefix << ", 2, build_ms, space, byte, s_stree_bpsupp_rank, " << st.size_in_bytes__rank << endl;
    }
}

void build_ms_ohleb(const string& prefix, string& t, InputSpec& S_rev, bvector& runs, bvector& ms, const InputFlags& flags){
    StreeOhleb<> st;
    sdsl::construct(st, S_rev.s_fname, 1);
    string s = S_rev.load_s();
    Bwt bwt(s);

    build_ms_from_st_and_bwt(st, bwt, t, prefix, runs, ms, flags.space_usage, flags.verbose);
    if(flags.space_usage){
        cout << prefix << ", 2, build_ms, space, byte, s, " << s.size() << endl;
        cout << prefix << ", 2, build_ms, space, byte, s_bwt_wtree, " << bwt.size_in_bytes__wtree << endl;
        cout << prefix << ", 2, build_ms, space, byte, s_bwt_bwt, " << bwt.size_in_bytes__bwt << endl;
        cout << prefix << ", 2, build_ms, space, byte, s_bwt_alp, " << bwt.size_in_bytes__C + bwt.size_in_bytes__char2int + bwt.size_in_bytes__Sigma << endl;
        cout << prefix << ", 2, build_runs_sada, space, byte, s_stree_bpsupp_select, NA" << endl;
    }
}

void build_runs_sada(const string& prefix, string& t, InputSpec& S_fwd, bvector& runs, const InputFlags& flags){
    string s = S_fwd.load_s();
    Bwt bwt(s);
    bvector bp = S_fwd.load_bps();
    bp_support_g<> bp_supp(&bp);
    StreeSada<bp_support_g<>> st(bp_supp, bwt);

    build_runs_from_st_and_bwt(st, bwt, t, runs, flags.verbose);
    if(flags.space_usage){
        cout << prefix << ", 2, build_runs_sada, space, byte, s, " << s.size() << endl;
        cout << prefix << ", 2, build_runs_sada, space, byte, s_bwt_wtree, " << bwt.size_in_bytes__wtree << endl;
        cout << prefix << ", 2, build_runs_sada, space, byte, s_bwt_bwt, " << bwt.size_in_bytes__bwt << endl;
        cout << prefix << ", 2, build_runs_sada, space, byte, s_bwt_alp, " << bwt.size_in_bytes__C + bwt.size_in_bytes__char2int + bwt.size_in_bytes__Sigma << endl;
        cout << prefix << ", 2, build_runs_sada, space, byte, s_stree_bp, " << sdsl::size_in_bytes(bp) << endl;
        cout << prefix << ", 2, build_runs_sada, space, byte, s_stree_bpsupp, " << bp_supp.size_in_bytes_total << endl;
        cout << prefix << ", 2, build_runs_sada, space, byte, s_stree_bpsupp_select, " << st.size_in_bytes__select << endl;
        cout << prefix << ", 2, build_runs_sada, space, byte, s_stree_bpsupp_rank, " << st.size_in_bytes__rank << endl;
    }
}

void build_runs_ohleb(const string& prefix, string& t, InputSpec& S_fwd, bvector& runs, const InputFlags& flags){
    StreeOhleb<> st;
    sdsl::construct(st, S_fwd.s_fname, 1);
    string s = S_fwd.load_s();
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

template<class t_bp_support>
void comp(const string& prefix, InputSpec& T, InputSpec& S_fwd, InputSpec& S_rev, const InputFlags& flags){
    string t = T.load_s();
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
    //build_runs_sada(prefix, t, S_fwd, runs, flags);
    build_runs_ohleb(prefix, t, S_fwd, runs, flags);
    auto runs_stop = timer::now();

    auto ms_start = timer::now();
    build_ms_sada(prefix, t, S_rev, runs, ms, flags);
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
        string prefix {"ab200_2"};
        InputFlags flags(false, false, false, true, true);
        InputSpec tspec(base_dir + prefix + "t.txt", string(""));
        InputSpec sfwd_spec(base_dir + prefix + "s_fwd.txt", base_dir + prefix + "s_fwd_bp.txt");
        InputSpec srev_spec(base_dir + prefix + "s_rev.txt", base_dir + prefix + "s_rev_bp.txt");
        comp<fdms::bp_support_g<>>(prefix, tspec, sfwd_spec, srev_spec, flags);
    } else {
        const string& base_dir = input.getCmdOption("-d");
        const string& prefix = input.getCmdOption("-p");
        InputFlags flags(input);
        InputSpec tspec(base_dir + "/" + prefix + "t.txt", string(""));
        InputSpec sfwd_spec(base_dir + "/" + prefix + "s_fwd.txt", base_dir + "/" + prefix + "s_fwd_bp.txt");
        InputSpec srev_spec(base_dir + "/" + prefix + "s_rev.txt", base_dir + "/" + prefix + "s_rev_bp.txt");
        comp<fdms::bp_support_g<>>(prefix, tspec, sfwd_spec, srev_spec, flags);
    }
    return 0;
}

