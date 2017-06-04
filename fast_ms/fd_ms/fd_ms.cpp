/*
 * fabio_djamal_ms.cpp
 *
 *  Created on: Oct 13, 2016
 *      Author: denas
 */

#include <iostream>
#include <string>
#include <future>
#include <tuple>
#include <thread>

#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/select_support.hpp>

#include <sdsl/iterators.hpp>
#include <sdsl/lcp.hpp>
#include <sdsl/bp_support.hpp>
#include <sdsl/csa_wt.hpp> // for std initialization of StreeOhleb
#include <sdsl/cst_iterators.hpp>
#include <sdsl/sdsl_concepts.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/suffix_tree_helper.hpp>
#include <sdsl/suffix_tree_algorithm.hpp>
#include <sdsl/util.hpp>

#include "utils.hpp"
#include "fd_ms.hpp"


using namespace std;
using namespace fdms;


typedef typename StreeOhleb<>::node_type node_type;


string t, s;
StreeOhleb<> st;
sdsl::bit_vector runs(1);
sdsl::bit_vector maxrep(1);
vector<sdsl::bit_vector> mses(1); // the ms vector for each thread
vector<Interval> ms_sizes(1);

Counter space_usage, time_usage;

runs_rt fill_runs_slice_thread(const size_type thread_id, const Interval slice, node_type v, InputFlags flags){
    // runs does not support laziness
    flags.lazy = false;
    return fill_runs_slice(t, st, flags.get_wl_method(), flags.get_rank_method(), flags.get_parent_seq_method(),
                           runs, v, slice.first, slice.second);
}

void build_runs_ohleb(const InputFlags& flags, const InputSpec &s_fwd){
    cerr << "building RUNS ... " << endl;

    /* build the CST */
    time_usage["runs_cst"]  = load_st<StreeOhleb<>>(st, s, s_fwd.fwd_cst_fname, flags.load_stree);
    cerr << "DONE (" << time_usage["runs_cst"] / 1000 << " seconds, " << st.size() << " nodes)" << endl;

    /* compute RUNS */
    cerr << flags.runs_strategy_string(1, true) << endl;

    auto runs_start = timer::now();
    std::vector<Interval> slices = slice_input(t.size(), flags.nthreads);
    std::vector<std::future<runs_rt>> results(flags.nthreads);
    for(size_type i=0; i<flags.nthreads; i++){
        cerr << " ** launching runs computation over : [" << slices[i].first << " .. " << slices[i].second << ")" << endl;
        node_type v = st.double_rank_nofail_wl(st.root(), t[slices[i].second - 1]); // stree node
        //fill_runs_slice(i, slices[i].first, slices[i].second);
        results[i] = std::async(std::launch::async, fill_runs_slice_thread, i, slices[i], v, flags);
	}
    vector<runs_rt> runs_results(flags.nthreads);
    for(size_type i=0; i<flags.nthreads; i++){
        runs_results[i] = results[i].get();
        //cerr << " *** [" << get<0>(runs_results[i]) << " .. " << get<1>(runs_results[i]) << ")" << endl;
    }

    results = std::vector<std::future<runs_rt>>(flags.nthreads);
    cerr << " ** merging over " << flags.nthreads - 1 << " threads ... " << endl;
    for(int i = (int) flags.nthreads - 1; i > 0; i--){
        cerr << " *** launching runs merge of slices " << i << " and " << i - 1 << " ... " << endl;
        results[i] = std::async(std::launch::async, fill_runs_slice_thread,
                                (size_type)i,
                                make_pair(i == 0 ? 0 : get<0>(runs_results[i - 1]), get<1>(runs_results[i])),
                                get<2>(runs_results[i]), flags);
                                //flags.rank_fail, flags.lca_parents);
    }
    for(int i = (int) flags.nthreads - 1; i > 0; i--)
        results[i].get();

    auto runs_stop = timer::now();
    time_usage["runs_bvector"]  = std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count();
    cerr << "DONE (" << time_usage["runs_bvector"] / 1000 << " seconds)" << endl;

    //cout << endl;
    //for(size_type i = 0; i < runs.size(); i++)
    //    cout << runs[i] << " ";
    //cout << endl;
}

Interval fill_ms_slice_thread(const size_type thread_id, const Interval slice, InputFlags flags){
    if(flags.use_maxrep){
        assert (!flags.lazy);
        return fill_ms_slice_nonlazy_fail(t, st, mses[thread_id], runs, maxrep, slice.first, slice.second);
    }


    return fill_ms_slice_lazy(t, st,
                              flags.get_wl_method(), flags.get_rank_method(), flags.get_parent_seq_method(),
                              mses[thread_id], runs, slice.first, slice.second);
}

void build_ms_ohleb(const InputFlags& flags, InputSpec &s_fwd){
    cerr << "building MS ... " << endl;

    /* build the CST */
    time_usage["ms_cst"] = load_st<StreeOhleb<>>(st, s, s_fwd.rev_cst_fname, flags.load_stree);
    cerr << "DONE (" << time_usage["ms_cst"] / 1000 << " seconds, " << st.size() << " nodes)" << endl;

    /* build the maxrep vector */
    if(flags.use_maxrep){
        time_usage["ms_maxrep"] = load_maxrep(maxrep, st, s, s_fwd.rev_maxrep_fname, flags.load_maxrep);
        cerr << "DONE (" << time_usage["ms_maxrep"] / 1000 << " seconds)" << endl;
    }

    /* build MS */
    cerr << flags.ms_strategy_string(1, true) << endl;
    auto runs_start = timer::now();
    std::vector<Interval> slices = slice_input(t.size(), flags.nthreads);
    std::vector<std::future<Interval>> results(flags.nthreads);
    for(size_type i=0; i<flags.nthreads; i++){
        cerr << " ** launching ms computation over : [" << slices[i].first << " .. " << slices[i].second << ")" << endl;
        //fill_ms_slice(i, slices[i].first, slices[i].second);
        results[i] = std::async(std::launch::async, fill_ms_slice_thread, i, slices[i], flags);
    }
    for(size_type i=0; i<flags.nthreads; i++){
        pair<size_type, size_type> rr = results[i].get();
        space_usage["ms_bvector_allocated" + std::to_string(i)] = rr.first;
        space_usage["ms_bvector_used" + std::to_string(i)]   = rr.second;
    }
    auto runs_stop = timer::now();
    time_usage["ms_bvector"] = std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count();

    size_type total_ms_length = 0;
    for(size_type i=0; i<flags.nthreads; i++)
        total_ms_length += mses[i].size();
    cerr << " * total ms length : " << total_ms_length << " (with |t| = " << t.size() << ")" << endl;
    cerr << "DONE (" << time_usage["ms_bvector"] / 1000 << " seconds)" << endl;
}

void comp(InputSpec& T, InputSpec& S_fwd, const string& out_path, InputFlags& flags){
    auto comp_start = timer::now();
    cerr << "loading input ";
    auto start = timer::now();
    t = T.load_s();
    cerr << ". ";
    s = S_fwd.load_s();
    cerr << ". ";
    cerr << "|s| = " << s.size() << ", |t| = " << t.size() << ". ";
    time_usage["loadstr"] = std::chrono::duration_cast<std::chrono::milliseconds>(timer::now() - start).count();
    cerr << "DONE (" << time_usage["loadstr"] / 1000 << " seconds)" << endl;
    space_usage["s"] = s.size();
    space_usage["t"] = t.size();


    /* prepare global data structures */
    // runs
    runs.resize(t.size());
    sdsl::util::set_to_value(runs, 0);
	// ms
    mses.resize(flags.nthreads);
    ms_sizes.resize(flags.nthreads);
    for(int i=0; i<flags.nthreads; i++){
        mses[i].resize(t.size() / flags.nthreads);
        sdsl::util::set_to_value(mses[i], 0);
    }
    // maxrep
    maxrep.resize(s.size() + 1); sdsl::util::set_to_value(maxrep, 0);
    // space usage
    space_usage["runs_bvector"] = runs.size();
    space_usage["ms_bvector"]   = (2 * t.size()) * flags.nthreads;
    space_usage["maxrep"] = maxrep.size();

    if(flags.ms_progress > t.size())
        flags.ms_progress = t.size() - 1;

    build_runs_ohleb(flags, S_fwd);

    start = timer::now();
    cerr << " * reversing string s of length " << s.size() << " ";
    reverse_in_place(s);
    time_usage["reverse_str"] = std::chrono::duration_cast<std::chrono::milliseconds>(timer::now() - start).count();
    cerr << "DONE (" << time_usage["reverse_str"] / 1000 << " seconds)" << endl;

    build_ms_ohleb(flags, S_fwd);
    time_usage["total_time"] = std::chrono::duration_cast<std::chrono::milliseconds>(timer::now() - comp_start).count();

    if(flags.space_usage || flags.time_usage){
        cerr << "dumping reports" << endl;
        cout << "len_s,len_t,measuring,item,value" << endl;
        if(flags.space_usage){
            for(auto item: space_usage)
                cout << s.size() << "," << t.size() << ",space," << item.first << "," << item.second << endl;
        }
        if(flags.time_usage){
            for(auto item : time_usage)
                cout << s.size() << "," << t.size() << ",time," << item.first << "," << item.second << endl;
        }
    }

    if(flags.answer){
        if(out_path == "0")
            for(size_type mses_idx=0; mses_idx < mses.size(); mses_idx++)
                dump_ms(mses[mses_idx]);
            cout << endl;
    }
}


int main(int argc, char **argv){
    OptParser input(argc, argv);
    if(argc == 1){
        const string base_dir = {"/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/"};
        InputFlags flags(false, // lazy_wl
                         false,  // rank-and-fail
                         true,  // use maxrep
                         false,  // lca_parents
                         false, // space
                         false, // time
                         true,  // ans
                         false, // verbose
                         10,    // nr. progress messages for runs construction
                         10,    // nr. progress messages for ms construction
                         false, // load CST
                         false, // load MAXREP
                         1      // nthreads
                         );
        InputSpec tspec(base_dir + "rnd_20_10.t");
        InputSpec sfwd_spec(base_dir + "rnd_20_10.s");
        const string out_path = "0";
        comp(tspec, sfwd_spec, out_path, flags);
    } else {
        InputFlags flags(input);
        InputSpec tspec(input.getCmdOption("-t_path"));
        InputSpec sfwd_spec(input.getCmdOption("-s_path"));
        const string out_path = input.getCmdOption("-out_path");
        comp(tspec, sfwd_spec, out_path, flags);
    }

    return 0;
}

