/*
 * fabio_djamal_ms.cpp
 *
 *  Created on: Oct 13, 2016
 *      Author: denas
 */


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <future>
#include <thread>


#include "input_spec.hpp"
#include "opt_parser.hpp"
#include "cmd_utils.hpp"
#include "stree_sct3.hpp"
#include "maxrep_vector.hpp"
#include "utils.hpp"
#include "runs_and_ms_algorithms.hpp"
#include "runs_ms.hpp"
#include "slices.hpp"

using namespace std;
using namespace fdms;


string t, s;
StreeOhleb<> st;
MsVectors<StreeOhleb<>, sdsl::bit_vector> ms_vec;
Maxrep<StreeOhleb<>, sdsl::bit_vector> maxrep;
Counter time_usage;


runs_rt fill_runs_slice_thread(const size_type thread_id, const Interval slice, node_type v, InputFlags flags){
    // runs does not support laziness
    flags.lazy = false;
    return fill_runs_slice(t, st, get_wl_method(flags), get_rank_method(flags), get_parent_seq_method(flags),
                           ms_vec, v, slice);
}

void build_runs_ohleb(const InputFlags& flags, const InputSpec &s_fwd){
    cerr << "building RUNS ... " << endl;

    /* build the CST */
    time_usage["runs_cst"]  = load_or_build(st, s, s_fwd.fwd_cst_fname, flags.load_stree);
    cerr << "DONE (" << time_usage["runs_cst"] / 1000 << " seconds, " << st.size() << " nodes)" << endl;

    /* compute RUNS */
    cerr << flags.runs_strategy_string(1, true) << endl;

    auto runs_start = timer::now();
    std::vector<std::future<runs_rt>> results(flags.nthreads);
    Slices<size_type> slices(t.size(), flags.nthreads);

    for(size_type i=0; i<flags.nthreads; i++){
        cerr << " ** launching runs computation over : " << slices.repr(i) << endl;
        node_type v = st.double_rank_nofail_wl(st.root(), t[slices[i].second - 1]); // stree node
        //fill_runs_slice_thread(i, slices[i], v, flags);
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
}

Interval fill_ms_slice_thread(const size_type thread_id, const Interval slice, InputFlags flags){
    if(!flags.use_maxrep)
        maxrep.set_to_one(st.size() + 1);

    return fill_ms_slice_maxrep(t, st,
                                get_rank_method(flags), get_parent_seq_method(flags),
                                ms_vec, maxrep, thread_id, slice);
}

void build_ms_ohleb(const InputFlags& flags, InputSpec &s_fwd){
    cerr << "building MS ... " << endl;

    /* build the CST */
    time_usage["ms_cst"] = load_or_build(st, s, s_fwd.rev_cst_fname, flags.load_stree);
    cerr << "DONE (" << time_usage["ms_cst"] / 1000 << " seconds, " << st.size() << " nodes)" << endl;

    /* build the maxrep vector */
    if(flags.use_maxrep){
        time_usage["ms_maxrep"] = Maxrep<StreeOhleb<>, sdsl::bit_vector>::load_or_build(maxrep, st, s_fwd.rev_maxrep_fname, flags.load_maxrep);
        cerr << "DONE (" << time_usage["ms_maxrep"] / 1000 << " seconds)" << endl;
    }

    /* build MS */
    cerr << flags.ms_strategy_string(1, true) << endl;
    auto runs_start = timer::now();
    Slices<size_type> slices(t.size(), flags.nthreads);
    std::vector<std::future<Interval>> results(flags.nthreads);
    for(size_type i=0; i<flags.nthreads; i++){
        cerr << " ** launching ms computation over : " << slices.repr(i) << endl;
        //fill_ms_slice_thread(i, slices[i], flags);
        results[i] = std::async(std::launch::async, fill_ms_slice_thread, i, slices[i], flags);
    }
    for(size_type i=0; i<flags.nthreads; i++){
        results[i].get();
    }
    auto runs_stop = timer::now();
    time_usage["ms_bvector"] = std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count();

    cerr << " * total ms length : " << ms_vec.ms_size()  << " (with |t| = " << t.size() << ")" << endl;
    cerr << "DONE (" << time_usage["ms_bvector"] / 1000 << " seconds)" << endl;
}

void comp(const InputSpec& tspec, InputSpec& S_fwd, const string& out_path, InputFlags& flags){
    auto comp_start = timer::now();
    cerr << "loading input ";
    auto start = timer::now();
    t = tspec.load_s();
    cerr << ". ";
    s = S_fwd.load_s();
    cerr << ". ";
    cerr << "|s| = " << s.size() << ", |t| = " << t.size() << ". ";
    time_usage["loadstr"] = std::chrono::duration_cast<std::chrono::milliseconds>(timer::now() - start).count();
    cerr << "DONE (" << time_usage["loadstr"] / 1000 << " seconds)" << endl;

    /* prepare global data structures */
    ms_vec = MsVectors<StreeOhleb<>, sdsl::bit_vector>(t.size(), flags.nthreads);

    build_runs_ohleb(flags, S_fwd);

    start = timer::now();
    cerr << " * reversing string s of length " << s.size() << " ";
    InputSpec::reverse_in_place(s);
    time_usage["reverse_str"] = std::chrono::duration_cast<std::chrono::milliseconds>(timer::now() - start).count();
    cerr << "DONE (" << time_usage["reverse_str"] / 1000 << " seconds)" << endl;

    build_ms_ohleb(flags, S_fwd);
    time_usage["total_time"] = std::chrono::duration_cast<std::chrono::milliseconds>(timer::now() - comp_start).count();

    if(flags.space_usage || flags.time_usage){
        cerr << "dumping reports" << endl;
        cout << "len_s,len_t,item,value" << endl;
        if(flags.time_usage){
            for(auto item : time_usage)
                cout << s.size() << "," << t.size() << "," << item.first << "," << item.second << endl;
        }
    }

    if(flags.answer){
        if(out_path == "0")
            for(size_type mses_idx=0; mses_idx < ms_vec.mses.size(); mses_idx++)
                dump_ms(ms_vec.mses[mses_idx]);
            cout << endl;
        
        size_type expected_array[200] = {4, 5, 4, 3, 4, 5, 4, 4, 4, 4, 4, 4, 3, 3, 5, 4, 5, 4, 3, 4, 4, 4, 3, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 5, 4, 3, 4, 3, 5, 4, 5, 4, 4, 4, 3, 6, 5, 4, 4, 5, 4, 4, 4, 4, 6, 5, 4, 4, 3, 4, 3, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 3, 5, 4, 3, 3, 3, 4, 5, 4, 4, 5, 5, 4, 5, 4, 5, 4, 3, 5, 4, 4, 5, 4, 3, 5, 4, 5, 4, 3, 4, 4, 5, 4, 6, 5, 5, 4, 4, 4, 4, 3, 4, 4, 6, 5, 4, 4, 4, 4, 5, 4, 4, 3, 4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 3, 5, 4, 5, 5, 4, 3, 4, 3, 4, 5, 6, 5, 5, 4, 4, 4, 3, 5, 4, 5, 4, 3, 3, 4, 5, 4, 4, 4, 4, 4, 4, 5, 4, 4, 3, 5, 4, 3, 6, 5, 4, 4, 4, 3, 5, 4, 3, 3, 4, 3, 4, 4, 4, 5, 4, 3, 4, 3, 2, 1};
        //for(size_type i = 0; i < 200; i++)
        //    cout << expected_array[i] << " ";
        //cout << endl;
/*
        size_type ii = 0, res = 0, exp_res = 0;
        for(size_type mses_idx=0; mses_idx < mses.size(); mses_idx++){
            size_type k = 0;
            for (size_type i = 0; i < mses[mses_idx].size(); i++){
                if(mses[mses_idx][i] == 1){
                    res = i - (2*k);
                    k += 1;
                    exp_res = expected_array[ii++];
                    cout << (res == exp_res ? "." : "*") << " ";
                }
            }

        }
        cout << endl;
*/
        
    }
}


int main(int argc, char **argv){
    OptParser input(argc, argv);
    InputSpec sfwd_spec, tspec;
    InputFlags flags;
    string out_path;

    if(argc == 1){
        const string base_dir = {"/Users/denas/projects/matching_statistics/indexed_ms/tests/datasets/testing/"};
        tspec = InputSpec(base_dir + "rnd_200_128.t");
        sfwd_spec = InputSpec(base_dir + "rnd_200_128.s");
        out_path = "0";
        flags = InputFlags(false, // lazy_wl
                           false,  // rank-and-fail
                           false,  // use maxrep
                           true,  // lca_parents
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
    } else {
        tspec = InputSpec(input.getCmdOption("-t_path"));
        sfwd_spec = InputSpec(input.getCmdOption("-s_path"));
        out_path = input.getCmdOption("-out_path");
        flags = InputFlags(input);
    }
    comp(tspec, sfwd_spec, out_path, flags);
    return 0;
}

