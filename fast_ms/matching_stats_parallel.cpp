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

//#define VERBOSE

#include "fd_ms/input_spec.hpp"
#include "fd_ms/opt_parser.hpp"
#include "fd_ms/stree_sct3.hpp"
#include "fd_ms/maxrep_vector.hpp"
#include "fd_ms/runs_and_ms_algorithms.hpp"
#include "fd_ms/slices.hpp"

using namespace std;
using namespace fdms;


string t, s;
StreeOhleb<> st;
MsVectors<StreeOhleb<>, sdsl::bit_vector> ms_vec;
Maxrep<StreeOhleb<>, sdsl::bit_vector> maxrep;
Counter time_usage;


class InputFlags{
private:
    void check() const {
        if(use_maxrep_rc && use_maxrep_vanilla){
            cerr << "use_maxrep_rc and use_maxrep_vanilla cannot be active at the same time" << endl;
            exit(1);
        }
        if (use_maxrep() && lazy){
            cerr << "lazy and use_maxrep_xx cannot be active at the same time" << endl;
            cerr << "use_maxrep_xx goes with double rank and fail" << endl;
            exit(1);
        }
        if (use_maxrep() && !double_rank){
            cerr << "single_rank and use_maxrep_xx cannot be active at the same time" << endl;
            cerr << "use_maxrep_xx goes with double rank and fail" << endl;
            exit(1);
        }
        if (rank_fail && !double_rank){
            cerr << "single_rank and rank_fail cannot be active at the same time" << endl;
            exit(1);
        }

        if(nthreads == 0){
            cerr << "nr. of threads (parallelism) should be a positive number (got " << nthreads << ")" << endl;
            exit(1);
        }
    }

public:
    bool double_rank, lazy, rank_fail, use_maxrep_vanilla, use_maxrep_rc, lca_parents;
    bool time_usage, answer;
    bool load_stree, load_maxrep;
    size_t nthreads;
    
    InputFlags(){}
    
    InputFlags(const InputFlags& f) :
        double_rank{f.double_rank},
		lazy{f.lazy},
		rank_fail{f.rank_fail},
		use_maxrep_vanilla{f.use_maxrep_vanilla}, use_maxrep_rc{f.use_maxrep_rc},
		lca_parents{f.lca_parents},
		time_usage{f.time_usage},
		answer{f.answer},
		load_stree{f.load_stree},
		load_maxrep{f.load_maxrep},
		nthreads{f.nthreads}{}
    
    InputFlags(bool double_rank, bool lazy_wl, bool use_rank_fail, bool use_maxrep_vanilla, bool use_maxrep_rc, bool lca_parents,
               bool time_, bool ans,
               bool load_stree, bool load_maxrep, size_t nthreads) :
        double_rank{double_rank},
		lazy{lazy_wl},
		rank_fail{use_rank_fail},
		use_maxrep_vanilla{use_maxrep_vanilla}, use_maxrep_rc{use_maxrep_rc},
		lca_parents{lca_parents},
		time_usage {time_},
		answer {ans},
		load_stree{load_stree},
		load_maxrep{load_maxrep},
		nthreads{nthreads}
	{ check(); }
    
    InputFlags (OptParser input) :
		double_rank {input.getCmdOption("-double_rank") == "1"},  // use double rank
		lazy {input.getCmdOption("-lazy_wl") == "1"},             // lazy winer links
		rank_fail {input.getCmdOption("-rank_fail") == "1"},      // use the rank-and-fail strategy
		use_maxrep_rc {input.getCmdOption("-use_maxrep_rc") == "1"},    // use the maxrep vector with rank_and_check
		use_maxrep_vanilla {input.getCmdOption("-use_maxrep_vanilla") == "1"},    // use the maxrep vector the vanilla way
		lca_parents {input.getCmdOption("-lca_parents") == "1"},  // use lca insted of conscutive parent calls
		time_usage {input.getCmdOption("-time_usage") == "1"},    // time usage
		answer {input.getCmdOption("-answer") == "1"},            // answer
		load_stree{input.getCmdOption("-load_cst") == "1"},       // load CST of S and S'
		load_maxrep{input.getCmdOption("-load_maxrep") == "1"},   // load MAXREP of S'
		nthreads{static_cast<size_t>(std::stoi(input.getCmdOption("-nthreads")))}
    {
        nthreads = (nthreads <= 0 ? 1 : nthreads);
        check();
    }

    bool use_maxrep() const { return (use_maxrep_rc || use_maxrep_vanilla); }

	wl_method_t1 get_wl_method() {
        if(double_rank){
    	    if(lazy)
	            return (rank_fail ?
	        		    &StreeOhleb<>::lazy_double_rank_fail_wl :
	        		    &StreeOhleb<>::lazy_double_rank_nofail_wl);
            return (rank_fail ?
                    &StreeOhleb<>::double_rank_fail_wl :
                    &StreeOhleb<>::double_rank_nofail_wl);
        } else {
            return (lazy ?
        		    &StreeOhleb<>::lazy_single_rank_wl :
        		    &StreeOhleb<>::single_rank_wl);
        }
	}

	wl_method_t2 get_mrep_wl_method() {
        return (use_maxrep_rc ?
        		&StreeOhleb<>::double_rank_fail_wl_mrep:
        		&StreeOhleb<>::double_rank_fail_wl_mrep);
	}
};


runs_rt fill_runs_slice_thread(const size_type thread_id, const Interval slice, node_type v, InputFlags flags){
    // runs does not support laziness
    flags.lazy = false;
    return fill_runs_slice(t, st, flags.get_wl_method(), flags.lca_parents, ms_vec, v, slice);
}


void build_runs_ohleb(const InputFlags& flags, const InputSpec &s_fwd){
    cerr << "building RUNS ... " << endl;

    /* build the CST */
    time_usage["runs_cst"]  = load_or_build(st, s, s_fwd.fwd_cst_fname, flags.load_stree);
    cerr << "DONE (" << time_usage["runs_cst"] / 1000 << " seconds, " << st.size() << " nodes)" << endl;

    /* compute RUNS */
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
    if(flags.use_maxrep_rc || flags.use_maxrep_vanilla)
        return fill_ms_slice_maxrep(t, st, flags.get_mrep_wl_method(), ms_vec, maxrep, thread_id, slice);
    return fill_ms_slice(t, st, flags.get_wl_method(), flags.lca_parents, ms_vec, thread_id, slice);
}

void build_ms_ohleb(const InputFlags& flags, InputSpec &s_fwd){
    cerr << "building MS ... " << endl;

    /* build the CST */
    time_usage["ms_cst"] = load_or_build(st, s, s_fwd.rev_cst_fname, flags.load_stree);
    cerr << "DONE (" << time_usage["ms_cst"] / 1000 << " seconds, " << st.size() << " nodes)" << endl;

    /* build the maxrep vector */
    if(flags.use_maxrep()){
        time_usage["ms_maxrep"] = Maxrep<StreeOhleb<>, sdsl::bit_vector>::load_or_build(maxrep, st, s_fwd.rev_maxrep_fname, flags.load_maxrep);
        cerr << "DONE (" << time_usage["ms_maxrep"] / 1000 << " seconds)" << endl;
    }

    /* build MS */
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

    if(flags.time_usage){
        cerr << "dumping reports" << endl;
        cout << "len_s,len_t,item,value" << endl;
        if(flags.time_usage){
            for(auto item : time_usage)
                cout << s.size() << "," << t.size() << "," << item.first << "," << item.second << endl;
        }
    }

    if(flags.answer){
        if(out_path == "0"){
            ms_vec.show_MS(cout);
            cout << endl;
        }
    }
}


int main(int argc, char **argv){
    OptParser input(argc, argv);
    InputSpec sfwd_spec, tspec;
    InputFlags flags;
    string out_path;

    if(argc == 1){
        const string base_dir = {"/home/brt/Documents/projects/matching_statistics/indexed_ms/tests/datasets/testing/"};
        tspec = InputSpec(base_dir + "rep_100000s_1000t.t");
        sfwd_spec = InputSpec(base_dir + "rep_100000s_1000t.s");
        out_path = "0";
        flags = InputFlags(true,  // use double rank
                           false, // lazy_wl
                           false, // rank-and-fail
                           false, // use maxrep vanilla
                           false, // use maxrep rank&check
                           false, // lca_parents
                           false, // time
                           true,  // ans
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

