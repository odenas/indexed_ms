#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <future>
#include <thread>

//#define VERBOSE

#include "fd_ms/opt_parser.hpp"
#include "fd_ms/input_spec.hpp"
#include "fd_ms/counter.hpp"
#include "fd_ms/query.hpp"
#include "fd_ms/slices.hpp"
#include "fd_ms/stree_sct3.hpp"
#include "fd_ms/maxrep_vector.hpp"
#include "fd_ms/runs_and_ms_algorithms.hpp"

using namespace std;
using namespace fdms;

typedef StreeOhleb<>               cst_t;
typedef typename cst_t::node_type  node_type;
typedef typename cst_t::size_type  size_type;
typedef typename cst_t::char_type  char_type;
typedef sdsl::bit_vector           bitvec_t;
typedef MsVectors<cst_t, bitvec_t> msvec_t;
typedef Maxrep<cst_t, bitvec_t>    maxrep_t;
typedef Counter<size_type>         counter_t;


cst_t st;
msvec_t ms_vec;
maxrep_t maxrep;

class InputFlags{
private:
    void check() const {
        if(use_maxrep_rc && use_maxrep_vanilla){
            cerr << "use_maxrep_rc and use_maxrep_vanilla cannot be active at the same time" << endl;
            exit(1);
        }
        if(use_maxrep() && !(double_rank && rank_fail)){
            cerr << "use_maxrep_xx goes with double rank and fail" << endl;
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
        if (use_maxrep() && !rank_fail){
            cerr << "no_fail and use_maxrep_xx cannot be active at the same time" << endl;
            cerr << "use_maxrep_xx goes with double rank and fail" << endl;
            exit(1);
        }
        if (rank_fail && !double_rank){
            cerr << "single_rank and rank_fail cannot be active at the same time" << endl;
            exit(1);
        }
        if (answer && avg){
            cerr << "answer and avg cannot be active at the same time" << endl;
            exit(1);
        }
    }

public:
    bool double_rank, lazy, rank_fail, use_maxrep_vanilla, use_maxrep_rc, lca_parents;
    bool time_usage, answer, avg;
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
        answer{f.answer}, avg{f.avg},
        load_stree{f.load_stree},
        load_maxrep{f.load_maxrep},
        nthreads{f.nthreads}{}
    
    InputFlags(bool double_rank, bool lazy_wl, bool use_rank_fail, bool use_maxrep_vanilla, bool use_maxrep_rc, bool lca_parents,
               bool time_, bool ans, bool avg,
               bool load_stree, bool load_maxrep, size_t nthreads) :
        double_rank{double_rank},
        lazy{lazy_wl},
        rank_fail{use_rank_fail},
        use_maxrep_vanilla{use_maxrep_vanilla}, use_maxrep_rc{use_maxrep_rc},
        lca_parents{lca_parents},
        time_usage {time_},
        answer {ans}, avg{avg}, 
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
        avg {input.getCmdOption("-avg") == "1"},                  // average matching statistics
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
                return (rank_fail ? &StreeOhleb<>::lazy_double_rank_fail_wl : &StreeOhleb<>::lazy_double_rank_nofail_wl);
            return (rank_fail ? &StreeOhleb<>::double_rank_fail_wl : &StreeOhleb<>::double_rank_nofail_wl);
        } else {
            return (lazy ? &StreeOhleb<>::lazy_single_rank_wl : &StreeOhleb<>::single_rank_wl);
        }
    }

    wl_method_t2 get_mrep_wl_method() {
        return (use_maxrep_rc ? &StreeOhleb<>::double_rank_fail_wl_mrep_vanilla: &StreeOhleb<>::double_rank_fail_wl_mrep_vanilla);
    }
};

StreeOhleb<>::size_type load_or_build(StreeOhleb<>& st, const InputSpec& ispec, const bool reverse, const bool load){
    string potential_stree_fname = (reverse ? ispec.rev_cst_fname : ispec.fwd_cst_fname);
    using timer = std::chrono::high_resolution_clock;
    auto start = timer::now();
    if(load){
        std::cerr << " * loading the CST from " << potential_stree_fname << " ";
        sdsl::load_from_file(st, potential_stree_fname);
    } else {
        std::cerr << " * loadding index string  from " << ispec.s_fname << " " << endl;
        string s = ispec.load_s(reverse);
        std::cerr << " * building the CST of length " << s.size() << " ";
        sdsl::construct_im(st, s, 1);
    }
    auto stop = timer::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
}

runs_rt fill_runs_slice_thread(const size_type thread_id, const Interval slice, node_type v, InputFlags flags, const string  t_fname){
    // runs does not support laziness
    flags.lazy = false;
    return fill_runs_slice(t_fname, st, flags.get_wl_method(), flags.lca_parents, ms_vec, v, slice);
}

vector<runs_rt> aa(const vector<runs_rt> v, const Slices<size_type> slices){
    vector<runs_rt> u;
    u.reserve(v.size());

    int i = 0, j = 1, n = v.size();
    while(j < n){
        runs_rt prev_slice = v[i], next_slice = v[j];
        while(get<0>(next_slice) == get<0>(slices.slices[j]) + 1){ // j-th slice had a full match
            j++;
        }
        next_slice = v[j - (j == n)];
        u.push_back(make_tuple(get<0>(prev_slice), get<1>(next_slice), get<2>(next_slice)));
        i = j;
        j += 1;
    }
    return u;

}


void build_runs_ohleb(const InputFlags& flags, const InputSpec &s_fwd, const InputSpec &tspec, counter_t &time_usage){
    cerr << "building RUNS ... " << endl;

    /* build the CST */
    time_usage.reg["runs_cst"]  = load_or_build(st, s_fwd, false, flags.load_stree);
    cerr << "DONE (" << time_usage.reg["runs_cst"] / 1000 << " seconds, " << st.size() << " leaves)" << endl;

    /* compute RUNS */
    auto runs_start = timer::now();
    std::vector<std::future<runs_rt>> results(flags.nthreads);
    Slices<size_type> slices(Query::query_length(tspec.s_fname), flags.nthreads);

    /* open connection to the query string */
    Query_rev t{tspec.s_fname, (size_t) 1e+5};
    for(size_type i=0; i<flags.nthreads; i++){
        node_type v = st.double_rank_nofail_wl(st.root(), t[slices[i].second - 1]); // stree node
        cerr << " ** launching runs computation over : " << slices.repr(i) << " ";
        cerr << "(" << v.i << ", " << v.j << ")" << endl;
        results[i] = std::async(std::launch::async, fill_runs_slice_thread, i, slices[i], v, flags, tspec.s_fname);
    }
    vector<runs_rt> runs_results(flags.nthreads);
    for(size_type i=0; i<flags.nthreads; i++){
        runs_results[i] = results[i].get();
        cerr << " *** [" << get<0>(runs_results[i]) << " .. " << get<1>(runs_results[i]) << ")" << endl;
    }
    vector<runs_rt> merge_idx = aa(runs_results, slices);
    //ms_vec.show_runs(cerr);

    cerr << " * merging over " << merge_idx.size() << " threads ... " << endl;
    for(int i = 0; i < (int) merge_idx.size(); i++){
        cerr << " ** ([" << get<0>(merge_idx[i]) << ", " << get<1>(merge_idx[i]) << "), " << "(" << get<2>(merge_idx[i]).i << ", " << get<2>(merge_idx[i]).j << ")) " << endl;
        results[i] = std::async(std::launch::async, fill_runs_slice_thread,
                                (size_type)i,
                                make_pair(get<0>(merge_idx[i]) - 1, get<1>(merge_idx[i])),
                                get<2>(merge_idx[i]), flags,
                                tspec.s_fname);
    }
    for(int i = 0; i < merge_idx.size(); i++)
        results[i].get();
    //ms_vec.show_runs(cerr);
    time_usage.register_now("runs_bvector", runs_start);
}

Interval fill_ms_slice_thread(const size_type thread_id, const Interval slice,
                           const node_type start_node, InputFlags flags, const InputSpec tspec){
    if(flags.use_maxrep_rc || flags.use_maxrep_vanilla)
        return fill_ms_slice_maxrep(tspec.s_fname, st, flags.get_mrep_wl_method(), ms_vec, maxrep, thread_id, start_node, slice);
    return fill_ms_slice(tspec.s_fname, st, flags.get_wl_method(), flags.lca_parents, ms_vec, thread_id, start_node, slice);
}

void build_ms_ohleb(const InputFlags& flags, InputSpec &s_fwd, const InputSpec &tspec, counter_t &time_usage){
    cerr << "building MS ... " << endl;

    /* build the CST */
    time_usage.reg["ms_cst"] = load_or_build(st, s_fwd, true, flags.load_stree);
    cerr << "DONE (" << time_usage.reg["ms_cst"] / 1000 << " seconds, " << st.size() << " nodes)" << endl;

    /* build the maxrep vector */
    if(flags.use_maxrep()){
        time_usage.reg["ms_maxrep"] = Maxrep<StreeOhleb<>, sdsl::bit_vector>::load_or_build(maxrep, st, s_fwd.rev_maxrep_fname, flags.load_maxrep);
        cerr << "DONE (" << time_usage.reg["ms_maxrep"] / 1000 << " seconds)" << endl;
    }

    Slices<size_type> slices(Query::query_length(tspec.s_fname), flags.nthreads);
    /* build MS */
    auto runs_start = timer::now();
    std::vector<std::future<Interval>> results(flags.nthreads);
    for(size_type i=0; i<flags.nthreads; i++){
        cerr << " ** launching ms computation over : " << slices.repr(i) << endl;
        results[i] = std::async(std::launch::async, fill_ms_slice_thread, i, slices[i], st.root(), flags, tspec);
    }
    for(size_type i=0; i<flags.nthreads; i++)
        results[i].get();
    
    auto runs_stop = timer::now();
    time_usage.reg["ms_bvector"] = std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count();

    cerr << " * total ms length : " << ms_vec.ms_size()  << " (with |t| = " << Query::query_length(tspec.s_fname) << ")" << endl;
    cerr << "DONE (" << time_usage.reg["ms_bvector"] / 1000 << " seconds)" << endl;
}



void comp(const InputSpec& tspec, InputSpec& S_fwd, counter_t& time_usage, InputFlags& flags){
    auto comp_start = timer::now();
    /* prepare global data structures */
    ms_vec = msvec_t(Query::query_length(tspec.s_fname), flags.nthreads);

    /* build runs */
    build_runs_ohleb(flags, S_fwd, tspec, time_usage);

    /* build ms */
    build_ms_ohleb(flags, S_fwd, tspec, time_usage);
    time_usage.register_now("total_time", comp_start);

    if(flags.time_usage){
        cerr << "dumping reports" << endl;
        cout << "len_s,len_t,item,value" << endl;
        if(flags.time_usage){
            for(auto item : time_usage.reg)
                cout << st.size() - 1 << "," << ms_vec.runs.size() << "," << item.first << "," << item.second << endl;
        }
    }

    if(flags.answer){
        ms_vec.show_MS(cout);
        cout << endl;
    }
}


int main(int argc, char **argv){
    OptParser input(argc, argv);
    InputSpec sfwd_spec, tspec;
    InputFlags flags;
    counter_t time_usage{};

    if(argc == 1){
        const string base_dir = {"/home/brt/code/matching_statistics/indexed_ms/fast_ms/tests/"};
        tspec = InputSpec(base_dir + "a.t");
        sfwd_spec = InputSpec(base_dir + "a.s");
        flags = InputFlags(true,  // use double rank
                           false, // lazy_wl
                           false, // rank-and-fail
                           false, // use maxrep vanilla
                           false, // use maxrep rank&check
                           false, // lca_parents
                           false, // time
                           true,  // ans
                           false, // avg
                           false, // load CST
                           false, // load MAXREP
                           2      // nthreads
                           );
    } else {
        tspec = InputSpec(input.getCmdOption("-t_path"));
        sfwd_spec = InputSpec(input.getCmdOption("-s_path"));
        flags = InputFlags(input);
    }
    comp(tspec, sfwd_spec, time_usage, flags);
    return 0;
}

