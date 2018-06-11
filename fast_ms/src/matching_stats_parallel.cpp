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
#include "fd_ms/p_runs_vector.hpp"

using namespace std;
using namespace fdms;

typedef StreeOhleb<> cst_t;
typedef typename cst_t::node_type node_type;
typedef typename cst_t::size_type size_type;
typedef typename cst_t::char_type char_type;
typedef sdsl::bit_vector bitvec_t;
typedef Maxrep<cst_t, bitvec_t> maxrep_t;
typedef Counter<size_type> counter_t;

typedef typename p_runs_vector<cst_t>::pseq_method_t pseq_method_t;
typedef typename p_runs_vector<cst_t>::wl_method_t1 wl_method_t1;
typedef typename p_runs_vector<cst_t>::wl_method_t2 wl_method_t2;
typedef typename p_runs_vector<cst_t>::pair_t pair_t;
typedef typename p_runs_vector<cst_t>::alg_state_t runs_rt;

cst_t st;
maxrep_t maxrep;

class InputFlags {
private:

    void check() const {
        if (use_maxrep_rc && use_maxrep_vanilla) {
            cerr << "use_maxrep_rc and use_maxrep_vanilla cannot be active at the same time" << endl;
            exit(1);
        }
        if (use_maxrep() && !(double_rank && rank_fail)) {
            cerr << "use_maxrep_xx goes with double rank and fail" << endl;
            exit(1);
        }
        if (use_maxrep() && lazy) {
            cerr << "lazy and use_maxrep_xx cannot be active at the same time" << endl;
            cerr << "use_maxrep_xx goes with double rank and fail" << endl;
            exit(1);
        }
        if (use_maxrep() && !double_rank) {
            cerr << "single_rank and use_maxrep_xx cannot be active at the same time" << endl;
            cerr << "use_maxrep_xx goes with double rank and fail" << endl;
            exit(1);
        }
        if (use_maxrep() && !rank_fail) {
            cerr << "no_fail and use_maxrep_xx cannot be active at the same time" << endl;
            cerr << "use_maxrep_xx goes with double rank and fail" << endl;
            exit(1);
        }
        if (rank_fail && !double_rank) {
            cerr << "single_rank and rank_fail cannot be active at the same time" << endl;
            exit(1);
        }
        if (answer && avg) {
            cerr << "answer and avg cannot be active at the same time" << endl;
            exit(1);
        }
    }

public:
    bool double_rank, lazy, rank_fail, use_maxrep_vanilla, use_maxrep_rc, lca_parents;
    bool time_usage, answer, avg;
    bool load_stree, load_maxrep;
    size_t nthreads;

    InputFlags() {
    }

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
    nthreads{f.nthreads}
    {
    }

    InputFlags(bool double_rank, bool lazy_wl, bool use_rank_fail, bool use_maxrep_vanilla, bool use_maxrep_rc, bool lca_parents,
            bool time_, bool ans, bool avg,
            bool load_stree, bool load_maxrep, size_t nthreads) :
    double_rank{double_rank},
    lazy{lazy_wl},
    rank_fail{use_rank_fail},
    use_maxrep_vanilla{use_maxrep_vanilla}, use_maxrep_rc{use_maxrep_rc},
    lca_parents{lca_parents},
    time_usage{time_},
    answer{ans}, avg{avg},
    load_stree{load_stree},
    load_maxrep{load_maxrep},
    nthreads{nthreads}
    {
        check();
    }

    InputFlags(OptParser input) :
    double_rank{input.getCmdOption("-double_rank") == "1"}, // use double rank
    lazy{input.getCmdOption("-lazy_wl") == "1"}, // lazy winer links
    rank_fail{input.getCmdOption("-rank_fail") == "1"}, // use the rank-and-fail strategy
    use_maxrep_rc{input.getCmdOption("-use_maxrep_rc") == "1"}, // use the maxrep vector with rank_and_check
    use_maxrep_vanilla{input.getCmdOption("-use_maxrep_vanilla") == "1"}, // use the maxrep vector the vanilla way
    lca_parents{input.getCmdOption("-lca_parents") == "1"}, // use lca insted of conscutive parent calls
    time_usage{input.getCmdOption("-time_usage") == "1"}, // time usage
    answer{input.getCmdOption("-answer") == "1"}, // answer
    avg{input.getCmdOption("-avg") == "1"}, // average matching statistics
    load_stree{input.getCmdOption("-load_cst") == "1"}, // load CST of S and S'
    load_maxrep{input.getCmdOption("-load_maxrep") == "1"}, // load MAXREP of S'
    nthreads{static_cast<size_t> (std::stoi(input.getCmdOption("-nthreads")))}
    {
        nthreads = (nthreads <= 0 ? 1 : nthreads);
        check();
    }

    bool use_maxrep() const {
        return (use_maxrep_rc || use_maxrep_vanilla);
    }

    wl_method_t1 get_wl_method() {
        if (double_rank) {
            if (lazy)
                return (rank_fail ?
                    &cst_t::lazy_double_rank_fail_wl :
                    &cst_t::lazy_double_rank_nofail_wl);
            return (rank_fail ?
                    &cst_t::double_rank_fail_wl :
                    &cst_t::double_rank_nofail_wl);
        } else {
            return (lazy ?
                    &cst_t::lazy_single_rank_wl :
                    &cst_t::single_rank_wl);
        }
    }

    wl_method_t2 get_mrep_wl_method() {
        return (use_maxrep_rc ?
                &cst_t::double_rank_fail_wl_mrep_rc :
                &cst_t::double_rank_fail_wl_mrep_vanilla);
    }
    
    pseq_method_t get_pseq_method(){
        return (lca_parents ? 
            &p_runs_vector<cst_t>::lca_parent : 
            &p_runs_vector<cst_t>::parent_sequence);
    }

};

runs_rt fill_runs_slice_thread1(const InputSpec& ispec,
                               const size_type thread_id,
                               const pair_t slice, node_type v,
                               InputFlags flags) {

    flags.lazy = false;  // runs does not support laziness
    return (p_runs_vector<cst_t>(flags.nthreads, thread_id, slice)
            .fill_slice(ispec, st, flags.get_wl_method(), flags.get_pseq_method(), v, 1024));
}

int fill_runs_slice_thread2(const InputSpec& ispec,
                           const size_type slice_idx,
                           pair_t slice, node_type v,
                           const Slices<size_type>& slices,
                           InputFlags flags) {

    flags.lazy = false;  // runs does not support laziness
    return (p_runs_vector<cst_t>(flags.nthreads, slice_idx, slice)
            .fill_inter_slice(ispec, st,
                              flags.get_wl_method(), flags.get_pseq_method(),
                              v, slice_idx, slices,
                              1024));
}

vector<runs_rt> aa(const vector<runs_rt> v, const Slices<size_type> slices) {
    vector<runs_rt> u;
    u.reserve(v.size());

    int i = 0, j = 1, n = v.size();
    while (j < n) {
        runs_rt prev_slice = v[i], next_slice = v[j];
        while (get<0>(next_slice) == get<0>(slices.slices[j]) + 1) { // j-th slice had a full match
            j++;
        }
        next_slice = v[j - (j == n)];
        u.push_back(make_tuple(get<0>(prev_slice), get<1>(next_slice), get<2>(next_slice)));
        i = j;
        j += 1;
    }
    return u;

}

void build_runs(const InputSpec& ispec, counter_t& time_usage, const InputFlags& flags) {
    size_type t_length = Query::query_length(ispec.t_fname);
    size_type buffer_size = (t_length > 1000 ? t_length / 100 : t_length);
    
    cerr << "building RUNS ... " << endl;

    /* build the CST */
    time_usage.reg["runs_cst"] = cst_t::load_or_build(st, ispec, false, flags.load_stree);
    cerr << "DONE (" << time_usage.reg["runs_cst"] / 1000 << " seconds, " << st.size() << " leaves)" << endl;

    /* compute RUNS */
    auto runs_start = timer::now();
    std::vector<std::future < runs_rt >> results(flags.nthreads);
    Slices<size_type> slices(t_length, flags.nthreads);

    /* open connection to the query string */
    Query_rev t{ispec.s_fname, (size_t) buffer_size};
    for (size_type i = 0; i < flags.nthreads; i++) {
        node_type v = st.double_rank_nofail_wl(st.root(), t[slices[i].second - 1]); // stree node
        cerr << " ** launching runs computation over : " << slices.repr(i) << " ";
        cerr << "(" << v.i << ", " << v.j << ")" << endl;
        results[i] = std::async(std::launch::async, fill_runs_slice_thread1,
                ispec, i, slices[i], v, flags);
    }
    vector<runs_rt> runs_results(flags.nthreads);
    for (size_type i = 0; i < flags.nthreads; i++) {
        runs_results[i] = results[i].get();
        cerr << " *** [" << get<0>(runs_results[i]) << " .. " << get<1>(runs_results[i]) << ")" << endl;
    }

    vector<runs_rt> merge_idx = aa(runs_results, slices);
    cerr << " * correcting edges over " << merge_idx.size() << " threads ... " << endl;
    std::vector<std::future < int >> results2(merge_idx.size());
    for (int i = 0; i < (int) merge_idx.size(); i++) {
        cerr << " ** ([" << get<0>(merge_idx[i]) << ", " << get<1>(merge_idx[i]) << "), " << "(" << get<2>(merge_idx[i]).i << ", " << get<2>(merge_idx[i]).j << ")) " << endl;
        results2[i] = std::async(std::launch::async, fill_runs_slice_thread2,
                ispec, 
                slices.slice_idx(get<1>(merge_idx[i])),
                make_pair(get<0>(merge_idx[i]) - 1, get<1>(merge_idx[i])), get<2>(merge_idx[i]), 
                slices,
                flags);
    }
    for (int i = 0; i < merge_idx.size(); i++)
        results2[i].get();

    cerr << " * merging into " << ispec.runs_fname << " ... " << endl;
    p_runs_vector<cst_t>::merge(ispec, slices, 1024);
    time_usage.register_now("runs_bvector", runs_start);
    
    p_runs_vector<cst_t>::show(ispec.runs_fname, cerr);
}


void comp(const InputSpec& ispec, counter_t& time_usage, const InputFlags& flags) {
    auto comp_start = timer::now();
    /* prepare global data structures */

    build_runs(ispec, time_usage, flags);
    return;
}

int main(int argc, char **argv) {
    OptParser input(argc, argv);
    InputSpec ispec;
    InputFlags flags;
    counter_t time_usage{};

    if (argc == 1) {
        const string base_dir = {"/home/brt/code/matching_statistics/indexed_ms/fast_ms/tests/"};
        ispec = InputSpec(base_dir + "a.s", base_dir + "a.t");
        flags = InputFlags(true, // use double rank
                false, // lazy_wl
                false, // rank-and-fail
                false, // use maxrep vanilla
                false, // use maxrep rank&check
                false, // lca_parents
                false, // time
                true, // ans
                false, // avg
                false, // load CST
                false, // load MAXREP
                2 // nthreads
                );
    } else {
        ispec = InputSpec(input.getCmdOption("-s_path"), input.getCmdOption("-t_path"));
        flags = InputFlags(input);
    }
    comp(ispec, time_usage, flags);
    return 0;
}

