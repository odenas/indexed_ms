#include <iostream>
#include <fstream>
#include <vector>
#include <string>


#include "fd_ms/help.hpp"
#include "fd_ms/opt_parser.hpp"
#include "fd_ms/input_spec.hpp"
#include "fd_ms/counter.hpp"
#include "fd_ms/query.hpp"
#include "fd_ms/counter.hpp"
#include "fd_ms/stree_sct3.hpp"
#include "fd_ms/maxrep_vector.hpp"
#include "fd_ms/runs_vector.hpp"
#include "fd_ms/ms_vector.hpp"

using namespace std;
using namespace fdms;

typedef StreeOhleb<> cst_t;
typedef typename cst_t::node_type node_type;
typedef typename cst_t::size_type size_type;
typedef typename cst_t::char_type char_type;
typedef sdsl::bit_vector bitvec_t;
typedef Maxrep<cst_t, sdsl::bit_vector> maxrep_t;
typedef Counter<size_type> counter_t;

typedef typename runs_vector<cst_t>::pseq_method_t pseq_method_t;
typedef typename runs_vector<cst_t>::wl_method_t1 wl_method_t1;
typedef typename runs_vector<cst_t>::wl_method_t2 wl_method_t2;

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
    load_maxrep{f.load_maxrep}
    {
    }

    InputFlags(bool double_rank, bool lazy_wl, bool use_rank_fail, bool use_maxrep_vanilla, bool use_maxrep_rc, bool lca_parents,
            bool time_, bool ans, bool avg,
            bool load_stree, bool load_maxrep) :
    double_rank{double_rank},
    lazy{lazy_wl},
    rank_fail{use_rank_fail},
    use_maxrep_vanilla{use_maxrep_vanilla}, use_maxrep_rc{use_maxrep_rc},
    lca_parents{lca_parents},
    time_usage{time_},
    answer{ans}, avg{avg},
    load_stree{load_stree},
    load_maxrep{load_maxrep}
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
    load_maxrep{input.getCmdOption("-load_maxrep") == "1"} // load MAXREP of S'
    {
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

    pseq_method_t get_pseq_method() {
        return (lca_parents ? &runs_vector<cst_t>::lca_parent : &runs_vector<cst_t>::parent_sequence);
    }

    size_type buffer_size(const InputSpec& ispec) const {
        size_type t_length = ispec.t_size();
        return (t_length > 1000 ? t_length / 100 : t_length);
    }
};

void build_runs(const InputSpec& ispec, counter_t& time_usage, InputFlags& flags) {
    cerr << "building RUNS ... " << endl;
    size_type buffer_size = flags.buffer_size(ispec);

    time_usage.reg["runs_cst"] = cst_t::load_or_build(st, ispec, false, flags.load_stree);
    cerr << "DONE (" << time_usage.reg["runs_cst"] / 1000 << " seconds, " << st.size() << " leaves)" << endl;

    /* compute RUNS */
    if (flags.use_maxrep()) {
        /* build the maxrep vector */
        time_usage.reg["ms_maxrep"] = maxrep_t::load_or_build(maxrep, st, ispec.fwd_maxrep_fname, flags.load_maxrep);
        cerr << "DONE (" << time_usage.reg["ms_maxrep"] / 1000 << " seconds)" << endl;

        auto start = timer::now();
        runs_vector<cst_t>::dump(ispec, st, flags.get_mrep_wl_method(), maxrep, buffer_size);
        time_usage.register_now("runs_bvector", start);
    } else {
        auto start = timer::now();
        runs_vector<cst_t>::dump(ispec, st, flags.get_wl_method(), flags.get_pseq_method(), buffer_size);
        time_usage.register_now("runs_bvector", start);
    }
}

void build_ms(const InputSpec& ispec, counter_t& time_usage, InputFlags& flags) {
    cerr << "building MS ... " << endl;
    size_type buffer_size = flags.buffer_size(ispec);

    /* build the CST */
    time_usage.reg["ms_cst"] = cst_t::load_or_build(st, ispec, true, flags.load_stree);
    cerr << "DONE (" << time_usage.reg["ms_cst"] / 1000 << " seconds, " << st.size() << " leaves)" << endl;

    /* compute MS */
    if (flags.use_maxrep()) {
        /* build the maxrep vector */
        time_usage.reg["ms_maxrep"] = maxrep_t::load_or_build(maxrep, st, ispec.rev_maxrep_fname, flags.load_maxrep);
        cerr << "DONE (" << time_usage.reg["ms_maxrep"] / 1000 << " seconds)" << endl;
        auto start = timer::now();
        ms_vector<cst_t>::dump(ispec, st,
                flags.get_mrep_wl_method(), maxrep, buffer_size);
        time_usage.register_now("ms_bvector", start);
    } else {
        auto start = timer::now();
        ms_vector<cst_t>::dump(ispec, st,
                flags.get_wl_method(), flags.get_pseq_method(), buffer_size);
        time_usage.register_now("ms_bvector", start);
    }
    (cerr << " * total ms length : " << ms_vector<cst_t>::size(ispec) <<
            " (with |t| = " << runs_vector<cst_t>::size(ispec) << ")" << endl);
}

void comp(const InputSpec& ispec, counter_t& time_usage, InputFlags& flags) {
    auto start = timer::now();
    build_runs(ispec, time_usage, flags);
    time_usage.register_now("runs_build", start, true);
    //runs_vector<cst_t>::show(ispec.runs_fname, cerr);

    start = timer::now();
    build_ms(ispec, time_usage, flags);
    time_usage.register_now("ms_build", start, true);
    //runs_vector<cst_t>::show(ispec.ms_fname, cerr);

    if (flags.time_usage) {
        cerr << "dumping reports ..." << endl;
        cout << "len_s,len_t,item,value" << endl;
        for (auto item : time_usage.reg)
            (cout << st.size() - 1 << ","
                << ms_vector<cst_t>::size(ispec) << ","
                << item.first << "," << item.second << endl);
    }

    if (flags.answer) {
        ms_vector<cst_t>::show_MS(ispec, cout);
        cout << endl;
    } else if (flags.avg)
        cout << ms_vector<cst_t>::avg_matching_statistics(ispec) << endl;
}

int main(int argc, char **argv) {
    OptParser input(argc, argv);
    InputSpec ispec;
    InputFlags flags;
    counter_t time_usage{};

    //(cerr << "Deprecated. Use 'matching_stats_parallel.x' instead" << endl); exit(0);
    if (argc == 1) {
        (cerr << "Compute the matching statistics of the given inputs.\n"
              << "Args:\n"
              << help__s_path << help__t_path
              << "\t-double_rank 1: use the double rank strategy\n"
              << "\t-lazy_wl 1: use the lazy weiner link strategy\n"
              << "\t-rank_fail 1: use the rank and fail strategy\n"
              << "\t-use_maxrep_rc 1: use the maxrep rank and check strategy\n"
              << "\t-use_maxrep_vanilla 1: use the maxrep bit vector\n"
              << "\t-lca_parents 1: use the lca parents strategy\n"
              << help__time_usage
              << "\t-answer 1: Dump the answer in the standard output\n"
              << help__load_cst
              << help__load_maxrep
              << endl);
        exit(0);

        const string base_dir = {"/home/brt/code/matching_statistics/indexed_ms/fast_ms/tests/"};
        ispec = InputSpec(base_dir + "a.s", base_dir + "a.t");
        flags = InputFlags(true, // use double rank
                false, // lazy_wl
                true, // rank-and-fail
                false, // use maxrep vanilla
                true, // use maxrep rank&check
                false, // lca_parents
                false, // time
                false, // ans
                false, // avg
                false, // load CST
                false // load MAXREP
                );
    } else {
        ispec = InputSpec(input.getCmdOption("-s_path"), input.getCmdOption("-t_path"));
        flags = InputFlags(input);
    }
    comp(ispec, time_usage, flags);
    return 0;
}

