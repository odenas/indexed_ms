#include <iostream>
#include <string>

#include "fd_ms/input_spec.hpp"
#include "fd_ms/opt_parser.hpp"
#include "fd_ms/stree_sct3.hpp"
#include "fd_ms/runs_vector.hpp"
#include "fd_ms/ms_vector.hpp"
#include "fd_ms/maxrep_vector.hpp"


using namespace std;
using namespace fdms;

typedef typename StreeOhleb<>::node_type node_type;
typedef typename StreeOhleb<>::size_type size_type;
typedef uint8_t char_type;

typedef StreeOhleb<> cst_t;
typedef sdsl::bit_vector bit_vector;
cst_t st;
Maxrep<cst_t, bit_vector> maxrep;

class InputFlags {
public:
    bool load_cst, load_maxrep;

    InputFlags() {
    }

    InputFlags(const bool load_cst_, const bool load_maxrep_) :
    load_cst{load_cst_}, load_maxrep{load_maxrep_}
    {
    }

    InputFlags(const InputFlags& i) :
    load_cst{i.load_cst}, load_maxrep{i.load_maxrep}
    {
    }

    InputFlags(const OptParser& args) :
    load_cst{(args.getCmdOption("-load_cst") == "1")}, load_maxrep{(args.getCmdOption("-load_maxrep") == "1")}
    {
    }
};

void comp(InputSpec& ispec, InputFlags& flags) {
    Stats<cst_t, Maxrep <cst_t, bit_vector >> stats(
            sdsl::util::file_size(ispec.s_fname),
            sdsl::util::file_size(ispec.t_fname));

    /* build runs */
    cerr << "build runs ... " << endl;
    cst_t::load_or_build(st, ispec, false, flags.load_cst);
    runs_vector<cst_t>::fill_runs(
            stats, ispec, st,
            &cst_t::double_rank_fail_wl, runs_vector<cst_t>::parent_sequence,
            1024 * 1024);

    /* build ms */
    cerr << "build ms ... " << endl;
    cst_t::load_or_build(st, ispec, true, flags.load_cst);
    Maxrep<cst_t, bit_vector>::load_or_build(maxrep, st, ispec.rev_maxrep_fname, flags.load_maxrep);
    ms_vector<cst_t>::fill_ms(
            stats, ispec, st,
            &cst_t::double_rank_fail_wl, ms_vector<cst_t>::parent_sequence,
            maxrep, 1024 * 1024);

    cerr << "dumping results ... " << endl;
    stats.dump_stats(cout);
}

int main(int argc, char **argv) {
    OptParser input(argc, argv);
    InputSpec ispec;
    InputFlags flags;

    (cerr << "Experimental." << endl); exit(0);

    if (argc == 1) {
        const string base_dir = {"/home/brt/code/matching_statistics/indexed_ms/tests/input_stats_data/exploration/index_based_query/"};
        ispec = InputSpec(base_dir + "rep_2.s", base_dir + "rep_2.t");
        flags = InputFlags(true, true);
    } else {
        ispec = InputSpec(input.getCmdOption("-s_path"), input.getCmdOption("-t_path"));
        flags = InputFlags(input);
    }
    comp(ispec, flags);
}
