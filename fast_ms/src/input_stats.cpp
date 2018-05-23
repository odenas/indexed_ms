#include <iostream>
#include <string>

#include <sdsl/suffix_trees.hpp>

#include "fd_ms/opt_parser.hpp"
#include "fd_ms/input_spec.hpp"
#include "fd_ms/stree_sct3.hpp"
#include "fd_ms/runs_ms.hpp"
#include "fd_ms/maxrep_vector.hpp"


using namespace std;
using namespace fdms;

typedef typename StreeOhleb<>::node_type node_type;
typedef typename StreeOhleb<>::size_type size_type;
typedef uint8_t char_type;

// global data structures
string t, s;
StreeOhleb<> st;
Maxrep<StreeOhleb<>, sdsl::bit_vector> maxrep;
MsVectors<StreeOhleb<>, sdsl::bit_vector> ms_vec;

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
    /* prepare global data structures */
    cerr << "loading input ... " << endl;
    t = tspec.load_s();
    s = s_fwd.load_s();
    ms_vec = MsVectors<StreeOhleb<>, sdsl::bit_vector>(t.size());
    Stats<StreeOhleb<>, Maxrep < StreeOhleb<>, sdsl::bit_vector >> stats(s.size(), t.size());

    /* build runs */
    cerr << "build runs ... " << endl;
    StreeOhleb<>::load_or_build(st, ispec, false, flags.load_cst);
    ms_vec.fill_runs(stats, t, st,
            &StreeOhleb<>::double_rank_fail_wl,
            &MsVectors<StreeOhleb<>, sdsl::bit_vector>::parent_sequence);

    /* build ms */
    cerr << "build ms ... " << endl;
    StreeOhleb<>::load_or_build(st, ispec, true, flags.load_cst);
    Maxrep<StreeOhleb<>, sdsl::bit_vector>::load_or_build(maxrep, st, s_fwd.rev_maxrep_fname, flags.load_maxrep);
    ms_vec.fill_ms(stats, t, st,
            &StreeOhleb<>::double_rank_fail_wl,
            &MsVectors<StreeOhleb<>, sdsl::bit_vector>::parent_sequence,
            maxrep);

    cerr << "dumping results ... " << endl;
    stats.dump_stats(cout);
}

int main(int argc, char **argv) {
    OptParser input(argc, argv);
    InputSpec ispec;
    InputFlags flags;

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
