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


class InputFlags{
public:
    bool load_cst, load_maxrep;

    InputFlags() {}

    InputFlags(const bool load_cst_, const bool load_maxrep_) :
    	load_cst{load_cst_}, load_maxrep{load_maxrep_} {}

    InputFlags(const InputFlags& i) :
    	load_cst{i.load_cst}, load_maxrep{i.load_maxrep} {}

    InputFlags(const OptParser& args) :
    	load_cst {(args.getCmdOption("-load_cst") == "1")}, load_maxrep{(args.getCmdOption("-load_maxrep") == "1")} {}
};


void comp(InputSpec& tspec, InputSpec& s_fwd, const string& out_path, InputFlags& flags){
    map<size_type, size_type> consecutive_runs_parent_calls, consecutive_ms_wl_calls, consecutive_ms_parent_calls;
    map<string, size_type> ms_wl_node_prop, runs_wl_node_prop, ms_rank_calls;

    /* prepare global data structures */
    cerr << "loading input ... " << endl;
    t = tspec.load_s();
    s = s_fwd.load_s();
    ms_vec = MsVectors<StreeOhleb<>, sdsl::bit_vector>(t.size());
    Stats<StreeOhleb<>, Maxrep<StreeOhleb<>, sdsl::bit_vector>> stats(s.size(), t.size());

    /* build runs */
    cerr << "build runs ... " << endl;
    load_or_build(st, s, s_fwd.fwd_cst_fname, flags.load_cst);
    ms_vec.fill_runs(stats, t, st,
    		&StreeOhleb<>::double_rank_fail_wl,
    		&MsVectors<StreeOhleb<>, sdsl::bit_vector>::parent_sequence);

    /* reverse s */
    cerr << "reversing index of length " << s.size() << " ... " << endl;
    InputSpec::reverse_in_place(s);

    /* build ms */
    cerr << "build ms ... " << endl;
    load_or_build(st, s, s_fwd.rev_cst_fname, flags.load_cst);
    Maxrep<StreeOhleb<>, sdsl::bit_vector>::load_or_build(maxrep, st, s_fwd.rev_maxrep_fname, flags.load_maxrep);
    ms_vec.fill_ms(stats, t, st,
    		&StreeOhleb<>::double_rank_fail_wl,
    		&MsVectors<StreeOhleb<>, sdsl::bit_vector>::parent_sequence,
    		maxrep);

    cerr << "dumping results ... " << endl;
	stats.dump_stats(cout);
}

int main(int argc, char **argv){
    OptParser input(argc, argv);
    InputSpec sfwd_spec, tspec;
    InputFlags flags;
    string out_path;

    if(argc == 1){
        const string base_dir = {"/home/brt/code/matching_statistics/indexed_ms/tests/datasets/testing/"};
        tspec = InputSpec(base_dir + "rnd_200_32.t");
        sfwd_spec = InputSpec(base_dir + "rnd_200_32.s");
        out_path = "0";
        flags = InputFlags(false, false);
    } else {
        tspec = InputSpec(input.getCmdOption("-t_path"));
        sfwd_spec = InputSpec(input.getCmdOption("-s_path"));
        flags = InputFlags(input);
    }
    comp(tspec, sfwd_spec, out_path, flags);
    return 0;
}
