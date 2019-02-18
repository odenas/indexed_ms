#include <iostream>
#include <fstream>
#include <vector>
#include <string>


#include "fd_ms/input_spec.hpp"
#include "fd_ms/opt_parser.hpp"

#include "fd_ms/stree_sct3.hpp"
#include "fd_ms/maxrep_vector.hpp"
#include "fd_ms/edge_list.hpp"


using namespace std;
using namespace fdms;

typedef StreeOhleb<>     cst_t;
typedef sdsl::bit_vector bitvector_t;

class InputFlags{
public:
    bool check, load_cst, load_maxrep;
    size_t sample_freq;

    InputFlags(){}

    InputFlags(const bool check, const bool load_cst, const bool load_maxrep, const size_t sample_freq) :
    check{check}, load_cst{load_cst}, sample_freq{sample_freq}, load_maxrep{load_maxrep} {}

    InputFlags(const InputFlags& i) : check{i.check}, load_cst{i.load_cst}, load_maxrep{i.load_maxrep}, sample_freq{i.sample_freq} {}

    InputFlags(const OptParser args){
        check = (args.getCmdOption("-check") == "1");
        load_cst = (args.getCmdOption("-load_cst") == "1");
        load_maxrep = (args.getCmdOption("-load_maxrep") == "1");
        sample_freq = (static_cast<size_t>(std::stoi(args.getCmdOption("-sample_freq"))));
    }
};


void comp(InputSpec& S_fwd, InputFlags& flags){
        /* build stree */
    cst_t st;
    size_type load_cst_time = cst_t::load_or_build(st, S_fwd, true, flags.load_cst);
    cerr << "DONE (" << load_cst_time / 1000 << "seconds)" << endl;

    EdgeList e{st, flags.sample_freq};
    if(flags.check){
        cerr << e.repr() << endl;
        Maxrep<cst_t, bitvector_t> maxrep;
        Maxrep<cst_t, bitvector_t>::load_or_build(maxrep, st, S_fwd.rev_maxrep_fname, flags.load_maxrep);
        cerr << "OK" << endl;
        e.check(maxrep, (size_t)st.csa.sigma);
    }
    e.write_bin(S_fwd.rev_elst_fname);
}

int main(int argc, char **argv){
    OptParser input(argc, argv);
    InputSpec ispec;
    InputFlags flags;

    if(argc == 1){
        //const string base_dir = {"/Users/denas/projects/matching_statistics/indexed_ms/tests/code_test/maxrep_inputs/"};
        //flags = InputFlags(true, false, false, 2);
        //ispec = InputSpec(base_dir + "rnd_20s_dis_10t_abcd.s", "");
        (cerr << "Experimental." << endl); exit(0);
    } else {
        flags = InputFlags(input);
        ispec = InputSpec(input.getCmdOption("-s_path"), "");
    }
    comp(ispec, flags);
    return 0;
}
