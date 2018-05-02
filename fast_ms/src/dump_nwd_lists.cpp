#include <iostream>
#include <fstream>
#include <vector>
#include <string>


#include "fd_ms/input_spec.hpp"
#include "fd_ms/opt_parser.hpp"

#include "fd_ms/stree_sct3.hpp"
#include "fd_ms/parent_depth_list.hpp"


using namespace std;
using namespace fdms;


class InputFlags{
public:
    bool check, load_cst;
    size_t sample_freq;
    
    InputFlags(){}
    
    InputFlags(const bool check, const bool load_cst, const size_t sample_freq) :
    check{check}, load_cst{load_cst}, sample_freq{sample_freq} {}
    
    InputFlags(const InputFlags& i) : check{i.check}, load_cst{i.load_cst}, sample_freq{i.sample_freq} {}
    
    InputFlags(const OptParser args){
        check = (args.getCmdOption("-check") == "1");
        load_cst = (args.getCmdOption("-load_cst") == "1");
        sample_freq = (static_cast<size_t>(std::stoi(args.getCmdOption("-sample_freq"))));
    }
};



void comp(InputSpec& S_fwd, InputFlags& flags){
    string s = S_fwd.load_s();
    
    /* build stree */
    StreeOhleb<> st;
    size_type load_cst_time = load_or_build(st, s, S_fwd.rev_cst_fname, flags.load_cst);
    cerr << "DONE (" << load_cst_time / 1000 << "seconds)" << endl;

    NwdList nwd_lst{st, flags.sample_freq};
    
    
    if(flags.check)
        nwd_lst.check(st, s, flags.sample_freq);
    nwd_lst.write_bin(S_fwd.fwd_nwdlst_fname);
}

int main(int argc, char **argv){
    OptParser input(argc, argv);
    InputSpec sfwd_spec;
    InputFlags flags;
    
    if(argc == 1){
        const string base_dir = {"/Users/denas/projects/matching_statistics/indexed_ms/tests/code_test/maxrep_inputs/"};
        flags = InputFlags(true, false, 2);
        sfwd_spec = InputSpec(base_dir + "rnd_20s_dis_10t_abcd.s");
    } else {
        flags = InputFlags(input);
        sfwd_spec = InputSpec(input.getCmdOption("-s_path"));
    }
    comp(sfwd_spec, flags);
    return 0;
}
