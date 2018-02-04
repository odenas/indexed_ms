/*
 dump index data structures (maxrep and cst) for given input string
 */


#include <iostream>
#include <fstream>
#include <vector>

#include "fd_ms/input_spec.hpp"
#include "fd_ms/opt_parser.hpp"
#include "fd_ms/stree_sct3.hpp"
#include "fd_ms/maxrep_vector.hpp"


using namespace std;
using namespace fdms;


class InputFlags{
public:
    bool load_fwd_cst, load_rev_cst;

    InputFlags(){}

    InputFlags(const bool fwd_cst, const bool rev_cst) :
        load_fwd_cst{fwd_cst}, load_rev_cst{rev_cst} {}

    InputFlags(const InputFlags& i) :
        load_fwd_cst{i.load_fwd_cst}, load_rev_cst{i.load_rev_cst} {}

    InputFlags(const OptParser args){
        load_fwd_cst = (args.getCmdOption("-load_fwd_cst") == "1");
        load_rev_cst = (args.getCmdOption("-load_rev_cst") == "1");
    }

    void help(const string exec_name) const {
        cerr << exec_name << endl << "\t";
        cerr << "[-load_fwd_cst 0/1]" << endl << "\t";
        cerr << "[-load_rev_cst 0/1]" << endl << "\t";
        cerr << "-s_path <s path>" << endl;
    }
};

void comp(const InputSpec& ispec){
    StreeOhleb<> st;

    cerr << "dumping rwd/rev cst/maxrep data structures ..." << endl;
    cerr << " * forward ... " << endl;
    cerr << " ** loading the index from  " << ispec.s_fname << " ..." << endl;
    string s = ispec.load_s(false);
    cerr << " ** building and dumping cst to " << ispec.fwd_cst_fname << " ..." << endl;
    sdsl::construct_im(st, s, 1);
    sdsl::store_to_file(st, ispec.fwd_cst_fname);

    cerr << " ** building and dumping maxrep to " << ispec.fwd_maxrep_fname << " ..." << endl;
    Maxrep<StreeOhleb<>, sdsl::bit_vector> mr(st, false);
    mr.dump_vec(ispec.fwd_maxrep_fname);

    cerr << " * reverse ... " << endl;
    cerr << " ** reversing the index string ..." << endl;
    ispec.reverse_in_place(s);
    cerr << " ** building and dumping cst to " << ispec.rev_cst_fname << " ..." << endl;
    sdsl::construct_im(st, s, 1);
    sdsl::store_to_file(st, ispec.rev_cst_fname);

    cerr << " ** building and dumping maxrep to " << ispec.rev_maxrep_fname << " ..." << endl;
    mr = Maxrep<StreeOhleb<>, sdsl::bit_vector>(st, false);
    mr.dump_vec(ispec.rev_maxrep_fname);
}


int main(int argc, char **argv){
    OptParser input(argc, argv);
    InputSpec s_spec;
    InputFlags flags;

    if(argc == 1){
        const string base_dir = {"/Users/denas/projects/matching_statistics/indexed_ms/tests/datasets/testing/"};
        s_spec = InputSpec(base_dir + "mut_200s_64t_15.s");
        flags = InputFlags(false, false);
    } else {
        flags = InputFlags(input);
        s_spec = InputSpec(input.getCmdOption("-s_path"));
    }
    comp(s_spec);
    return 0;
}
