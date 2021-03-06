/*
 dump a CST to a file. to be used for later use
 */

#include <iostream>
#include <fstream>
#include <vector>

#include "fd_ms/input_spec.hpp"
#include "fd_ms/opt_parser.hpp"
#include "fd_ms/stree_sct3.hpp"
#include "fd_ms/help.hpp"


using namespace std;
using namespace fdms;
using timer = std::chrono::high_resolution_clock;
typedef StreeOhleb<>                cst_t;


void dump(const StreeOhleb<>& st, const string fname){
    auto start = timer::now();
    cerr << " * dumping the CST to " << fname << " ";
    sdsl::store_to_file(st, fname);
    auto stop = timer::now();
    cerr << "DONE (" << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << "seconds)" << endl;
}

void comp(const InputSpec& s_spec){
    StreeOhleb<> st;
    StreeOhleb<>::size_type load_cst_time = cst_t::load_or_build(st, s_spec, false, false);
    cerr << "DONE (" << load_cst_time / 1000 << "seconds)" << endl;
    dump(st, s_spec.fwd_cst_fname);

    load_cst_time = cst_t::load_or_build(st, s_spec, true, false);
    cerr << "DONE (" << load_cst_time / 1000 << "seconds)" << endl;
    dump(st, s_spec.rev_cst_fname);
}


int main(int argc, char **argv){
    OptParser input(argc, argv);
    InputSpec ispec;

    if(argc == 1){
        (cerr << "Dump the cst of the forward and backward version of given input.\n"
              << "Creates files <s_path>.fwd.stree and <s_path>.rev.stree in the dir of <s_path>\n"
              << "Args:\n"
              << help__s_path
              << endl);
        exit(0);
    } else {
        ispec = InputSpec(input.getCmdOption("-s_path"), "");
    }
    comp(ispec);
    return 0;
}
