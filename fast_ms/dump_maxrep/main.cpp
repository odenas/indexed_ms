//
//  main.cpp
//  dump_maxrep
//
//  Created by denas on 5/22/17.
//  Copyright © 2017 denas. All rights reserved.
//

#include <iostream>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>


#include "opt_parser.hpp"
#include "input_spec.hpp"
#include "stree_sct3.hpp"
#include "maxrep_vector.hpp"


using namespace std;
using namespace fdms;


enum class OutFormat{txt, bin};

class InputFlags{
public:
    OutFormat out_format;
    bool load_cst;

    InputFlags(){}

    InputFlags(const OutFormat of, const bool load_cst) :
    out_format{of}, load_cst{load_cst} {}

    InputFlags(const InputFlags& i) :
    out_format{i.out_format}, load_cst{i.load_cst} {}

    InputFlags(const OptParser args){
        out_format = (args.getCmdOption("-txt_format") == "1" ? OutFormat::txt : OutFormat::bin);
        load_cst = (args.getCmdOption("-load_cst") == "1");
    }

    void help(const string exec_name) const {
        cerr << exec_name << endl << "\t";
        cerr << "[-load_cst 0/1]" << endl << "\t";
        cerr << "[-txt_format 0/1]" << endl << "\t";
        cerr << "-s_path <s path>" << endl << "\t";
        cerr << "-out_path <output path>" << endl;
    }
};


typedef typename StreeOhleb<>::size_type size_type;


void comp(const InputSpec& s_spec, const InputFlags& flags){
    auto start = timer::now();
    cerr << " * loading and reversing the string " << s_spec.s_fname << " ";
    string s = s_spec.load_s(true);
    auto stop = timer::now();
    cerr << "DONE (" << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << "seconds)" << endl;

    StreeOhleb<> st;
    size_type load_cst_time = load_or_build(st, s, s_spec.rev_cst_fname, flags.load_cst);
    cerr << "DONE (" << load_cst_time << " milliseconds)" << endl;

    Maxrep mr(st, true);

    cerr << " * dumping MAXREP to " << s_spec.rev_maxrep_fname << " ";
    start = timer::now();
    mr.dump_vec(s_spec.rev_maxrep_fname);
    stop = timer::now();
    cerr << "DONE (" << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << " milliseconds)" << endl;
}


int main(int argc, char **argv){
    OptParser input(argc, argv);
    InputSpec s_spec;
    InputFlags flags;

    if(argc == 1){
        const string base_dir = {"/Users/denas/projects/matching_statistics/indexed_ms/tests/datasets/testing/"};
        s_spec = InputSpec(base_dir + "mut_200s_64t_15.s");
        flags = InputFlags(OutFormat::txt, // lazy_wl
                           false);         // load CST
    } else {
        flags = InputFlags(input);
        s_spec = InputSpec(input.getCmdOption("-s_path"));
    }
    comp(s_spec, flags);
    return 0;
}
