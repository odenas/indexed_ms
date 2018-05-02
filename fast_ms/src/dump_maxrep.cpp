#include <iostream>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>


#include "fd_ms/opt_parser.hpp"
#include "fd_ms/input_spec.hpp"
#include "fd_ms/stree_sct3.hpp"
#include "fd_ms/maxrep_vector.hpp"


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
        cerr << "-s_path <s path>" << endl;
    }
};


typedef typename StreeOhleb<>::size_type size_type;

StreeOhleb<>::size_type load_or_build(StreeOhleb<>& st, const InputSpec& ispec, const bool load){
    string potential_stree_fname = ispec.rev_cst_fname;
    using timer = std::chrono::high_resolution_clock;
    auto start = timer::now();
    if(load){
        std::cerr << " * loading the CST from " << potential_stree_fname << " ";
        sdsl::load_from_file(st, potential_stree_fname);
    } else {
        std::cerr << " * loadding reversed index string from " << ispec.s_fname << " " << endl;
        string s = ispec.load_s(true);
        std::cerr << " * building the CST of length " << s.size() << " ";
        sdsl::construct_im(st, s, 1);
    }
    auto stop = timer::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
}

void comp(const InputSpec& s_spec, const InputFlags& flags){
    StreeOhleb<> st;
    size_type load_cst_time = load_or_build(st, s_spec, flags.load_cst);
    cerr << "DONE (" << load_cst_time << " milliseconds)" << endl;

    Maxrep<StreeOhleb<>, sdsl::bit_vector> mr(st, true);

    if(flags.out_format == OutFormat::txt){
        cerr << " * dumping MAXREP to stodut" << endl;
        for(size_type i = 0; i < mr.size(); i++)
            cout << mr[i] << " ";
        cout << endl;
    } else {
        cerr << " * dumping MAXREP to " << s_spec.rev_maxrep_fname << " ";
        auto start = timer::now();
        mr.dump_vec(s_spec.rev_maxrep_fname);
        auto stop = timer::now();
        cerr << "DONE (" << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << " milliseconds)" << endl;
    }
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
