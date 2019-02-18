#include <iostream>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>


#include "fd_ms/opt_parser.hpp"
#include "fd_ms/input_spec.hpp"
#include "fd_ms/stree_sct3.hpp"
#include "fd_ms/maxrep_vector.hpp"
#include "fd_ms/help.hpp"


using namespace std;
using namespace fdms;

typedef StreeOhleb<>              cst_t;
typedef typename cst_t::size_type size_type;

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

void comp(const InputSpec& s_spec, const InputFlags& flags){
    StreeOhleb<> st;
    size_type load_cst_time = cst_t::load_or_build(st, s_spec, true, flags.load_cst);
    cerr << "DONE (" << load_cst_time/1000 << " seconds)" << endl;

    Maxrep<cst_t, sdsl::bit_vector> mr(st, true);

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
        cerr << "DONE (" << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << " seconds)" << endl;
    }
}


int main(int argc, char **argv){
    OptParser input(argc, argv);
    InputSpec s_spec;
    InputFlags flags;

    if(argc == 1){
        const string base_dir = {"/Users/denas/projects/matching_statistics/indexed_ms/tests/datasets/testing/"};
        s_spec = InputSpec(base_dir + "mut_200s_64t_15.s", "");
        flags = InputFlags(OutFormat::txt, // lazy_wl
                           false);         // load CST
        (cerr << "Dump the cst of the forward and backward version of given input.\n"
              << "Creates file <s_path>.rev.maxrep in the dir of <s_path>\n"
              << "Args:\n"
              << help__s_path
              << help__load_cst
              << endl);
        exit(0);
    } else {
        flags = InputFlags(input);
        s_spec = InputSpec(input.getCmdOption("-s_path"), "'");
    }
    comp(s_spec, flags);
    return 0;
}
