#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <future>
#include <thread>


#include "fd_ms/input_spec.hpp"
#include "fd_ms/opt_parser.hpp"
#include "fd_ms/stree_sct3.hpp"
#include "fd_ms/maxrep_vector.hpp"
#include "fd_ms/runs_and_ms_algorithms.hpp"
#include "fd_ms/slices.hpp"

using namespace std;
using namespace fdms;


string t, s;
StreeOhleb<> st;
MsVectors<StreeOhleb<>, sdsl::bit_vector> ms_vec;
Maxrep<StreeOhleb<>, sdsl::bit_vector> maxrep;
Counter time_usage;


class InputFlags{
private:
    void check() const { }

public:
    bool load_stree, load_maxrep;
    size_t nthreads;
    
    InputFlags(){}
    
    InputFlags(const InputFlags& f) : load_stree{f.load_stree}, load_maxrep{f.load_maxrep} {}
    
    InputFlags(bool load_stree, bool load_maxrep) : load_stree{load_stree}, load_maxrep{load_maxrep}
	{ check(); }
    
    InputFlags (OptParser input) :
		load_stree{input.getCmdOption("-load_cst") == "1"},       // load CST of S and S'
		load_maxrep{input.getCmdOption("-load_maxrep") == "1"}    // load MAXREP of S'
    { check(); }

    wl_method_t1 wl_method() const { return &StreeOhleb<>::single_rank_wl; }
};


void build_runs_ohleb(const InputFlags& flags, const InputSpec &s_fwd){
    cerr << "building RUNS ... " << endl;
	Interval slice = std::make_pair(0, t.size());

    /* build the CST */
    time_usage["runs_cst"]  = load_or_build(st, s, s_fwd.fwd_cst_fname, flags.load_stree);
    cerr << "DONE (" << time_usage["runs_cst"] / 1000 << " seconds, " << st.size() << " nodes)" << endl;

    /* compute RUNS */
    auto runs_start = timer::now();
    fill_runs_slice(t, st, flags.wl_method(), false, ms_vec,
                    st.double_rank_nofail_wl(st.root(), t[slice.second - 1]), slice);
    auto runs_stop = timer::now();
    time_usage["runs_bvector"]  = std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count();
    cerr << "DONE (" << time_usage["runs_bvector"] / 1000 << " seconds)" << endl;
}

void build_ms_ohleb(const InputFlags& flags, InputSpec &s_fwd){
    cerr << "building MS ... " << endl;
	Interval slice = std::make_pair(0, t.size());

    /* build the CST */
    time_usage["ms_cst"] = load_or_build(st, s, s_fwd.rev_cst_fname, flags.load_stree);
    cerr << "DONE (" << time_usage["ms_cst"] / 1000 << " seconds, " << st.size() << " nodes)" << endl;

    /* build MS */
    auto runs_start = timer::now();
    fill_ms_slice(t, st, flags.wl_method(), false, ms_vec, 0, slice); 
    auto runs_stop = timer::now();
    time_usage["ms_bvector"] = std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count();

    cerr << " * total ms length : " << ms_vec.ms_size()  << " (with |t| = " << t.size() << ")" << endl;
    cerr << "DONE (" << time_usage["ms_bvector"] / 1000 << " seconds)" << endl;
}

void comp(const InputSpec& tspec, InputSpec& S_fwd, const string& out_path, InputFlags& flags){
    auto comp_start = timer::now();
    cerr << "loading input ";
    auto start = timer::now();
    t = tspec.load_s();
    cerr << ". ";
    s = S_fwd.load_s();
    cerr << ". ";
    cerr << "|s| = " << s.size() << ", |t| = " << t.size() << ". ";
    time_usage["loadstr"] = std::chrono::duration_cast<std::chrono::milliseconds>(timer::now() - start).count();
    cerr << "DONE (" << time_usage["loadstr"] / 1000 << " seconds)" << endl;

    /* prepare global data structures */
    ms_vec = MsVectors<StreeOhleb<>, sdsl::bit_vector>(t.size(), 1);

    build_runs_ohleb(flags, S_fwd);

    start = timer::now();
    cerr << " * reversing string s of length " << s.size() << " ";
    InputSpec::reverse_in_place(s);
    time_usage["reverse_str"] = std::chrono::duration_cast<std::chrono::milliseconds>(timer::now() - start).count();
    cerr << "DONE (" << time_usage["reverse_str"] / 1000 << " seconds)" << endl;

    build_ms_ohleb(flags, S_fwd);
    time_usage["total_time"] = std::chrono::duration_cast<std::chrono::milliseconds>(timer::now() - comp_start).count();

    cerr << "dumping reports" << endl;
    cout << "len_s,len_t,item,value" << endl;
}


int main(int argc, char **argv){
    OptParser input(argc, argv);
    InputSpec sfwd_spec, tspec;
    InputFlags flags;
    string out_path;

    if(argc == 1){
        const string base_dir = {"/home/brt/Documents/projects/matching_statistics/indexed_ms/tests/datasets/testing/"};
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


