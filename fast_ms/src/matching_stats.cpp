#include <iostream>
#include <fstream>
#include <vector>
#include <string>

//#define VERBOSE

#include "fd_ms/opt_parser.hpp"
#include "fd_ms/input_spec.hpp"
#include "fd_ms/counter.hpp"
#include "fd_ms/query.hpp"
#include "fd_ms/counter.hpp"
#include "fd_ms/stree_sct3.hpp"
#include "fd_ms/maxrep_vector.hpp"
#include "fd_ms/runs_ms.hpp"


using namespace std;
using namespace fdms;

typedef StreeOhleb<>                cst_t;
typedef typename cst_t::node_type   node_type;
typedef typename cst_t::size_type   size_type;
typedef typename cst_t::char_type   char_type;
typedef sdsl::bit_vector            bitvec_t;
typedef MsVectors<cst_t, bitvec_t>  msvec_t;
typedef Maxrep<cst_t, bitvec_t>     maxrep_t;
typedef Counter<size_type>          counter_t;


cst_t st;
msvec_t ms_vec;
maxrep_t maxrep;

typedef node_type (cst_t::*wl_method_t1) (const node_type& v, const char_type c) const;
typedef node_type (cst_t::*wl_method_t2) (const node_type& v, const char_type c, const bool is_max) const;
typedef node_type (msvec_t::*pseq_method_t) (const cst_t& st, wl_method_t1 wl_f_ptr, const node_type& v, const char_type c) const;

class InputFlags{
private:
    void check() const {
        if(use_maxrep_rc && use_maxrep_vanilla){
            cerr << "use_maxrep_rc and use_maxrep_vanilla cannot be active at the same time" << endl;
            exit(1);
        }
        if(use_maxrep() && !(double_rank && rank_fail)){
            cerr << "use_maxrep_xx goes with double rank and fail" << endl;
            exit(1);
        }
        if (use_maxrep() && lazy){
            cerr << "lazy and use_maxrep_xx cannot be active at the same time" << endl;
            cerr << "use_maxrep_xx goes with double rank and fail" << endl;
            exit(1);
        }
        if (use_maxrep() && !double_rank){
            cerr << "single_rank and use_maxrep_xx cannot be active at the same time" << endl;
            cerr << "use_maxrep_xx goes with double rank and fail" << endl;
            exit(1);
        }
        if (use_maxrep() && !rank_fail){
            cerr << "no_fail and use_maxrep_xx cannot be active at the same time" << endl;
            cerr << "use_maxrep_xx goes with double rank and fail" << endl;
            exit(1);
        }
        if (rank_fail && !double_rank){
            cerr << "single_rank and rank_fail cannot be active at the same time" << endl;
            exit(1);
        }
        if (answer && avg){
            cerr << "answer and avg cannot be active at the same time" << endl;
            exit(1);
        }
    }

public:
    bool double_rank, lazy, rank_fail, use_maxrep_vanilla, use_maxrep_rc, lca_parents;
    bool time_usage, answer, avg;
    bool load_stree, load_maxrep;

    InputFlags(){}

    InputFlags(const InputFlags& f) :
        double_rank{f.double_rank},
		lazy{f.lazy},
		rank_fail{f.rank_fail},
		use_maxrep_vanilla{f.use_maxrep_vanilla}, use_maxrep_rc{f.use_maxrep_rc},
		lca_parents{f.lca_parents},
		time_usage{f.time_usage},
		answer{f.answer}, avg{f.avg},
		load_stree{f.load_stree},
		load_maxrep{f.load_maxrep}{}

    InputFlags(bool double_rank, bool lazy_wl, bool use_rank_fail, bool use_maxrep_vanilla, bool use_maxrep_rc, bool lca_parents,
               bool time_, bool ans, bool avg,
               bool load_stree, bool load_maxrep) :
        double_rank{double_rank},
		lazy{lazy_wl},
		rank_fail{use_rank_fail},
		use_maxrep_vanilla{use_maxrep_vanilla}, use_maxrep_rc{use_maxrep_rc},
		lca_parents{lca_parents},
		time_usage {time_},
		answer {ans}, avg{avg},
		load_stree{load_stree},
		load_maxrep{load_maxrep}
	{ check(); }

    InputFlags (OptParser input) :
		double_rank {input.getCmdOption("-double_rank") == "1"},  // use double rank
		lazy {input.getCmdOption("-lazy_wl") == "1"},             // lazy winer links
		rank_fail {input.getCmdOption("-rank_fail") == "1"},      // use the rank-and-fail strategy
		use_maxrep_rc {input.getCmdOption("-use_maxrep_rc") == "1"},    // use the maxrep vector with rank_and_check
		use_maxrep_vanilla {input.getCmdOption("-use_maxrep_vanilla") == "1"},    // use the maxrep vector the vanilla way
		lca_parents {input.getCmdOption("-lca_parents") == "1"},  // use lca insted of conscutive parent calls
		time_usage {input.getCmdOption("-time_usage") == "1"},    // time usage
		answer {input.getCmdOption("-answer") == "1"},            // answer
		avg {input.getCmdOption("-avg") == "1"},                  // average matching statistics
		load_stree{input.getCmdOption("-load_cst") == "1"},       // load CST of S and S'
		load_maxrep{input.getCmdOption("-load_maxrep") == "1"}    // load MAXREP of S'
    { check(); }

    bool use_maxrep() const { return (use_maxrep_rc || use_maxrep_vanilla); }

	wl_method_t1 get_wl_method() {
        if(double_rank){
    	    if(lazy)
	            return (rank_fail ?
	        		    &cst_t::lazy_double_rank_fail_wl :
	        		    &cst_t::lazy_double_rank_nofail_wl);
            return (rank_fail ?
                    &cst_t::double_rank_fail_wl :
                    &cst_t::double_rank_nofail_wl);
        } else {
            return (lazy ?
        		    &cst_t::lazy_single_rank_wl :
        		    &cst_t::single_rank_wl);
        }
	}

	wl_method_t2 get_mrep_wl_method() {
        return (use_maxrep_rc ?
        		&cst_t::double_rank_fail_wl_mrep_rc:
        		&cst_t::double_rank_fail_wl_mrep_vanilla);
	}

	pseq_method_t get_pseq_method() {
		return (lca_parents ? &msvec_t::lca_parent : &msvec_t::parent_sequence);
	}
};

StreeOhleb<>::size_type load_or_build(StreeOhleb<>& st, const InputSpec& ispec, const bool reverse, const bool load){
    string potential_stree_fname = (reverse ? ispec.rev_cst_fname : ispec.fwd_cst_fname);
    using timer = std::chrono::high_resolution_clock;
    auto start = timer::now();
    if(load){
        std::cerr << " * loading the CST from " << potential_stree_fname << " ";
        sdsl::load_from_file(st, potential_stree_fname);
    } else {
        std::cerr << " * loadding index string  from " << ispec.s_fname << " " << endl;
        string s = ispec.load_s(reverse);
        std::cerr << " * building the CST of length " << s.size() << " ";
        sdsl::construct_im(st, s, 1);
    }
    auto stop = timer::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
}


void comp(const InputSpec& tspec, InputSpec& s_fwd, counter_t& time_usage, InputFlags& flags){
    size_type t_length = Query::query_length(tspec.s_fname);

    /* prepare global data structures */
    ms_vec = msvec_t(t_length);

    /* build runs */
    cerr << "building RUNS ... " << endl;
    /* build the CST */
    time_usage.reg["runs_cst"]  = load_or_build(st, s_fwd, false, flags.load_stree);
    cerr << "DONE (" << time_usage.reg["runs_cst"] / 1000 << " seconds, " << st.size() << " leaves)" << endl;
    /* compute RUNS */
    if(flags.use_maxrep()){
		/* build the maxrep vector */
        time_usage.reg["ms_maxrep"] = maxrep_t::load_or_build(maxrep, st, s_fwd.fwd_maxrep_fname, flags.load_maxrep);
        cerr << "DONE (" << time_usage.reg["ms_maxrep"] / 1000 << " seconds)" << endl;
        auto start = timer::now();
		ms_vec.fill_runs(tspec.s_fname, st, flags.get_mrep_wl_method(), maxrep);
		time_usage.register_now("runs_bvector", start);
	} else {
		auto start = timer::now();
		ms_vec.fill_runs(tspec.s_fname, st, flags.get_wl_method(), flags.get_pseq_method());
		time_usage.register_now("runs_bvector", start);
	}

    /* build ms */
    cerr << "building MS ... " << endl;
    /* build the CST */
    time_usage.reg["ms_cst"] = load_or_build(st, s_fwd, true, flags.load_stree);
    cerr << "DONE (" << time_usage.reg["ms_cst"] / 1000 << " seconds, " << st.size() << " leaves)" << endl;
    /* compute MS */
    if(flags.use_maxrep()){
		/* build the maxrep vector */
        time_usage.reg["ms_maxrep"] = maxrep_t::load_or_build(maxrep, st, s_fwd.rev_maxrep_fname, flags.load_maxrep);
        cerr << "DONE (" << time_usage.reg["ms_maxrep"] / 1000 << " seconds)" << endl;
        auto start = timer::now();
        ms_vec.fill_ms(tspec.s_fname, st, flags.get_mrep_wl_method(), maxrep);
		time_usage.register_now("ms_bvector", start);
    } else {
        auto start = timer::now();
        ms_vec.fill_ms(tspec.s_fname, st, flags.get_wl_method(), flags.get_pseq_method());
		time_usage.register_now("ms_bvector", start);
    }
    cerr << " * total ms length : " << ms_vec.ms.size()  << " (with |t| = " << t_length << ")" << endl;

    if(flags.time_usage){
        cerr << "dumping reports ..." << endl;
        cout << "len_s,len_t,item,value" << endl;
        for(auto item : time_usage.reg)
            cout << st.size() - 1<< "," << ms_vec.runs.size() << "," << item.first << "," << item.second << endl;
    }

    if(flags.answer){
        ms_vec.show_MS(cout);
        cout << endl;
    }
    else if(flags.avg)
        cout << ms_vec.avg_matching_statistics() << endl;
}


int main(int argc, char **argv){
    OptParser input(argc, argv);
    InputSpec sfwd_spec, tspec;
    InputFlags flags;
    counter_t time_usage{};
 
    if(argc == 1){
        const string base_dir = {"/home/brt/code/matching_statistics/indexed_ms/fast_ms/tests/"};
        tspec = InputSpec(base_dir + "a.t");
        sfwd_spec = InputSpec(base_dir + "a.s");
        flags = InputFlags(true,  // use double rank
                           false, // lazy_wl
                           true, // rank-and-fail
                           false, // use maxrep vanilla
                           true, // use maxrep rank&check
                           false, // lca_parents
                           false, // time
                           false,  // ans
                           false, // avg
                           false, // load CST
                           false  // load MAXREP
                           );
    } else {
        tspec = InputSpec(input.getCmdOption("-t_path"));
        sfwd_spec = InputSpec(input.getCmdOption("-s_path"));
        flags = InputFlags(input);
    }
    comp(tspec, sfwd_spec, time_usage, flags);
    return 0;
}

