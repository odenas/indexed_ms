#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <future>
#include <thread>


#include "fd_ms/opt_parser.hpp"
#include "fd_ms/input_spec.hpp"
#include "fd_ms/counter.hpp"
#include "fd_ms/query.hpp"
#include "fd_ms/slices.hpp"
#include "fd_ms/stree_sct3.hpp"
#include "fd_ms/maxrep_vector.hpp"
#include "fd_ms/p_runs_vector.hpp"
#include "fd_ms/p_ms_vector.hpp"
#include "fd_ms/ms_vector.hpp"

using namespace std;
using namespace fdms;

typedef StreeOhleb<> cst_t;
typedef typename cst_t::node_type node_type;
typedef typename cst_t::size_type size_type;
typedef typename cst_t::char_type char_type;
typedef sdsl::bit_vector bitvec_t;
typedef Maxrep<cst_t, bitvec_t> maxrep_t;
typedef Counter<size_type> counter_t;

typedef typename p_runs_vector<cst_t>::pseq_method_t pseq_method_t;
typedef typename p_runs_vector<cst_t>::wl_method_t1 wl_method_t1;
typedef typename p_runs_vector<cst_t>::wl_method_t2 wl_method_t2;
typedef typename p_runs_vector<cst_t>::pair_t pair_t;
typedef typename p_runs_vector<cst_t>::p_runs_state runs_state_t;

#define VERBOSE
//#define VVERBOSE
#define PARALLEL_POLICY std::launch::async
//#define SEQUENTIAL

cst_t st;
maxrep_t maxrep;
p_runs_vector<cst_t> runs;
p_ms_vector<cst_t> ms;
Slices<size_type> runs_slices;

vector<runs_state_t> runs_results;
vector<runs_state_t> merge_slices;
size_type available_slice_idx = 0;
std::mutex res_mutex;


class InputFlags {
private:

    void check() const {
        if (use_maxrep_rc && use_maxrep_vanilla) {
            cerr << "use_maxrep_rc and use_maxrep_vanilla cannot be active at the same time" << endl;
            exit(1);
        }
        if (use_maxrep() && !(double_rank && rank_fail)) {
            cerr << "use_maxrep_xx goes with double rank and fail" << endl;
            exit(1);
        }
        if (use_maxrep() && lazy) {
            cerr << "lazy and use_maxrep_xx cannot be active at the same time" << endl;
            cerr << "use_maxrep_xx goes with double rank and fail" << endl;
            exit(1);
        }
        if (use_maxrep() && !double_rank) {
            cerr << "single_rank and use_maxrep_xx cannot be active at the same time" << endl;
            cerr << "use_maxrep_xx goes with double rank and fail" << endl;
            exit(1);
        }
        if (use_maxrep() && !rank_fail) {
            cerr << "no_fail and use_maxrep_xx cannot be active at the same time" << endl;
            cerr << "use_maxrep_xx goes with double rank and fail" << endl;
            exit(1);
        }
        if (rank_fail && !double_rank) {
            cerr << "single_rank and rank_fail cannot be active at the same time" << endl;
            exit(1);
        }
        if (answer && avg) {
            cerr << "answer and avg cannot be active at the same time" << endl;
            exit(1);
        }
    }

public:
    bool double_rank, lazy, rank_fail, use_maxrep_vanilla, use_maxrep_rc, lca_parents;
    bool time_usage, answer, avg;
    bool load_stree, load_maxrep;
    size_t nthreads, nslices;

    InputFlags() {
    }

    InputFlags(const InputFlags& f) :
    double_rank{f.double_rank},
    lazy{f.lazy},
    rank_fail{f.rank_fail},
    use_maxrep_vanilla{f.use_maxrep_vanilla}, use_maxrep_rc{f.use_maxrep_rc},
    lca_parents{f.lca_parents},
    time_usage{f.time_usage},
    answer{f.answer}, avg{f.avg},
    load_stree{f.load_stree},
    load_maxrep{f.load_maxrep},
    nthreads{f.nthreads}, nslices{f.nslices}
    {
    }

    InputFlags(bool double_rank, bool lazy_wl, bool use_rank_fail,
            bool use_maxrep_vanilla, bool use_maxrep_rc, bool lca_parents,
            bool time_, bool ans, bool avg,
            bool load_stree, bool load_maxrep, size_t nthreads, size_t nslices) :
    double_rank{double_rank},
    lazy{lazy_wl},
    rank_fail{use_rank_fail},
    use_maxrep_vanilla{use_maxrep_vanilla}, use_maxrep_rc{use_maxrep_rc},
    lca_parents{lca_parents},
    time_usage{time_},
    answer{ans}, avg{avg},
    load_stree{load_stree},
    load_maxrep{load_maxrep},
    nthreads{nthreads}, nslices{nslices}
    {
        check();
    }

    InputFlags(OptParser input) :
    double_rank{input.getCmdOption("-double_rank") == "1"}, // use double rank
    lazy{input.getCmdOption("-lazy_wl") == "1"}, // lazy winer links
    rank_fail{input.getCmdOption("-rank_fail") == "1"}, // use the rank-and-fail strategy
    use_maxrep_rc{input.getCmdOption("-use_maxrep_rc") == "1"}, // use the maxrep vector with rank_and_check
    use_maxrep_vanilla{input.getCmdOption("-use_maxrep_vanilla") == "1"}, // use the maxrep vector the vanilla way
    lca_parents{input.getCmdOption("-lca_parents") == "1"}, // use lca insted of conscutive parent calls
    time_usage{input.getCmdOption("-time_usage") == "1"}, // time usage
    answer{input.getCmdOption("-answer") == "1"}, // answer
    avg{input.getCmdOption("-avg") == "1"}, // average matching statistics
    load_stree{input.getCmdOption("-load_cst") == "1"}, // load CST of S and S'
    load_maxrep{input.getCmdOption("-load_maxrep") == "1"}, // load MAXREP of S'
    nthreads{static_cast<size_t> (std::stoi(input.getCmdOption("-nthreads")))},
    nslices{static_cast<size_t> (std::stoi(input.getCmdOption("-nslices")))}
    { 
        nslices = (nslices > 0 ? nslices : 1);
        nthreads = (nthreads > 0 ? nthreads : 1);
        check(); 
    }

    bool use_maxrep() const {
        return (use_maxrep_rc || use_maxrep_vanilla);
    }

    wl_method_t1 get_wl_method() const {
        if (double_rank) {
            if (lazy)
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

    wl_method_t2 get_mrep_wl_method() const {
        return (use_maxrep_rc ?
                &cst_t::double_rank_fail_wl_mrep_rc :
                &cst_t::double_rank_fail_wl_mrep_vanilla);
    }

    pseq_method_t get_pseq_method() const {
        return (lca_parents ?
            &p_runs_vector<cst_t>::lca_parent :
            &p_runs_vector<cst_t>::parent_sequence);
    }

    size_type buffer_size(const InputSpec& ispec) const {
        size_type t_length = ispec.t_size();
        return (t_length > 1000 ? t_length / 100 : t_length);
    }

};

int fill_runs_slice_thread1(const InputSpec& ispec) {
    while(true){
        size_type thread_id;
        {
#ifndef SEQUENTIAL
            std::lock_guard<std::mutex> l(res_mutex);
#endif
            if(available_slice_idx < runs_slices.nslices){
                thread_id = available_slice_idx++;
#ifdef VERBOSE
                cerr << " *** " << runs_slices.repr(thread_id) << endl;
#endif
            }
            else{
                return 0;
            }
        }
        try{
            runs_results[thread_id] = runs.fill_slice(ispec, st, thread_id);
        } catch (string s) {
            throw string{"runs.fill_slice on slice " + to_string(thread_id) + 
                " failed with message: " + s};
        }
    }
    return 1;
}

int fill_runs_slice_thread2(const InputSpec& ispec, const int i){
    while(true){
        size_type thread_id;
        {
#ifndef SEQUENTIAL
            std::lock_guard<std::mutex> l(res_mutex);
#endif
            if(available_slice_idx < merge_slices.size()){
                thread_id = available_slice_idx++;
#ifdef VERBOSE
                runs_state_t st = merge_slices[thread_id];
                cerr << " *** [" << i << "]"
                     << st.repr()
                     << "     intervals "
                     << runs.m_slices.slice_idx(st.ff_index)
                     << " - "
                     << runs.m_slices.slice_idx(st.lf_index)
                     << endl;
#endif
            }
            else {
                return 0;
            }
        } // lock released at the end of block
        runs.fill_inter_slice(ispec, st, merge_slices[thread_id]);
    }
    return 1;
}

void build_runs(const InputSpec& ispec, counter_t& time_usage, InputFlags& flags) {
    if(flags.lazy)
        throw string{"lazy mode not supported"};

    cerr << "building RUNS ... " << endl;

    time_usage.reg["runs_cst"] = cst_t::load_or_build(st, ispec, false, flags.load_stree);
    (cerr << "DONE (" << time_usage.reg["runs_cst"] / 1000 << " seconds, "
            << st.size() << " leaves)" << endl);

    runs_slices = Slices<size_type>(ispec.t_size(), flags.nslices);
    runs_results = vector<runs_state_t>(runs_slices.nslices);
    runs = p_runs_vector<cst_t>(1024, runs_slices, 
                                flags.get_wl_method(), flags.get_pseq_method());

    available_slice_idx = 0;
    auto runs_start = timer::now();
    {
        (cerr << " ** filling " << runs_slices.nslices << " slices with : " 
                << flags.nthreads  << " threads ..." << endl);

#ifdef SEQUENTIAL
        std::vector<int> thread_st(flags.nthreads);
        for (size_type i = 0; i < flags.nthreads; i++) {
            thread_st[i] = fill_runs_slice_thread1(ispec);
        }
#else
        std::vector<std::future<int>> thread_st(flags.nthreads);
        for (size_type i = 0; i < flags.nthreads; i++) {
            thread_st[i] = std::async(PARALLEL_POLICY, fill_runs_slice_thread1, ispec);
        }
#endif

        int sum = 0;
        for (size_type i = 0; i < flags.nthreads; i++) {
#ifdef SEQUENTIAL
            sum += thread_st[i]; 
#else
            sum += thread_st[i].get();
#endif
        }
        cerr << " ** DONE: " << sum << endl;
    }
    time_usage.register_now("runs_build", runs_start);

    available_slice_idx = 0;
    merge_slices = runs.reduce(runs_results);
    runs_start = timer::now();
    {
        (cerr << " * correcting " << merge_slices.size() << " intervals over " 
                << flags.nthreads << " threads ... " << endl);
#ifdef SEQUENTIAL
        std::vector<int>thread_st(flags.nthreads);
        for (int i = 0; i < (int) flags.nthreads; i++) {
            thread_st[i] = fill_runs_slice_thread2(ispec, i);
        }
#else
        std::vector<std::future<int>> thread_st(flags.nthreads);
        for (int i = 0; i < (int) flags.nthreads; i++) {
            thread_st[i] = std::async(PARALLEL_POLICY, fill_runs_slice_thread2, ispec, i);
        }
#endif

        int sum = 0;
        for (size_type i = 0; i < flags.nthreads; i++) {
#ifdef SEQUENTIAL
            sum += thread_st[i]; 
#else
            sum += thread_st[i].get();
#endif
        }
        cerr << " ** DONE: " << sum << endl;
    }
    time_usage.register_now("runs_correct", runs_start);

    runs_start = timer::now();
    cerr << " * merging into " << ispec.runs_fname << " ... " << endl;
    runs.merge(ispec, merge_slices);
    time_usage.register_now("runs_merge", runs_start);
}

int fill_ms_slice_thread(const size_type thread_id, const pair_t slice, const InputSpec& ispec){
    return ms.fill_slice(ispec, st, slice, thread_id);
}

void build_ms(const InputSpec& ispec, counter_t& time_usage, const InputFlags& flags) {
    cerr << "building MS ... " << endl;
    /* build the CST */
    time_usage.reg["ms_cst"] = cst_t::load_or_build(st, ispec, true, flags.load_stree);
    cerr << "DONE (" << time_usage.reg["ms_cst"] / 1000 << " seconds, " << st.size() << " leaves)" << endl;
    Slices<size_type> slices(Query::query_length(ispec.t_fname), flags.nthreads);
    ms = p_ms_vector<cst_t>(flags.nthreads, flags.buffer_size(ispec),
            flags.get_wl_method(), flags.get_pseq_method());

    auto ms_start = timer::now();
    /* compute MS */
    {
#ifdef SEQUENTIAL
        std::vector<int> results(flags.nthreads);
#else
        std::vector<std::future<int>> results(flags.nthreads);
#endif
        for(size_type i=0; i<flags.nthreads; i++){
            cerr << " ** launching ms computation over : " << slices.repr(i) << endl;
#ifdef SEQUENTIAL
            results[i] = fill_ms_slice_thread(i, slices[i], ispec);
#else
            results[i] = std::async(PARALLEL_POLICY, fill_ms_slice_thread, i, slices[i], ispec);
#endif
        }
        for(size_type i=0; i<flags.nthreads; i++){
            double ms_max = slices.input_size * 2;
#ifdef SEQUENTIAL
            size_type filled = results[i];
#else
            size_type filled = results[i].get();
#endif
            assert(filled <= ms_max);
            (cerr << " *** [" << i << "]" << "filled " << filled <<
                    " of " << ms_max << " entries " << endl);
        }
    }
    time_usage.register_now("ms_build", ms_start);

    cerr << " * merging into " << ispec.ms_fname << " ... " << endl;
    ms_start = timer::now();
    ms.merge(ispec, slices);
    time_usage.register_now("ms_merge", ms_start);
}

void comp(const InputSpec& ispec, counter_t& time_usage, InputFlags& flags) {
    auto comp_start = timer::now();
    try{
        build_runs(ispec, time_usage, flags);
    } catch(string s) {
    	cerr << "ERROR from buiold_runs: " << s << endl;
        throw string{"build_runs failed with message: \n" + s};
    }
    time_usage.register_now("runs_total", comp_start, true);
#ifdef VVERBOSE
    p_runs_vector<cst_t>::show(ispec.runs_fname, cerr);
#endif

    comp_start = timer::now();
    build_ms(ispec, time_usage, flags);
    time_usage.register_now("ms_total", comp_start, true);
#ifdef VVERBOSE
    p_runs_vector<cst_t>::show(ispec.ms_fname, cerr);
#endif

    comp_start = timer::now();
    if (flags.answer) {
        ms_vector<cst_t>::show_MS(ispec, cout);
        cout << endl;
    } else if (flags.avg)
        cout << ms_vector<cst_t>::avg_matching_statistics(ispec) << endl;
    time_usage.register_now("answer", comp_start);
}

int main(int argc, char **argv) {
    OptParser input(argc, argv);
    InputSpec ispec;
    InputFlags flags;
    counter_t time_usage{};

    if (argc == 1) {
        const string base_dir = {"/home/brt/code/matching_statistics/indexed_ms/fast_ms/tests/"};
        ispec = InputSpec(base_dir + "a.s", base_dir + "a.t");
        flags = InputFlags(true, // use double rank
                false, // lazy_wl
                false, // rank-and-fail
                false, // use maxrep vanilla
                false, // use maxrep rank&check
                false, // lca_parents
                false, // time
                true, // ans
                false, // avg
                false, // load CST
                false, // load MAXREP
                2, // nthreads
                2 // nslices
                );
    } else {
        ispec = InputSpec(input.getCmdOption("-s_path"), input.getCmdOption("-t_path"));
        flags = InputFlags(input);
    }
    
    if(ispec.t_size() < flags.nthreads * 2){
        (cerr << "max parallelization is " << ispec.t_size() / 2 
                << "but nthreads = " << flags.nthreads << endl);
        return 1;
    }

    auto comp_start = timer::now();
    comp(ispec, time_usage, flags);
    time_usage.register_now("comp_total", comp_start);

    if (flags.time_usage) {
        cerr << "dumping reports ..." << endl;
        cout << "len_s,len_t,item,value" << endl;
        for (auto item : time_usage.reg)
            (cout << st.size() - 1 << ","
                    << ms_vector<cst_t>::size(ispec) << ","
                    << item.first << ","
                    << item.second << endl);
    }

    return 0;
}

