/*
 * Define & test a serializable bit vector.
 */

#include <cstdlib>
#include <future>
#include <thread>

#include <sdsl/int_vector.hpp>

#include "fd_ms/opt_parser.hpp"
#include "fd_ms/input_spec.hpp"
#include "fd_ms/runs_ms.hpp"
#include "fd_ms/stree_sct3.hpp"
#include "fd_ms/query.hpp"
#include "fd_ms/counter.hpp"

#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))

using namespace std;
using namespace fdms;

typedef sdsl::bit_vector            bitvec_t;
typedef StreeOhleb<> cst_t;
typedef typename cst_t::size_type size_type;
typedef typename cst_t::node_type node_type;
typedef typename cst_t::char_type char_type;
typedef MsVectors<cst_t, bitvec_t>  msvec_t;

//Declare various wl strategies
typedef node_type (cst_t::*wl_method_t1) (const node_type& v, const char_type c) const;
typedef node_type (cst_t::*wl_method_t2) (const node_type& v, const char_type c, const bool is_max) const;
// Declare a parent-sequence strategy type
typedef node_type (MsVectors<cst_t, bitvec_t>::*pseq_method_t)(const cst_t& st, wl_method_t1 wl_f_ptr, const node_type& v, const char_type c) const;
typedef tuple<size_type, size_type, node_type> runs_rt;
typedef std::pair<size_type, size_type>        Interval;
typedef Counter<size_type>                     counter_t;
cst_t st;


class InputFlags{
public:
    bool double_rank, rank_fail;
    bool time_usage;
    bool load_cst;
    size_t nthreads;

    InputFlags(){}

    InputFlags(const bool load_cst, const size_t nthreads) :
    load_cst{load_cst}, nthreads{nthreads}
    {}

    InputFlags(const InputFlags& i) :
    load_cst{i.load_cst}, nthreads{i.nthreads}
    {}

    InputFlags(const OptParser args) :
    load_cst{args.getCmdOption("-load_cst") == "1"},
    nthreads{static_cast<size_t>(std::stoi(args.getCmdOption("-nthreads")))}
    {}

    wl_method_t1 get_wl_method() const {
        if(double_rank)
            return (rank_fail ? &cst_t::double_rank_fail_wl : &cst_t::double_rank_nofail_wl);
        return &cst_t::single_rank_wl;
    }

};


size_type load_or_build(const InputSpec& ispec, const bool load){
    string potential_stree_fname = ispec.fwd_cst_fname;
    using timer = std::chrono::high_resolution_clock;
    auto start = timer::now();
    if(load){
        std::cerr << " * loading the CST from " << potential_stree_fname << " ";
        sdsl::load_from_file(st, potential_stree_fname);
    } else {
        std::cerr << " * loadding  index string from " << ispec.s_fname << " " << endl;
        string s = ispec.load_s(false);
        std::cerr << " * building the CST of length " << s.size() << " ";
        sdsl::construct_im(st, s, 1);
    }
    auto stop = timer::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
}


/*
* call parent(v) in sequece until reaching a node u for which wl(u, c) exists
*/
node_type parent_sequence(wl_method_t1 wl_f_ptr, const node_type& v, const char_type c) {
    node_type vv = v, u = st.root();

    if(!st.has_complete_info(vv))
        st.lazy_wl_followup(vv);

    bool has_wl = false;
    do{ // remove suffixes of t[k..] until you can extend by 'c'
        vv = st.parent(vv);
        u = CALL_MEMBER_FN(st, wl_f_ptr)(vv, c);
        has_wl = !st.is_root(u);
    } while(!has_wl && !st.is_root(vv));
    return vv;
}

runs_rt dump_runs_slice(const string& t_fname, const string& o_fname,
        const Interval slice, node_type v, wl_method_t1 wl_f_ptr){

    Query_rev t{t_fname, (size_t) 1024*1024};
    size_type first_fail = 0, last_fail = 0, from = slice.first, to = slice.second;
    node_type last_fail_node = v, u = v;
    bool lca_parents = false;  // TODO: changeme

    size_type k = to;
    char_type c = t[k - 1];
    bool idx_set = false;

    sdsl::int_vector_buffer<1> runs(o_fname, std::ios::out);
    while(--k > from){
        assert (k > from && k < to);
        c = t[k-1];

        u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
        if(st.is_root(u)){
            if(!idx_set){ // first failing wl()
                first_fail = k;
                idx_set = true;
            }
            runs[k] = 0;
            // remove suffixes of t[k..] until you can extend by 'c'
            if(lca_parents){
#ifdef VERBOSE
                cerr << "runs_lca" << endl;
#endif
                v = st.maxrep_ancestor(v, c);
                u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
            } else {
#ifdef VERBOSE
                cerr << "runs_pseq" << endl;
#endif
                if(!st.has_complete_info(v))
                    st.lazy_wl_followup(v);
                bool has_wl = false;
                u = st.root();
                do{ // remove suffixes of t[k..] until you can extend by 'c'
                    v = st.parent(v);
                    u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
                    has_wl = !st.is_root(u);
                } while(!has_wl && !st.is_root(v));
            }

            // idx of last 0 in runs - 1 (within this block) and corresponding wl(node)
            last_fail_node = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);// given, parent_sequence() above, this has a wl()
            last_fail = k;
        } else {
            runs[k] = 1;
        }
        v = CALL_MEMBER_FN(st, wl_f_ptr)(v, c); // update v
    }
    if(!idx_set){
        first_fail = last_fail = from + 1;
        last_fail_node = v;
    }
    return make_tuple(first_fail, last_fail, last_fail_node);
}


runs_rt fill_runs_slice_thread(const Interval slice, const node_type v,
        const InputFlags& flags, const InputSpec& tspec){
    return dump_runs_slice(tspec.s_fname, tspec.runs_fname,
            slice, v, flags.get_wl_method());
}

// this should go to Slices
vector<runs_rt> aa(const vector<runs_rt> v, const Slices<size_type> slices){
    vector<runs_rt> u;
    u.reserve(v.size());

    int i = 0, j = 1, n = v.size();
    while(j < n){
        runs_rt prev_slice = v[i], next_slice = v[j];
        while(get<0>(next_slice) == get<0>(slices.slices[j]) + 1){ // j-th slice had a full match
            j++;
        }
        next_slice = v[j - (j == n)];
        u.push_back(make_tuple(get<0>(prev_slice), get<1>(next_slice), get<2>(next_slice)));
        i = j;
        j += 1;
    }
    return u;

}

void dump_runs(const InputFlags& flags, const InputSpec& sspec, const InputSpec& tspec, counter_t& time_usage){
    cerr << "building RUNS ... " << endl;

    /* build the CST */
    time_usage.reg["runs_cst"]  = load_or_build(sspec, flags.load_cst);
    cerr << "DONE (" << time_usage.reg["runs_cst"] / 1000 << " seconds, " << st.size() << " leaves)" << endl;

    /* compute RUNS */
    auto runs_start = timer::now();
    std::vector<std::future<runs_rt>> results(flags.nthreads);
    Slices<size_type> slices(Query::query_length(tspec.s_fname), flags.nthreads);

    /* open connection to the query string */
    Query_rev t{tspec.s_fname, (size_t) 1e+5};
    for(size_type i=0; i<flags.nthreads; i++){
        node_type v = st.double_rank_nofail_wl(st.root(), t[slices[i].second - 1]); // stree node
        cerr << " ** launching runs computation over : " << slices.repr(i) << " ";
        cerr << "(" << v.i << ", " << v.j << ")" << endl;
        results[i] = std::async(std::launch::async, fill_runs_slice_thread, slices[i], v, flags, tspec);
    }
    vector<runs_rt> runs_results(flags.nthreads);
    for(size_type i=0; i<flags.nthreads; i++){
        runs_results[i] = results[i].get();
        cerr << " *** [" << get<0>(runs_results[i]) << " .. " << get<1>(runs_results[i]) << ")" << endl;
    }
    vector<runs_rt> merge_idx = aa(runs_results, slices);
    //ms_vec.show_runs(cerr);

    cerr << " * merging over " << merge_idx.size() << " threads ... " << endl;
    for(int i = 0; i < (int) merge_idx.size(); i++){
        cerr << " ** ([" << get<0>(merge_idx[i]) << ", " << get<1>(merge_idx[i]) << "), " << "(" << get<2>(merge_idx[i]).i << ", " << get<2>(merge_idx[i]).j << ")) " << endl;
        results[i] = std::async(std::launch::async, fill_runs_slice_thread,
                make_pair(get<0>(merge_idx[i]) - 1, get<1>(merge_idx[i])),
                get<2>(merge_idx[i]),
                flags, tspec);
    }
    for(int i = 0; i < merge_idx.size(); i++)
        results[i].get();
    //ms_vec.show_runs(cerr);
    time_usage.register_now("runs_bvector", runs_start);
}


void comp(const InputSpec& t_spec, const InputSpec& s_spec, const InputFlags& flags){
    /* build runs */
    counter_t time_usage{};
    dump_runs(flags, s_spec, t_spec, time_usage);

    if(flags.time_usage){
        cerr << "dumping reports" << endl;
        cout << "len_s,len_t,item,value" << endl;
        if(flags.time_usage){
            for(auto item : time_usage.reg)
                cout << st.size() - 1 << "," << item.first << "," << item.second << endl;
        }
    }
}

int main(int argc, char **argv){
    OptParser input(argc, argv);
    InputSpec s_spec, t_spec;
    InputFlags flags;

    if(argc == 1){
        const string base_dir = {"/home/denas/projects/matching_statistics/indexed_ms/tests/datasets/testing/"};
        s_spec = InputSpec(base_dir + "mut_200s_64t_15.s");
        t_spec = InputSpec(base_dir + "mut_200s_64t_15.t");
        flags = InputFlags(false, 1);
    } else {
        s_spec = InputSpec(input.getCmdOption("-s_path"));
        t_spec = InputSpec(input.getCmdOption("-t_path"));
        flags = InputFlags(input);
    }
    comp(t_spec, s_spec, flags);
    return 0;
}

