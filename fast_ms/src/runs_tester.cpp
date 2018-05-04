/* 
 * Define & test a serializable bit vector.
 */

#include <cstdlib>
#include <sdsl/int_vector.hpp>

#include "fd_ms/opt_parser.hpp"
#include "fd_ms/input_spec.hpp"
#include "fd_ms/runs_ms.hpp"
#include "fd_ms/stree_sct3.hpp"
#include "fd_ms/query.hpp"

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


class InputFlags{
public:
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
};


size_type load_or_build(cst_t& st, const InputSpec& ispec, const bool load){
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


#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))
#define STREAM_BUFFER_SIZE 1e+5


/*
* call parent(v) in sequece until reaching a node u for which wl(u, c) exists
*/
node_type parent_sequence(const cst_t& st, wl_method_t1 wl_f_ptr, const node_type& v, const char_type c) {
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

void fill_runs(const string& t_fname, const string& o_fname, const cst_t& st, wl_method_t1 wl_f_ptr){
    sdsl::int_vector_buffer<1> runs(o_fname, std::ios::out);
    Query_rev t{t_fname, (size_t) STREAM_BUFFER_SIZE};
    size_type k = t.size();
    char_type c = t[k - 1];
    node_type v = CALL_MEMBER_FN(st, wl_f_ptr)(st.root(), c),
              u = v;

    while(--k > 0){
        c = t[k-1];

        u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
        if(st.is_root(u)){
            runs[k] = 0;
            v = parent_sequence(st, wl_f_ptr, v, c);
            v = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
        } else {
            runs[k] = 1;
            v = u;
        }
    }
}

void comp(const InputSpec& s_spec, const InputFlags& flags){
    StreeOhleb<> st;
    size_type load_cst_time = load_or_build(st, s_spec, flags.load_cst);
    cerr << "DONE (" << load_cst_time << " milliseconds)" << endl;

    //fill_runs(s_spec.s_fname, runs_fname, st, &cst_t::double_rank_nofail_wl);
    bitvec_t runs;
    sdsl::load_from_file(runs, s_spec.runs_fname);
    
    sdsl::int_vector_buffer<1> q_runs(s_spec.runs_fname, std::ios::in);
    
    for(int i=0; i<runs.size(); i++){
        cout << i << " : " << runs[i] - q_runs[i] << endl;
    }
    
}


int main(int argc, char **argv){
    OptParser input(argc, argv);
    InputSpec s_spec;
    InputFlags flags;

    if(argc == 1){
        const string base_dir = {"/home/denas/projects/matching_statistics/indexed_ms/tests/datasets/testing/"};
        s_spec = InputSpec(base_dir + "mut_200s_64t_15.s");
        flags = InputFlags(false);
    } else {
        flags = InputFlags(input);
        s_spec = InputSpec(input.getCmdOption("-s_path"));
    }
    comp(s_spec, flags);
    return 0;
}

