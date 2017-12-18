#include <iostream>
#include <fstream>
#include <vector>
#include <string>


#include "fd_ms/input_spec.hpp"
#include "fd_ms/opt_parser.hpp"
#include "fd_ms/stree_sct3.hpp"
#include "fd_ms/maxrep_vector.hpp"
#include "fd_ms/runs_and_ms_algorithms.hpp"

using namespace std;
using namespace fdms;


string t, s;
StreeOhleb<> st;
MsVectors<StreeOhleb<>, sdsl::bit_vector> ms_vec;
Maxrep<StreeOhleb<>, sdsl::bit_vector> maxrep;
map<size_type, size_type> runs_stats, ms_stats;


class InputFlags{
public:
    bool load_stree;

    InputFlags(){}
    
    InputFlags(const InputFlags& f) : load_stree{f.load_stree} {}
    
    InputFlags(bool load_stree) : load_stree{load_stree} {}
    
    InputFlags (OptParser input) : load_stree{input.getCmdOption("-load_cst") == "1"} {}
};

size_type fill_runs(){
    auto runs_start = timer::now();
    size_type k = t.size();
    char_type c = t[k - 1];
    node_type v = st.double_rank_nofail_wl(st.root(), c), u = v;
    while(--k > 0){
        c = t[k-1];

        u = st.double_rank_nofail_wl(v, c);
        if(st.is_root(u)){
            ms_vec.runs[k] = 0;
            
            if(!st.has_complete_info(v))
				st.lazy_wl_followup(v);
            bool has_wl = false;
            u = st.root();
            size_type seq_len = 0;
            do{ // remove suffixes of t[k..] until you can extend by 'c'
                v = st.parent(v);
				u = st.double_rank_nofail_wl(v, c);
                has_wl = !st.is_root(u);
                seq_len += 1;
            } while(!has_wl && !st.is_root(v));
            runs_stats[seq_len] += 1;
        } else {
            ms_vec.runs[k] = 1;
        }
        v = st.double_rank_nofail_wl(v, c);
    }
    auto runs_stop = timer::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count();
}


size_type fill_ms(){
    auto runs_start = timer::now();
    size_type k = 0, h_star = k + 1, h = h_star, ms_idx = 0, ms_size = t.size();
    char_type c = t[k];
    node_type v = st.double_rank_nofail_wl(st.root(), c), u = v;

    while(k < t.size()){
        h = h_star;
        
        while(h_star < ms_size){
            c = t[h_star];
			u = st.double_rank_nofail_wl(v, c);
            if(!st.is_root(u)){
                v = u;
                h_star += 1;
            } else
                break;
        }
        ms_vec.set_next_ms_values1(0, ms_idx, h, h_star, t.size() * 2);

        if(h_star < ms_size){ // remove prefixes of t[k..h*] until you can extend by 'c'
            if(!st.has_complete_info(v))
                st.lazy_wl_followup(v);

            bool has_wl = false;
            size_type seq_len = 0;
            u = st.root();
            do{ // remove suffixes of t[k..] until you can extend by 'c'
                v = st.parent(v);
				u = st.double_rank_nofail_wl(v, c);
                has_wl = !st.is_root(u);
                seq_len += 1;
            } while(!has_wl && !st.is_root(v));
            ms_stats[seq_len] += 1;
            h_star += 1;
        }
        k = ms_vec.set_next_ms_values2(0, ms_idx, k, t.size(), t.size() * 2);
        v = u;
    }
    auto runs_stop = timer::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count();
}

void comp(const InputSpec& tspec, InputSpec& s_fwd, const string& out_path, InputFlags& flags){
	size_type t_ms = 0;

	/* load input */
    cerr << "loading input ";
    auto start = timer::now();
    t = tspec.load_s();
    s = s_fwd.load_s();
    cerr << "|s| = " << s.size() << ", |t| = " << t.size() << ". ";
    t_ms = std::chrono::duration_cast<std::chrono::milliseconds>(timer::now() - start).count();
    cerr << "DONE (" << t_ms / 1000 << " seconds)" << endl;

    /* prepare global data structures */
    ms_vec = MsVectors<StreeOhleb<>, sdsl::bit_vector>(t.size(), 1);

	/* build runs */
    cerr << "building RUNS ... " << endl;
    t_ms  = load_or_build(st, s, s_fwd.fwd_cst_fname, flags.load_stree);
    cerr << "DONE (" << t_ms / 1000 << " seconds, " << st.size() << " nodes)" << endl;
    t_ms = fill_runs();
    cerr << "DONE (" << t_ms / 1000 << " seconds)" << endl;

    /* reverse input */
    InputSpec::reverse_in_place(s);

    /* build ms */
    cerr << "building MS ... " << endl;
    t_ms = load_or_build(st, s, s_fwd.rev_cst_fname, flags.load_stree);
    cerr << "DONE (" << t_ms / 1000 << " seconds, " << st.size() << " nodes)" << endl;
    t_ms = fill_ms();
    cerr << " * total ms length : " << ms_vec.ms_size()  << " (with |t| = " << t.size() << ")" << endl;
    cerr << "DONE (" << t_ms / 1000 << " seconds)" << endl;


    cerr << "dumping reports" << endl;
    cout << "method,seq_len,cnt" << endl;
    for(auto item : runs_stats)
    	cout << "runs," << item.first << "," << item.second << endl;
    for(auto item : ms_stats)
    	cout << "ms," << item.first << "," << item.second << endl;
}


int main(int argc, char **argv){
    OptParser input(argc, argv);
    InputSpec sfwd_spec, tspec;
    InputFlags flags;
    string out_path;

    if(argc == 1){
        const string base_dir = {"/home/brt/code/matching_statistics/indexed_ms/tests/datasets/testing/"};
        tspec = InputSpec(base_dir + "rnd_200_32.t");
        sfwd_spec = InputSpec(base_dir + "rnd_200_32.s");
        out_path = "0";
        flags = InputFlags(false);
    } else {
        tspec = InputSpec(input.getCmdOption("-t_path"));
        sfwd_spec = InputSpec(input.getCmdOption("-s_path"));
        flags = InputFlags(input);
    }
    comp(tspec, sfwd_spec, out_path, flags);
    return 0;
}


