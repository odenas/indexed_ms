#include <iostream>
#include <string>

#include <sdsl/suffix_trees.hpp>

#include "fd_ms/opt_parser.hpp"
#include "fd_ms/input_spec.hpp"
#include "fd_ms/stree_sct3.hpp"
#include "fd_ms/runs_ms.hpp"
#include "fd_ms/maxrep_vector.hpp"


using namespace std;
using namespace fdms;

typedef typename StreeOhleb<>::node_type node_type;
typedef typename StreeOhleb<>::size_type size_type;
typedef uint8_t char_type;

// global data structures
string t, s;
StreeOhleb<> st;
Maxrep<StreeOhleb<>, sdsl::bit_vector> maxrep;
MsVectors<StreeOhleb<>, sdsl::bit_vector> ms_vec;


class InputFlags{
public:
    bool load_cst, load_maxrep;

    InputFlags() {}

    InputFlags(const bool load_cst_, const bool load_maxrep_) :
    	load_cst{load_cst_}, load_maxrep{load_maxrep_} {}

    InputFlags(const InputFlags& i) :
    	load_cst{i.load_cst}, load_maxrep{i.load_maxrep} {}

    InputFlags(const OptParser& args) :
    	load_cst {(args.getCmdOption("-load_cst") == "1")}, load_maxrep{(args.getCmdOption("-load_maxrep") == "1")} {}
};

class NodeProperty
{
private:
	Maxrep<StreeOhleb<>, sdsl::bit_vector> maxrep_;
    StreeOhleb<> st_;
public:
    NodeProperty(Maxrep<StreeOhleb<>, sdsl::bit_vector> maxrep, const StreeOhleb<>& st){
        maxrep_ = maxrep;
        st_ = st;
    }

    bool is_max(const node_type v) const { return maxrep_.is_maximal(v); }

    bool is_wide(const node_type v) const { return (((v.i)>>8) != ((v.j)>>8)); }

    bool has_wl(const node_type v, char c) const {
        node_type u = st_.double_rank_nofail_wl(v, c);
        return !st_.is_root(u);
    }

    string runs_node_label(const node_type v, const char_type c) const {
        string ch = {(char)c};
        string key = (ch + "_" +                               // char
                      "na" + "_" +                             // maximality
                      (has_wl(v, c) ? "wl" : "nowl") + "_" +   // has wl(c)
                      (is_wide(v) ? "wide" : "narrow"));       // interval width
        return key;
    }

    string ms_node_label(const node_type v, const char_type c) const {
        string ch = {(char)c};
        string key = (ch + "_" +                                   // char
                      (is_max(v) ? "maxrep" : "nomaxrep") + "_" +  // maximality
                      (has_wl(v, c) ? "wl" : "nowl") + "_" +       // has wl(c)
                      (is_wide(v) ? "wide" : "narrow"));           // interval width
        return key;
    }
};

class Stats
{
	typedef map<size_type, size_type> hist_type;
	typedef map<string, size_type> counter_type;

	private:
		size_type len_s, len_t;

		template<typename map_t>
		void dump_map(std::ostream& out, string measuring, string where, map_t& values){
			for(auto item: values)
				out << len_s << "," << len_t << "," << measuring << "," << where << "," << item.first << "," << item.second << endl;
		}

	public:
		hist_type runs_pcalls_seq, ms_pcalls_seq, ms_wlcalls_seq;
		counter_type runs_wl_calls, ms_wl_calls;

		Stats(size_type sl, size_type tl) : len_s{sl}, len_t{tl} {};

		void dump_stats(std::ostream &out){
			out << "len_s,len_t,measuring,where,key,value" << endl;
			dump_map<counter_type>(out, "wl_calls", "runs", runs_wl_calls);
			dump_map<counter_type>(out, "wl_calls", "ms", ms_wl_calls);

			dump_map<hist_type>(out, "pseq", "runs", runs_pcalls_seq);
			dump_map<hist_type>(out, "pseq", "ms", ms_pcalls_seq);
			dump_map<hist_type>(out, "wlseq", "ms", ms_wlcalls_seq);
		}
};


size_type build_runs(Stats& stats){
	NodeProperty NP(maxrep, st);
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
            size_type parent_seq_len = 0;
            do{ // remove suffixes of t[k..] until you can extend by 'c'
                v = st.parent(v);
				u = st.double_rank_nofail_wl(v, c);
				stats.runs_wl_calls[NP.runs_node_label(v, c)] += 1;  // record
                has_wl = !st.is_root(u);
                parent_seq_len += 1;
            } while(!has_wl && !st.is_root(v));
            stats.runs_pcalls_seq[parent_seq_len] += 1;  // record
        } else {
            ms_vec.runs[k] = 1;
        }
        v = st.double_rank_nofail_wl(v, c);
    }
}

void build_ms(Stats& stats){
    NodeProperty NP(maxrep, st);
    size_type k = 0, h_star = k + 1, h = h_star, h_star_prev = h_star, ms_idx = 0, ms_size = t.size();
    char_type c = t[k];
    node_type v = st.double_rank_nofail_wl(st.root(), c), u = v;
    bool is_maximal = true;

    while(k < t.size()){
        h = h_star;
        
        h_star_prev = h_star;
        while(h_star < ms_size){
            c = t[h_star];
			u = st.double_rank_nofail_wl(v, c);
			stats.ms_wl_calls[NP.ms_node_label(v,  c)] += 1; // record
            if(!st.is_root(u)){
                v = u;
                h_star += 1;
            } else
                break;
        }
        stats.ms_wlcalls_seq[(size_type) (h_star - h_star_prev)] += 1;  // record
        ms_vec.set_next_ms_values1(0, ms_idx, h, h_star, t.size() * 2);

        if(h_star < ms_size){ // remove prefixes of t[k..h*] until you can extend by 'c'
            if(!st.has_complete_info(v))
                st.lazy_wl_followup(v);

            bool has_wl = false;
            size_type parent_seq_len = 0;
            u = st.root();
            do{ // remove suffixes of t[k..] until you can extend by 'c'
                v = st.parent(v);
				u = st.double_rank_nofail_wl(v, c);
				stats.ms_wl_calls[NP.ms_node_label(v, c)] += 1;  // record

                has_wl = !st.is_root(u);
                parent_seq_len += 1;
            } while(!has_wl && !st.is_root(v));
            stats.ms_pcalls_seq[parent_seq_len] += 1;  // record
            h_star += 1;
        }
        k = ms_vec.set_next_ms_values2(0, ms_idx, k, t.size(), t.size() * 2);
        v = u;
    }
}

void comp(InputSpec& tspec, InputSpec& s_fwd, const string& out_path, InputFlags& flags){
    map<size_type, size_type> consecutive_runs_parent_calls, consecutive_ms_wl_calls, consecutive_ms_parent_calls;
    map<string, size_type> ms_wl_node_prop, runs_wl_node_prop, ms_rank_calls;

    /* prepare global data structures */
    cerr << "loading input ... " << endl;
    t = tspec.load_s();
    s = s_fwd.load_s();
    ms_vec = MsVectors<StreeOhleb<>, sdsl::bit_vector>(t.size(), 1);
    Stats stats(s.size(), t.size());

    /* build runs */
    cerr << "build runs ... " << endl;
    load_or_build(st, s, s_fwd.fwd_cst_fname, flags.load_cst);
    build_runs(stats);

    /* reverse s */
    cerr << "reversing index of length " << s.size() << " ... " << endl;
    InputSpec::reverse_in_place(s);

    /* build ms */
    cerr << "build ms ... " << endl;
    load_or_build(st, s, s_fwd.rev_cst_fname, flags.load_cst);
    Maxrep<StreeOhleb<>, sdsl::bit_vector>::load_or_build(maxrep, st, s_fwd.rev_maxrep_fname, flags.load_maxrep);
    build_ms(stats);


    cerr << "dumping results ... " << endl;
	stats.dump_stats(cout);
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
        flags = InputFlags(false, false);
    } else {
        tspec = InputSpec(input.getCmdOption("-t_path"));
        sfwd_spec = InputSpec(input.getCmdOption("-s_path"));
        flags = InputFlags(input);
    }
    comp(tspec, sfwd_spec, out_path, flags);
    return 0;
}
