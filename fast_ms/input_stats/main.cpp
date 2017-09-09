//
//  main.cpp
//  input_stats
//
//  Created by denas on 4/5/17.
//  Copyright © 2017 denas. All rights reserved.
//

#include <iostream>
#include <string>

#include <sdsl/suffix_trees.hpp>

#include "utils.hpp"
#include "stree_sct3.hpp"
#include "fd_ms.hpp"


#define IS_INTNODE_MAXIMAL(node)                    ((maxrep[(node).i] == 1) && (maxrep[(node).j] == 1) )
#define IS_MAXIMAL(node) ( ((node).i != (node).j) && (maxrep[(node).i] == 1) && (maxrep[(node).j] == 1) )

//#define IS_WIDE(i,j)  (((i)>>8) != ((j)>>8))
//#define IS_MAXREP(v)  ((maxrep[v.i] == 1) && (maxrep[v.j] == 1))

using namespace std;
using namespace fdms;

typedef typename StreeOhleb<>::node_type node_type;
enum class IntervalWidth {same_block, different_blocks};


class NodeProperty
{
private:
    sdsl::bit_vector maxrep_;
    StreeOhleb<> st_;
public:
    NodeProperty(const sdsl::bit_vector& maxrep, const StreeOhleb<>& st){
        maxrep_ = maxrep;
        st_ = st;
    }

    bool is_max(const node_type v) const {
        return ((v.i < v.j) && (maxrep_[v.i] == 1) && (maxrep_[v.j] == 1));
    }

    bool is_wide(const node_type v) const {
        return (((v.i)>>8) != ((v.j)>>8));
        //return v.i < v.j;
    }

    bool has_wl(const node_type v, char c) const {
        node_type u = st_.double_rank_nofail_wl(v, c);
        return !st_.is_root(u);
    }

    string node_label(const node_type v, const char_type c) const {
        string ch = {(char)c};
        string key = (ch + "_" +                              // char
                      (is_max(v) ? "maxrep" : "nomaxrep") + "_" +  // maximality
                      (has_wl(v, c) ? "wl" : "nowl") + "_" + // has wl(c)
                      (is_wide(v) ? "wide" : "narrow"));      // interval width
        return key;
    }
};


Interval bstep(const StreeOhleb<> &st_, Interval &I, char_type c){
    int cc = st_.csa.char2comp[c];
    return std::make_pair(st_.csa.C[cc] + st_.csa.bwt.rank(I.first, c),
                          st_.csa.C[cc] + st_.csa.bwt.rank(I.second + 1, c) - 1);
}

Interval bvector_composition(sdsl::bit_vector v){
    size_type ones = 0;
    for(size_type i = 0; i < v.size(); i++)
        if(v[i] == 1)
            ones += 1;
    Interval result =  std::make_pair(ones, v.size() - ones);
    assert(result.first + result.second == v.size());
    return result;
}

template <typename map_t>
void dump_map(size_type s_size, size_type t_size, string measuring, string where, map_t values){
    for(auto item: values)
        cout << s_size << "," << t_size << "," << measuring << "," << where << "," << item.first << "," << item.second << endl;
}



void build_runs(const string& t, const StreeOhleb<>& st, sdsl::bit_vector& runs,
                map<size_type, size_type>& consecutive_parent_calls){
    size_type k = t.size();
    char_type c = t[k - 1];
    Interval I = make_pair(st.csa.C[st.csa.char2comp[c]], st.csa.C[st.csa.char2comp[c] + 1] - 1);
    node_type v = st.double_rank_nofail_wl(st.root(), c);

    while(--k > 0){
        c = t[k-1];
        I = bstep(st, I, c);
        if(I.first > I.second){
            runs[k] = 0;
            size_type consecutive_parent_calls_cnt = 0;
            do{ // update I to the parent of the proper locus of w until we can extend by 'c'
                v = st.parent(v);

                consecutive_parent_calls_cnt += 1;
                I = make_pair(v.i, v.j);
                I = bstep(st, I, c);
            } while(I.first > I.second);
            consecutive_parent_calls[consecutive_parent_calls_cnt] += 1;
        } else
            runs[k] = 1;
        v = st.double_rank_nofail_wl(v, c); // update v
    }
}

void build_ms(const string& t, const StreeOhleb<>& st, sdsl::bit_vector& runs, sdsl::bit_vector& ms, sdsl::bit_vector& maxrep,
              map<size_type, size_type>& consecutive_wl_calls, map<size_type, size_type>& consecutive_parent_calls,
              map<string, size_type>& maximal_visits, map<string, size_type>& rank_calls){

    NodeProperty NP(maxrep, st);

    size_type k = 0, h_star = k + 1, h = h_star, h_star_prev = h_star, k_prim, ms_idx = 0, ms_size = t.size() ;
    uint8_t c = t[k];
    rank_calls[NP.node_label(st.root(), c)] += 1;
    node_type v = st.double_rank_fail_wl(st.root(), c), u = v;
    
    bool is_maximal = true;

    while(k < ms_size){
        h = h_star;

        h_star_prev = h_star;
        while(h_star < ms_size){
            c = t[h_star];
            is_maximal = IS_MAXIMAL(v);
            u = st.double_rank_fail_wl_mrep(v, c, is_maximal);
            rank_calls[NP.node_label(v, c)] += 1;
            maximal_visits[NP.node_label(v, c)] += 1;
            if(!st.is_root(u)){
                v = u;
                h_star += 1;
            } else
                break;
        }

        assert(h_star >= h_star_prev);
        consecutive_wl_calls[(size_type) (h_star - h_star_prev)] += 1;

        ms_idx += (h_star -  h + 1);
        if(h_star - h + 1 > 0)
            ms[ms_idx++] = 1;

        h_star_prev = h_star;
        if(h_star < ms_size){
            size_type consecutive_parent_calls_cnt = 0;
            is_maximal = false;
            bool has_wl = false;
            u = st.root();
            do{ // remove suffixes of t[k..] until you can extend by 'c'
                v = st.parent(v);
                consecutive_parent_calls_cnt += 1;
                
                if(!is_maximal)
                    is_maximal = IS_INTNODE_MAXIMAL(v); //since parent of a maximal is a maximal
                if(is_maximal){
                    rank_calls[NP.node_label(v, c)] += 1;
                    u = st.double_rank_fail_wl_mrep(v, c, is_maximal);
                    has_wl = !st.is_root(u);
                }
            } while(!has_wl && !st.is_root(v)); // since !maximal => no wl
            h_star += 1;

            h_star +=  1;
            consecutive_parent_calls[consecutive_parent_calls_cnt] += 1;
        }
        assert(h_star >= h_star_prev);

        // k_prim: index of the first zero to the right of k in runs
        k_prim = find_k_prim_(k, ms_size, runs);

        for(size_type i = k + 1; i <= k_prim - 1; i++)
            ms[ms_idx++] = 1;
        k = k_prim;
        v = u;
    }

}

void comp(InputSpec& T, InputSpec& S_fwd, const string& out_path, InputFlags& flags){
    string t, s;
    map<size_type, size_type> consecutive_runs_parent_calls, consecutive_ms_wl_calls, consecutive_ms_parent_calls;
    //map<IntervalWidth, size_type> ms_interval_width, runs_interval_width;
    map<string, size_type> ms_wl_node_prop, runs_wl_node_prop, ms_rank_calls;


    /* load input */
    cerr << "loading input ... ";
    t = T.load_s();
    s = S_fwd.load_s();
    cerr << "DONE" << endl;

    /* build stree */
    StreeOhleb<> st;
    load_st(st, s, S_fwd.fwd_cst_fname, flags.load_stree);
    cerr << "DONE" << endl;

    /* prepare global data structures */
    sdsl::bit_vector runs(t.size());
    sdsl::util::set_to_value(runs, 0);

    sdsl::bit_vector ms(2 * t.size());
    sdsl::util::set_to_value(ms, 0);

    sdsl::bit_vector maxrep(st.size());
    sdsl::util::set_to_value(maxrep, 0);


    cerr << "build runs ... ";
    //runs_interval_width[IntervalWidth::different_blocks] = runs_interval_width[IntervalWidth::same_block] = 0;
    //runs_maximal_visits[Maximality::maximal] = runs_maximal_visits[Maximality::non_maximal] = 0;
    build_runs(t, st, runs, consecutive_runs_parent_calls);
    cerr << "DONE" << endl;

    /* reverse s */
    cerr << "reversing string s of length " << s.size() << " ... ";
    reverse_in_place(s);
    cerr << "DONE" << endl;

    /* build the cst */
    load_st(st, s, S_fwd.rev_cst_fname, flags.load_stree);
    cerr << "DONE" << endl;

    /* build the maxrep vector */
    size_type time_usage_ms_maxrep = load_maxrep(maxrep, st, s, S_fwd.rev_maxrep_fname, flags.load_maxrep);
    cerr << "DONE (" << time_usage_ms_maxrep / 1000 << " seconds)" << endl;


    cerr << "build ms ... ";
    //ms_interval_width[IntervalWidth::different_blocks] = ms_interval_width[IntervalWidth::same_block] = 0;
    build_ms(t, st, runs, ms, maxrep, consecutive_ms_wl_calls, consecutive_ms_parent_calls, ms_wl_node_prop, ms_rank_calls);
    cerr << "DONE" << endl;


    cout << "len_s,len_t,measuring,where,key,value" << endl;

    //for(auto item: runs_interval_width)
    //    cout << s.size() << "," << t.size() << "," << "interval_width" << "," << "runs" << "," << (item.first == IntervalWidth::different_blocks ? "large" : "small") << "," << item.second << endl;
    //for(auto item: ms_interval_width)
    //   cout << s.size() << "," << t.size() << "," << "interval_width" << "," << "ms" << "," << (item.first == IntervalWidth::different_blocks ? "large" : "small") << "," << item.second << endl;
    for(auto item: ms_wl_node_prop)
        cout << s.size() << "," << t.size() << "," << "wlnode_prop" << "," << "ms" << "," << item.first << "," << item.second << endl;

    for(auto item: ms_rank_calls)
        cout << s.size() << "," << t.size() << "," << "rank_call" << "," << "ms" << "," << item.first << "," << item.second << endl;

    Interval runs_comp = bvector_composition(runs);
    cout << s.size() << "," << t.size() << ",vector_composition,runs,0," << runs_comp.second << endl;
    cout << s.size() << "," << t.size() << ",vector_composition,runs,1," << runs_comp.first << endl;

    Interval maxrep_comp = bvector_composition(maxrep);
    cout << s.size() << "," << t.size() << ",vector_composition,maxrep,maximal," << maxrep_comp.first << endl;
    cout << s.size() << "," << t.size() << ",vector_composition,maxrep,non_maximal," << maxrep_comp.second << endl;

    dump_map(s.size(), t.size(), "consecutive_parent_calls", "runs", consecutive_runs_parent_calls);
    dump_map(s.size(), t.size(), "consecutive_parent_calls", "ms", consecutive_ms_parent_calls);
    dump_map(s.size(), t.size(), "consecutive_wl_calls", "ms", consecutive_ms_wl_calls);

}

int main(int argc, char **argv){
    OptParser input(argc, argv);
    if(argc == 1){
        const string base_dir = {"/Users/denas/projects/matching_statistics/indexed_ms/tests/"};
        InputFlags flags(false, // lazy_wl
                         false, // sada cst
                         false, // maxrep
                         true,  // lca_parents
                         true,  // space
                         false, // time
                         true,  // ans
                         false, // verbose
                         10,    // nr. progress messages for runs construction
                         10,    // nr. progress messages for ms construction
                         false, // load CST
                         false, // load MAXREP
                         1      // nthreads
                         );
        InputSpec tspec(base_dir + "t");
        InputSpec sfwd_spec(base_dir + "s");
        const string out_path = "0";
        comp(tspec, sfwd_spec, out_path, flags);
    } else {
        InputFlags flags(input);
        InputSpec tspec(input.getCmdOption("-t_path"));
        InputSpec sfwd_spec(input.getCmdOption("-s_path"));
        const string out_path = input.getCmdOption("-out_path");
        comp(tspec, sfwd_spec, out_path, flags);
    }
    
    return 0;
}
