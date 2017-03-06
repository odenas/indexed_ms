/*
 * fabio_djamal_ms.cpp
 *
 *  Created on: Oct 13, 2016
 *      Author: denas
 */

#include <iostream>
#include <string>

#include <sdsl/suffix_trees.hpp>

#include "basic.hpp"
#include "fd_ms.hpp"
#include "stree_sada.hpp"
#include "stree_sct3.hpp"


//#define STREE_SADA
#define NR_REPORTS 10

using namespace std;
using namespace fdms;
using timer = std::chrono::high_resolution_clock;


bvector construct_bp(string& _s){
    sdsl::cst_sada<> temp_st;
    sdsl::construct_im(temp_st, _s, 1);
    bvector _bp(temp_st.bp.size());
    for(int i = 0; i < _bp.size(); i++)
        _bp[i] = temp_st.bp[i];
    return _bp;
}

void reverse_in_place(string& s){
    size_type n = s.size();

    for(int i = 0; i < n / 2; i++){
        char c = s[i];
        s[i] = s[n - 1 - i];
        s[n - 1 - i] = c;
    }
}


/* find k': index of the first zero to the right of k in runs */
size_type find_k_prim_(size_type __k, size_type max__k, bvector& __runs){
    while(++__k < max__k && __runs[__k] != 0)
        ;
    return __k;
}

void report_progress(timer::time_point start_time, size_type curr_idx, size_type total){
    timer::time_point stop_time = timer::now();
    size_type elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time).count() + 1;
    cerr << endl << "[" << elapsed / 1000 << " s] " << 100.0 * curr_idx / total << "% @ " << (1.0 * curr_idx / elapsed) << " KHz";
}

performance_monitor build_runs_ohleb(string& t, string& s, bvector& runs, const InputFlags& flags, const InputSpec &s_fwd){
    typedef pair<size_type, size_type> IInterval;
    monitor::size_dict space_usage, time_usage;
    StreeOhleb<> st;

    auto init_interval = [] (StreeOhleb<>& st_, char c) -> IInterval {
        int cc = st_.csa.char2comp[c];
        return std::make_pair(st_.csa.C[cc], st_.csa.C[cc + 1] - 1);
    };

    auto bstep_interval = [] (StreeOhleb<>& st_, IInterval& cur_i, char c) -> IInterval{
        int cc = st_.csa.char2comp[c];
        return std::make_pair(st_.csa.C[cc] + st_.csa.bwt.rank(cur_i.first, c),
                              st_.csa.C[cc] + st_.csa.bwt.rank(cur_i.second + 1, c) - 1);
    };


    /* build the CST */
    auto runs_start = timer::now();
    if(flags.load_stree){
        cerr << "loading the CST T(s) from " << s_fwd.s_fname + ".fwd.stree" << "... ";
        sdsl::load_from_file(st, s_fwd.s_fname + ".fwd.stree");
    } else {
        cerr << "building the CST T(s) of lentgth " << s.size() << "... ";
        sdsl::construct_im(st, s, 1);
    }
    auto runs_stop = timer::now();
    time_usage["dstruct"]       = std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count();
    cerr << "DONE (" << time_usage["dstruct"] / 1000 << " seconds, " << st.size() << " nodes)" << endl;

    space_usage["stree_csa"]    = sdsl::size_in_bytes(st.csa);
    space_usage["stree_bp"]     = sdsl::size_in_bytes(st.bp);
    space_usage["stree_bpsupp"] = sdsl::size_in_bytes(st.bp_support);

    /* compute RUNS */
    runs_start = timer::now();
    cerr << "building RUNS ... ";
    size_type k = t.size(), c = t[k - 1];
    IInterval I = init_interval(st, static_cast<char>(c)); //Interval I{st.csa, static_cast<char>(c)};
    typedef typename StreeOhleb<>::node_type node_type;

    node_type v = st.wl(st.root(), c); // stree node
    while(--k > 0){
        c = t[k-1];
        I = bstep_interval(st, I, c); //I.bstep(c);
        if(I.first > I.second){ // empty
            runs[k] = 0;
            // update I to the parent of the proper locus of w until we can extend by 'c'
            do{
                v = st.parent(v);
                I.first = st.lb(v); I.second = st.rb(v); //I.set(st.lb(v), st.rb(v));
                I = bstep_interval(st, I, c); //I.bstep(c);
            } while(I.first > I.second);
        } else {
            runs[k] = 1;
        }
        v = st.wl(v, c); // update v

        if (flags.runs_progress > 0 && k % (t.size() / flags.runs_progress) == 0)
            report_progress(runs_start, t.size() - k, t.size());
    }

    runs_stop = timer::now();
    time_usage["alg"]  = std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count();
    space_usage["alg"] = 0;
    cerr << "DONE (" << time_usage["alg"] / 1000 << " seconds)" << endl;
    return std::make_pair(space_usage, time_usage);
}


performance_monitor build_ms_ohleb(string& t, string& s_rev, bvector& runs, bvector& ms, const InputFlags& flags, InputSpec &s_fwd){
    monitor::size_dict space_usage, time_usage;
    typedef pair<size_type, size_type> IInterval;
    typedef typename StreeOhleb<>::node_type node_type;

    auto init_interval = [] (StreeOhleb<>& st_, char c) -> IInterval {
        int cc = st_.csa.char2comp[c];
        return std::make_pair(st_.csa.C[cc], st_.csa.C[cc + 1] - 1);
    };

    auto bstep_interval = [] (StreeOhleb<>& st_, IInterval& cur_i, char c) -> IInterval{
        int cc = st_.csa.char2comp[c];
        return std::make_pair(st_.csa.C[cc] + st_.csa.bwt.rank(cur_i.first, c),
                              st_.csa.C[cc] + st_.csa.bwt.rank(cur_i.second + 1, c) - 1);
    };

    StreeOhleb<> st;

    auto runs_start = timer::now();
    if(flags.load_stree){
        cerr << "loading the CST T(s') from " << s_fwd.s_fname + ".rev.stree" << "... ";
        sdsl::load_from_file(st, s_fwd.s_fname + ".rev.stree");
    } else {
        cerr << "building the CST T(s') of lentgth " << s_rev.size() << "... ";
        sdsl::construct_im(st, s_rev, 1);
    }
    auto runs_stop = timer::now();
    time_usage["dstruct"] = std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count();
    cerr << "DONE (" << time_usage["dstruct"] / 1000 << " seconds, " << st.size() << " nodes)" << endl;

    /* build MS */
    cerr << "building MS in " << (flags.lazy ? "" : "non-") << "lazy mode ... ";
    runs_start = timer::now();
    size_type k = 0, h_star = k + 1, h = h_star, h_star_prev = h_star, k_prim, ms_idx = 0, ms_size = t.size() ;
    uint8_t c = t[k];
    IInterval I = init_interval(st, static_cast<char>(c)); //Interval I{bwt, static_cast<char>(c)};

    std::map<int, int> consecutive_lazy_wl_calls;
    node_type v = st.wl(st.root(), c); // stree node

    while(k < ms_size){
        h = h_star;
        h_star_prev = h_star;


        if(flags.lazy){ // lazy weiner links
            for(; I.first <= I.second && h_star < ms_size; ){
                c = t[h_star];
                I = bstep_interval(st, I, c);
                if(I.first <= I.second){
                    v = st.lazy_wl(v, c);
                    h_star++;
                }
            }
            if(h_star > h_star_prev) // we must have called lazy_wl(). complete the node
                st.lazy_wl_followup(v);
        } else { // non-lazy weiner links
            for(; I.first <= I.second && h_star < ms_size; ){
                c = t[h_star];
                I = bstep_interval(st, I, c);
                if(I.first <= I.second){
                    v = st.wl(v, c);
                    h_star++;
                }
            }
        }

        consecutive_lazy_wl_calls[(int) (h_star - h_star_prev)] += 1;
        ms_idx += (h_star -  h + 1);
        if(h_star - h + 1 > 0)
            ms[ms_idx++] = 1;

        if(h_star < ms_size){
            do {
                v = st.parent(v);
                I.first = st.lb(v); I.second = st.rb(v); //I.set(st.lb(v), st.rb(v));
                I = bstep_interval(st, I, t[h_star]); //I.bstep(t[h_star]);
            } while(I.first > I.second);
            h_star = h_star + 1;
        }
        // k_prim: index of the first zero to the right of k in runs
        k_prim = find_k_prim_(k, ms_size, runs);

        for(size_type i = k + 1; i <= k_prim - 1; i++)
            ms[ms_idx++] = 1;

        if (flags.ms_progress > 0 &&  k % (ms_size / flags.ms_progress) > k_prim % (ms_size / flags.ms_progress)){
            report_progress(runs_start, k_prim, ms_size);
            cerr << " (k  --> k', h*): " << k << ", " << k_prim << "," << h_star;
        }

        // update v
        v = st.wl(v, c);
        k = k_prim;
    }

    runs_stop = timer::now();
    time_usage["alg"]            = std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count();
    for(auto item: consecutive_lazy_wl_calls)
        time_usage["consecutive_lazy_wl_calls" + std::to_string(item.first)] = item.second;
    space_usage["stree_csa"]     = sdsl::size_in_bytes(st.csa);
    space_usage["stree_bp"]      = sdsl::size_in_bytes(st.bp);
    space_usage["stree_bpsupp"]  = sdsl::size_in_bytes(st.bp_support);
    cerr << "DONE (" << time_usage["alg"] / 1000 << " seconds)" << endl;

    return std::make_pair(space_usage, time_usage);
}


void comp(InputSpec& T, InputSpec& S_fwd, const string& out_path, InputFlags& flags){
    monitor::size_dict space_usage, time_usage;
    performance_monitor runs_usage, ms_usage;

    string t = T.load_s();
    string s = S_fwd.load_s();
    bvector runs(t.size());
    bvector ms(t.size() * 2);

    if(flags.ms_progress > t.size())
        flags.ms_progress = t.size() - 1;

    runs_usage = build_runs_ohleb(t, s, runs, flags, S_fwd);
    reverse_in_place(s);
    ms_usage = build_ms_ohleb(t, s, runs, ms, flags, S_fwd);

    // merge space usage to compute peak usage
    space_usage["s"]    = s.size();
    space_usage["t"]    = t.size();
    space_usage["runs"] = runs.size();
    space_usage["ms"]   = ms.size();
    {
        unsigned long rss, vs;
        getmem(&rss, &vs);
        space_usage["resident_mem"] = rss;
        space_usage["virtual_mem"]  = vs;
    }
    // space usage is max(space in runs, space in ms)
    for(auto item : ms_usage.first)
        space_usage[item.first] = max(ms_usage.first[item.first], runs_usage.first[item.first]);


    if(flags.space_or_time_usage){
        cout << "len_s,len_t,measuring,item,value" << endl;
        if(flags.space_usage){
            for(auto item: space_usage)
                cout << s.size() << "," << t.size() << ",space," << item.first << "," << item.second << endl;
        }
        if(flags.time_usage){
            for(auto item : runs_usage.second)
                cout << s.size() << "," << t.size() << ",time,runs_" << item.first << "," << item.second << endl;

            for(auto item : ms_usage.second)
                cout << s.size() << "," << t.size() << ",time,ms_" << item.first << "," << item.second << endl;
        }
    }

    if(flags.answer){
        if(out_path == "0")
            dump_ms(ms);
        else{
            sdsl::int_vector<32> MS(t.size());
            size_type k = 0;
            size_type j = 0;
            for(size_type i = 0; i < ms.size(); i++){
                if(ms[i] == 1){
                    MS[j++] = (uint32_t)(i - (2*k));
                    k += 1;
                }
            }
            cerr << "dumping binary MS array : " << "(j, k): (" << ", " << j << ", " << k << ")" << endl;
            sdsl::store_to_file(MS, out_path);
        }
    }

}


int main(int argc, char **argv){
    InputParser input(argc, argv);
    if(argc == 1){
        const string base_dir = {"/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/lazy_vs_nonlazy_data/input_data/"};
        InputFlags flags(false, // lazy_wl
                         false, // sada cst
                         false, // space
                         true, // time
                         false, // ans
                         false, // verbose
                         0,     // nr. progress messages for runs construction
                         0,     // nr. progress messages for ms construction
                         true  // load CST
                         );
        //InputSpec tspec(base_dir + "../Homo_sapiens.GRCh38.dna.chromosome.22.juststring");
        //InputSpec sfwd_spec(base_dir + "Mus_musculus.GRCm38.dna.chromosome.MT.juststring");
        InputSpec tspec(base_dir + "mut_100Ms_5Mt_10.t");
        InputSpec sfwd_spec(base_dir + "mut_100Ms_5Mt_10.s");
        const string out_path = "";
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

