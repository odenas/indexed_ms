/*
 * fabio_djamal_ms.cpp
 *
 *  Created on: Oct 13, 2016
 *      Author: denas
 */

#include <iostream>
#include <string>
#include <future>
#include <thread>

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

typedef pair<size_type, size_type> IInterval;
typedef typename StreeOhleb<>::node_type node_type;

string t, s;
StreeOhleb<> st;
bvector runs(1);
//bvector ms(1);
vector<bvector> mses(1); // the ms vector for each thread

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

IInterval init_interval(StreeOhleb<>& st_, const char c){
    int cc = st_.csa.char2comp[c];
    return std::make_pair(st_.csa.C[cc], st_.csa.C[cc + 1] - 1);
}

IInterval bstep_interval(StreeOhleb<>& st_, IInterval& cur_i, char c){
    int cc = st_.csa.char2comp[c];
    return std::make_pair(st_.csa.C[cc] + st_.csa.bwt.rank(cur_i.first, c),
                          st_.csa.C[cc] + st_.csa.bwt.rank(cur_i.second + 1, c) - 1);
}

std::vector<pair<size_type, size_type>> slice_input(const size_type input_size, const size_type nthreads){
    size_type chunk = input_size / nthreads;
    size_type extra = input_size % nthreads;

    std::vector<pair<size_type, size_type>> slices (nthreads);
    for(size_type i=0, from = 0; i<nthreads; i++){
        slices[i] = std::make_pair(from, from + chunk + (i < extra ? 1 : 0));
        from += (from + chunk + (i < extra ? 1 : 0));
    }
    return slices;
}

node_type parent_sequence(StreeOhleb<>& st_, node_type v, IInterval& I, const size_type c){
    do{
        v = st_.parent(v);
        I = std::make_pair(v.i, v.j);
        I = bstep_interval(st_, I, c);
    } while(I.first > I.second);
    return v;
}

int fill_runs_slice(const size_type from, const size_type to){
    // [from, to)
    size_type k = to, c = t[k - 1];
    IInterval I = init_interval(st, static_cast<char>(c));
    node_type v = st.wl(st.root(), c); // stree node

    while(--k > 0){
        c = t[k-1];
        I = bstep_interval(st, I, c);
        if(I.first > I.second){ // empty
            runs[k] = 0;
            // remove suffixes of t[k..] until you can extend by 'c'
            v = parent_sequence(st, v, I, c);
        } else {
            runs[k] = 1;
        }
        v = st.wl(v, c); // update v
    }
    return 0;
}

// TODO: refactor string and stree initialization

performance_monitor build_runs_ohleb(const InputFlags& flags, const InputSpec &s_fwd){
    monitor::size_dict space_usage, time_usage;

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

    /* compute RUNS slices */
    runs_start = timer::now();
    cerr << "building RUNS over " << 1 << " threads ...";
    std::future<int> result = std::async(std::launch::async, fill_runs_slice, 0, t.size());
    result.get();

    runs_stop = timer::now();
    time_usage["alg"]  = std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count();
    space_usage["alg"] = 0;
    cerr << "DONE (" << time_usage["alg"] / 1000 << " seconds)" << endl;
    return std::make_pair(space_usage, time_usage);
}

size_type fill_ms_slice(const size_type mses_idx, const size_type from, const size_type to){
    size_type k = from, h_star = k + 1, h = h_star, h_star_prev = h_star, k_prim, ms_idx = 0, ms_size = t.size();
    uint8_t c = t[k];
    IInterval I = init_interval(st, static_cast<char>(c));
    node_type v = st.wl(st.root(), c);

    while(k < to){
        h = h_star;
        h_star_prev = h_star;

        for(; I.first <= I.second && h_star < ms_size; ){
            c = t[h_star];
            I = bstep_interval(st, I, c);
            if(I.first <= I.second){
                v = st.wl(v, c);
                h_star++;
            }
        }
        ms_idx += (h_star -  h + 1);
        if(h_star - h + 1 > 0)
            mses[mses_idx][ms_idx++] = 1;

        if(h_star < ms_size){ // remove prefixes of t[k..h*] until you can extend by 'c'
            v = parent_sequence(st, v, I, t[h_star]);
            h_star += 1;
        }
        // k_prim: index of the first zero to the right of k in runs
        k_prim = find_k_prim_(k, ms_size, runs);

        for(size_type i = k + 1; i <= k_prim - 1 && i < to; i++)
            mses[mses_idx][ms_idx++] = 1;
        v = st.wl(v, c);
        k = k_prim;
    }
    mses[mses_idx].resize(ms_idx);
    return 0;
}

performance_monitor build_ms_ohleb(const InputFlags& flags, InputSpec &s_fwd){
    monitor::size_dict space_usage, time_usage;

    auto runs_start = timer::now();
    if(flags.load_stree){
        cerr << "loading the CST T(s') from " << s_fwd.s_fname + ".rev.stree" << "... ";
        sdsl::load_from_file(st, s_fwd.s_fname + ".rev.stree");
    } else {
        cerr << "building the CST T(s') of lentgth " << s.size() << "... ";
        sdsl::construct_im(st, s, 1);
    }
    auto runs_stop = timer::now();
    time_usage["dstruct"] = std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count();
    cerr << "DONE (" << time_usage["dstruct"] / 1000 << " seconds, " << st.size() << " nodes)" << endl;

    /* build MS */
    cerr << "building MS in " << (flags.lazy ? "" : "non-") << "lazy mode over " << flags.nthreads << " threads ..." << endl;
    runs_start = timer::now();
    std::vector<IInterval> slices = slice_input(t.size(), flags.nthreads);
    std::vector<std::future<size_type>> results(flags.nthreads);
    for(size_type i=0; i<flags.nthreads; i++){
        cerr << " * launching ms computation over : [" << slices[i].first << " .. " << slices[i].second << ")" << endl;
        //fill_ms_slice(i, slices[i].first, slices[i].second);
        results[i] = std::async(std::launch::async, fill_ms_slice, i, slices[i].first, slices[i].second);
    }
    for(size_type i=0; i<flags.nthreads; i++)
        results[i].get();
    size_type total_ms_length = 0;
    for(size_type i=0; i<flags.nthreads; i++)
        total_ms_length += mses[i].size();
    cerr << " * total ms length : " << total_ms_length << " (with |t| = " << t.size() << ")" << endl;
    runs_stop = timer::now();
    time_usage["alg"]            = std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count();
    space_usage["stree_csa"]     = sdsl::size_in_bytes(st.csa);
    space_usage["stree_bp"]      = sdsl::size_in_bytes(st.bp);
    space_usage["stree_bpsupp"]  = sdsl::size_in_bytes(st.bp_support);
    cerr << "DONE (" << time_usage["alg"] / 1000 << " seconds)" << endl;

    return std::make_pair(space_usage, time_usage);
}


void comp(InputSpec& T, InputSpec& S_fwd, const string& out_path, InputFlags& flags){
    monitor::size_dict space_usage, time_usage;
    performance_monitor runs_usage, ms_usage;

    t = T.load_s();
    s = S_fwd.load_s();
    runs.resize(t.size()); sdsl::util::set_to_value(runs, 0); //bvector runs(t.size());
    mses.resize(flags.nthreads);
    for(int i=0; i<flags.nthreads; i++){
        mses[i].resize(2 * t.size());
        sdsl::util::set_to_value(mses[i], 0);
    }

    if(flags.ms_progress > t.size())
        flags.ms_progress = t.size() - 1;

    runs_usage = build_runs_ohleb(flags, S_fwd);
    reverse_in_place(s);
    ms_usage = build_ms_ohleb(flags, S_fwd);

    // merge space usage to compute peak usage
    space_usage["s"]    = s.size();
    space_usage["t"]    = t.size();
    space_usage["runs"] = runs.size();
    space_usage["ms"]   = (2 * t.size()) * flags.nthreads;
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
        if(out_path == "0"){
            for(size_type mses_idx=0; mses_idx < mses.size(); mses_idx++){
                dump_ms(mses[mses_idx]);
            }
            cout << endl;
        }
        //dump_ms(ms);
        //else{
        //    sdsl::int_vector<32> MS(t.size());
        //    size_type k = 0;
        //    size_type j = 0;
        //    for(size_type i = 0; i < ms.size(); i++){
        //        if(ms[i] == 1){
        //            MS[j++] = (uint32_t)(i - (2*k));
        //            k += 1;
        //        }
        //    }
        //    cerr << "dumping binary MS array : " << "(j, k): (" << ", " << j << ", " << k << ")" << endl;
        //    sdsl::store_to_file(MS, out_path);
        //}
    }

}


int main(int argc, char **argv){
    InputParser input(argc, argv);
    if(argc == 1){
        const string base_dir = {"/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/"};
        InputFlags flags(false, // lazy_wl
                         false, // sada cst
                         false, // space
                         false,  // time
                         true,  // ans
                         false, // verbose
                         10,    // nr. progress messages for runs construction
                         10,    // nr. progress messages for ms construction
                         false, // load CST
                         2      // nthreads
                         );
        InputSpec tspec(base_dir + "abcde200_256t.txt");
        InputSpec sfwd_spec(base_dir + "abcde200_256s.txt");
        const string out_path = "0";
        comp(tspec, sfwd_spec, out_path, flags);
    } else {
        InputFlags flags(input);
        InputSpec tspec(input.getCmdOption("-t_path"));
        InputSpec sfwd_spec(input.getCmdOption("-s_path"));
        const string out_path = input.getCmdOption("-out_path");
        //cerr << "**** nthreads = " << flags.nthreads << " ****" << endl;
        comp(tspec, sfwd_spec, out_path, flags);
    }
    return 0;
}

