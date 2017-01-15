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

/* compute MS[k] from a select-index of the ms-vector*/
size_type get_ms(sdsl::select_support_mcl<1,1>& __ms_select1, size_type __k) {
    if(__k == -1)
        return (size_type) 1;
    return __ms_select1(__k + 1) - (2 * __k);
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

performance_monitor build_ms_sada(const string& prefix, string& t, string& s_rev, bvector& runs, bvector& ms, const InputFlags& flags){
    typedef fdms::bp_support_g<> t_bp_support;
    monitor::size_dict space_usage, time_usage;

    auto runs_start = timer::now();
    Bwt bwt(s_rev);
    bvector bp = construct_bp(s_rev);
    t_bp_support bp_supp(&bp);
    StreeSada<t_bp_support> st(bp_supp, bwt);
    auto runs_stop = timer::now();
    time_usage["dstruct"] = std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count();

    runs_start = timer::now();
    size_type size_in_bytes_alg = build_ms_from_st_and_bwt(st, bwt, t, prefix, runs, ms, flags.lazy);
    runs_stop = timer::now();
    time_usage["alg"] = std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count();

    space_usage["bwt_wtree"]           = bwt.size_in_bytes__wtree;
    space_usage["bwt_bwt"]             = bwt.size_in_bytes__bwt;
    space_usage["bwt_alp"]             = bwt.size_in_bytes__C + bwt.size_in_bytes__char2int + bwt.size_in_bytes__Sigma;
    space_usage["stree_bp"]            = sdsl::size_in_bytes(bp);
    space_usage["stree_bpsupp"]        = bp_supp.size_in_bytes_total;
    space_usage["stree_bpsupp_select"] = st.size_in_bytes__select;
    space_usage["stree_bpsupp_rank"]   = st.size_in_bytes__rank;
    space_usage["alg"]                 = size_in_bytes_alg;

    return std::make_pair(space_usage, time_usage);
}

performance_monitor build_runs_sada(const string& prefix, string& t, string& s, bvector& runs, const InputFlags& flags){
    typedef fdms::bp_support_g<> t_bp_support;
    typedef pair<size_type, size_type> IInterval;
    monitor::size_dict space_usage, time_usage;

    auto init_interval = [] (Bwt& bwt_, char c) -> IInterval {
        int cc = bwt_.char2int[c];
        return std::make_pair(bwt_.C[cc], bwt_.C[cc + 1]  - 1);
    };

    auto bstep_interval = [] (Bwt& bwt_, IInterval& cur_i, char c) -> IInterval{
        int cc = bwt_.char2int[c];
        return std::make_pair(bwt_.C[cc] + bwt_.rank(cur_i.first, c),
                              bwt_.C[cc] + bwt_.rank(cur_i.second + 1, c) - 1);
    };


    /* build the CST */
    cerr << "building the BWT and CST ...";
    auto runs_start = timer::now();

    Bwt bwt(s);
    bvector bp = construct_bp(s);
    t_bp_support bp_supp(&bp);
    StreeSada<t_bp_support> st(bp_supp, bwt);
    auto runs_stop = timer::now();
    time_usage["dstruct"]        = std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count();
    space_usage["bwt"]           = (bwt.size_in_bytes__wtree + bwt.size_in_bytes__bwt +
                                    bwt.size_in_bytes__C + bwt.size_in_bytes__char2int + bwt.size_in_bytes__Sigma);
    space_usage["stree_bp"]      = sdsl::size_in_bytes(bp);
    space_usage["stree_bpsupp"]  = (bp_supp.size_in_bytes_total + st.size_in_bytes__select + st.size_in_bytes__rank);
    cerr << "DONE (" << time_usage["dstruct"] / 1000 << " seconds)" << endl;

    /* compute RUNS */
    cerr << "computing runs ...";
    runs_start = timer::now();
    size_type k = t.size(), c = t[k - 1];
    IInterval I = init_interval(bwt, static_cast<char>(c)); //Interval<Bwt> I{bwt, static_cast<char>(c)};
    typedef typename StreeSada<t_bp_support>::node_type node_type;

    node_type v = st.wl(st.root(), c); // stree node

    while(--k > 0){
        c = t[k-1];
        I = bstep_interval(bwt, I, c); //I.bstep(c);
        if(I.first > I.second){
            runs[k] = 0;
            // update I to the parent of the proper locus of w until we can extend by 'c'
            do{
                v = st.parent(v);
                I.first = st.lb(v); I.second = st.rb(v); //I.set(st.lb(v), st.rb(v));
                I = bstep_interval(bwt, I, c); //I.bstep(c);
            } while(I.first > I.second);
        } else {
            runs[k] = 1;
        }
        v = st.wl(v, c); // update v
    }
    //size_type size_in_bytes_alg = build_runs_from_st_and_bwt(st, bwt, t, runs);
    runs_stop = timer::now();
    time_usage["alg"]  = std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count();
    space_usage["alg"] = 0;

    return std::make_pair(space_usage, time_usage);
}


performance_monitor build_ms_ohleb(const string& prefix, string& t, string& s_rev, bvector& runs, bvector& ms, const InputFlags& flags){
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

    auto get_ms = [] (sdsl::select_support_mcl<1,1>& __ms_select1, size_type __k) -> size_type {
        if(__k == -1)
            return (size_type) 1;
        return __ms_select1(__k + 1) - (2 * __k);
    };

    StreeOhleb<> st;

    cerr << "building the CST T(s') of lentgth " << s_rev.size() << "... ";
    auto runs_start = timer::now();
    sdsl::construct_im(st, s_rev, 1);
    auto runs_stop = timer::now();
    time_usage["dstruct"] = std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count();
    cerr << "DONE (" << time_usage["dstruct"] / 1000 << " seconds)" << endl;

    /* build MS */
    cerr << "building MS ... ";
    runs_start = timer::now();
    //size_type size_in_bytes_alg = build_ms_from_st_and_bwt(st, bwt, t, prefix, runs, ms, flags.lazy);
    size_type size_in_bytes_ms_select1 = 0;
    size_type k = 0, h_star = k + 1, k_prim, ms_idx = 0, ms_size = t.size() ;
    uint8_t c = t[k];
    IInterval I = init_interval(st, static_cast<char>(c)); //Interval I{bwt, static_cast<char>(c)};

    node_type v = st.wl(st.root(), c); // stree node
    while(k < ms_size){
        sdsl::select_support_mcl<1,1> ms_select1(&ms);
        size_in_bytes_ms_select1 = (size_in_bytes_ms_select1 < sdsl::size_in_bytes(ms_select1) ?
                                    sdsl::size_in_bytes(ms_select1) : size_in_bytes_ms_select1);

        if(flags.lazy){
            bool followed_wl = false;
            for(; I.first <= I.second && h_star < ms_size; ){
                c = t[h_star];
                I = bstep_interval(st, I, c); //I.bstep(c);
                if(I.first <= I.second){
                    v = st.lazy_wl(v, c);
                    followed_wl = true;
                    h_star++;
                }
            }
            if (followed_wl) // finish completing the new node
                st.lazy_wl_followup(v);
        } else { // non-lazy weiner links
            for(; I.first <= I.second && h_star < ms_size; ){
                c = t[h_star];
                I = bstep_interval(st, I, c); //I.bstep(c);
                if(I.first <= I.second){
                    v = st.wl(v, c);
                    h_star++;
                }
            }
        }

        for(int i = 0; i < h_star - k - get_ms(ms_select1, k - 1) + 1; i++) // probably not needed
            ms[ms_idx++] = 0;
        if(h_star - k - get_ms(ms_select1, k - 1) + 1 > 0)
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
            cerr << " (k, h*, k') = " << k << ", " << h_star << ", "  << k_prim << " : " << k_prim - k;
        }

        // update v
        v = st.wl(v, c);
        k = k_prim;
    }

    runs_stop = timer::now();
    time_usage["alg"]            = std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count();
    space_usage["stree_csa"]     = sdsl::size_in_bytes(st.csa);
    space_usage["stree_bp"]      = sdsl::size_in_bytes(st.bp);
    space_usage["stree_bpsupp"]  = sdsl::size_in_bytes(st.bp_support);
    space_usage["alg"]           = size_in_bytes_ms_select1;
    cerr << "DONE (" << time_usage["alg"] / 1000 << " seconds)" << endl;

    return std::make_pair(space_usage, time_usage);
}


performance_monitor build_runs_ohleb(const string& prefix, string& t, string& s, bvector& runs, const InputFlags& flags){
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
    cerr << "building the CST T(s) of lentgth " << s.size() << "... ";
    auto runs_start = timer::now();
    sdsl::construct_im(st, s, 1);
    auto runs_stop = timer::now();
    time_usage["dstruct"]       = std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count();
    cerr << "DONE (" << time_usage["dstruct"] / 1000 << " seconds)" << endl;
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
        if (flags.runs_progress > 0 && flags.runs_progress * k % t.size() == 0)
            report_progress(runs_start, t.size() - k, t.size());
    }
    //size_type size_in_bytes_alg = build_runs_from_st_and_bwt(st, bwt, t, runs);
    runs_stop = timer::now();
    time_usage["alg"]  = std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count();
    space_usage["alg"] = 0;
    cerr << "DONE (" << time_usage["alg"] / 1000 << " seconds)" << endl;
    return std::make_pair(space_usage, time_usage);
}

void comp(const string& prefix, InputSpec& T, InputSpec& S_fwd, const InputFlags& flags){
    monitor::size_dict space_usage, time_usage;
    performance_monitor runs_usage, ms_usage;

    string t = T.load_s();
    string s = S_fwd.load_s();
    bvector runs(t.size());
    bvector ms(t.size() * 2);

    if(flags.sada){
        runs_usage = build_runs_sada(prefix, t, s, runs, flags);
        reverse_in_place(s);
        ms_usage = build_ms_sada(prefix, t, s, runs, ms, flags);
    } else {
        runs_usage = build_runs_ohleb(prefix, t, s, runs, flags);
        reverse_in_place(s);
        ms_usage = build_ms_ohleb(prefix, t, s, runs, ms, flags);
    }

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
        cout << "prefix,len_s,len_t,measuring,unit,item,value" << endl;
        if(flags.space_usage){
            for(auto item: space_usage)
                cout << prefix << "," << s.size() << "," << t.size() << ",space,byte," << item.first << "," << item.second << endl;
        }
        if(flags.time_usage){
            for(auto item : runs_usage.second)
                cout << prefix << "," << s.size() << "," << t.size() << ",time,runs," << item.first << "," << item.second << endl;

            for(auto item : ms_usage.second)
                cout << prefix << "," << s.size() << "," << t.size() << ",time,ms," << item.first << "," << item.second << endl;
        }
    }

    if(flags.answer){
        cout << prefix << " ";
        dump_ms(ms);
    }

}


int main(int argc, char **argv){
    InputParser input(argc, argv);
    if(argc == 1){
        const string base_dir = {"/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/lazy_vs_nonlazy_data/input_data/"};
        string prefix {"acgt200000_100000"};
        InputFlags flags(false, // lazy_wl
                         false, // sada cst
                         false, // space
                         false, // time
                         true,  // ans
                         true,  // verbose
                         0,     //nr. progress messages for runs construction
                         5000      //nr. progress messages for ms construction
                         );
        InputSpec tspec(base_dir + prefix + "t.txt");
        InputSpec sfwd_spec(base_dir + prefix + "s.txt");
        comp(prefix, tspec, sfwd_spec, flags);
    } else {
        const string& base_dir = input.getCmdOption("-d");
        const string& prefix = input.getCmdOption("-p");
        InputFlags flags(input);
        InputSpec tspec(base_dir + "/" + prefix + "t.txt");
        InputSpec sfwd_spec(base_dir + "/" + prefix + "s.txt");
        comp(prefix, tspec, sfwd_spec, flags);
    }
    return 0;
}

