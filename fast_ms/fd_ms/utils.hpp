//
//  utils.h
//  fast_ms
//
//  Created by denas on 4/3/17.
//  Copyright © 2017 denas. All rights reserved.
//

#ifndef utils_h
#define utils_h

#include <iostream>
#include <string>
#include <vector>

#include "maxrep_construction.hpp"
#include "stree_sct3.hpp"

#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))


using namespace std;


namespace fdms{
    typedef uint8_t char_type;
    typedef unsigned long long size_type;
    typedef std::pair<size_type, size_type> Interval;
    typedef map<std::string, size_type> Counter;

    
    //Declare various wl strategies
    typedef StreeOhleb<>::node_type (StreeOhleb<>::*wl_method_t1) (const StreeOhleb<>::node_type& v, const StreeOhleb<>::char_type c) const;
    typedef StreeOhleb<>::node_type (StreeOhleb<>::*wl_method_t2) (const StreeOhleb<>::node_type& v, const StreeOhleb<>::char_type c, const bool is_max) const;
    typedef std::pair<StreeOhleb<>::size_type, StreeOhleb<>::size_type> (sdsl::bwt_of_csa_wt<sdsl::csa_wt<>>::*double_rank_method)(const StreeOhleb<>::size_type i, const StreeOhleb<>::size_type j, const StreeOhleb<>::char_type c)const;
    typedef StreeOhleb<>::node_type (StreeOhleb<>::*parent_seq_method) (const StreeOhleb<>::node_type& v, const StreeOhleb<>::char_type c) const;

    
    template<class flags_type>
    parent_seq_method get_parent_seq_method(flags_type &flags) {
        return (flags.lca_parents ? &StreeOhleb<>::maxrep_ancestor : &StreeOhleb<>::parent_sequence);
    }

    template<class flags_type>
    double_rank_method get_rank_method(flags_type &flags) {
        return (flags.rank_fail ?
                &sdsl::bwt_of_csa_wt<sdsl::csa_wt<>>::double_rank_and_fail :
                &sdsl::bwt_of_csa_wt<sdsl::csa_wt<>>::double_rank);
    }
    
    template<class flags_type>
    wl_method_t1 get_wl_method(flags_type &flags) {
        if(flags.lazy)
            return (flags.rank_fail ? &StreeOhleb<>::lazy_double_rank_fail_wl : &StreeOhleb<>::lazy_double_rank_wl);
        return (flags.rank_fail ? &StreeOhleb<>::double_rank_fail_wl : &StreeOhleb<>::double_rank_nofail_wl);
    }
    
    using timer = std::chrono::high_resolution_clock;

    void dump_ms(sdsl::bit_vector& ms){
        size_type k = 0;
        for (size_type i = 0; i < ms.size(); i++){
            if(ms[i] == 1){
                cout << i - (2*k) << " ";
                k += 1;
            }
        }
        //cout << endl;
    }

    void reverse_in_place(string& s){
        size_type n = s.size();

        for(int i = 0; i < n / 2; i++){
            char c = s[i];
            s[i] = s[n - 1 - i];
            s[n - 1 - i] = c;
        }
    }


    std::vector<pair<size_type, size_type>> slice_input(const size_type input_size, const size_type nthreads){
        size_type chunk = input_size / nthreads;
        size_type extra = input_size % nthreads;
        size_type step = 0;

        std::vector<pair<size_type, size_type>> slices (nthreads);
        for(size_type i=0, from = 0; i<nthreads; i++){
            step = chunk + (i < extra ? 1 : 0);
            slices[i] = std::make_pair(from, from + step);
            from += step;
        }
        return slices;
    }

    size_type load_st(StreeOhleb<>& st, const string &s, const string potential_stree_fname, const bool load){
        auto start = timer::now();
        if(load){
            cerr << " * loading the CST from " << potential_stree_fname << " ";
            sdsl::load_from_file(st, potential_stree_fname);
        } else {
            cerr << " * building the CST of length " << s.size() << " ";
            sdsl::construct_im(st, s, 1);
        }
        auto stop = timer::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    }


    /*
    size_type load_ohleb_st(StreeOhleb<> &st, const string &s, const string potential_stree_fname, const bool load){
        auto start = timer::now();
        if(load){
            cerr << " * loading the CST T(s') from " << potential_stree_fname << " ";
            sdsl::load_from_file(st, potential_stree_fname);
        } else {
            cerr << " * building the CST T(s') of length " << s.size() << " ";
            sdsl::construct_im(st, s, 1);
        }
        auto stop = timer::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    }
    */

    
    size_type load_maxrep(sdsl::bit_vector& vec, const StreeOhleb<>& st, const string& s, const string& load_fname, const bool load){
        auto start = timer::now();
        if(load){
            cerr << " * loading MAXREP from " << load_fname << " ";
            sdsl::load_from_file(vec, load_fname);
        } else {
            cerr << " * building MAXREP of length " << s.size() << " ";
            vec.resize(s.size() + 1);
            sdsl::util::set_to_value(vec, 0);
            build_maxrep_ohleb<StreeOhleb<>, sdsl::bit_vector, StreeOhleb<>::node_type>(st, vec);
        }
        auto stop = timer::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    }

    void report_progress(timer::time_point start_time, size_type curr_idx, size_type total){
        timer::time_point stop_time = timer::now();
        size_type elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time).count() + 1;
        cerr << endl << "[" << elapsed / 1000 << " s] " << 100.0 * curr_idx / total << "% @ " << (1.0 * curr_idx / elapsed) << " KHz";
    }
    
    sdsl::bit_vector parse_bitstr(string& s){
        sdsl::bit_vector b(s.size());
        
        for(size_t i = 0; i < s.size(); i++)
            b[i] = ((unsigned char)s[i] - 48);
        return b;
    }


}

#endif /* utils_h */
