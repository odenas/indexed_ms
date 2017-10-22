//
//  utils.h
//  fast_ms
//
//  Created by denas on 4/3/17.
//  Copyright © 2017 denas. All rights reserved.
//

#ifndef utils_h
#define utils_h

#include "stree_sct3.hpp"

#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))


using namespace std;


namespace fdms{
    typedef uint8_t char_type;
    typedef unsigned long long size_type;
    typedef std::pair<size_type, size_type> Interval;
    typedef map<std::string, size_type> Counter;
    using timer = std::chrono::high_resolution_clock;

    
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
