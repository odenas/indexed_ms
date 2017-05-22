//
//  fd_ms.hpp
//  fast_ms
//
//  Created by denas on 11/2/16.
//  Copyright © 2016 denas. All rights reserved.
//

#ifndef fd_ms_h
#define fd_ms_h

#include <iostream>
#include <string>

//#include "utils.hpp"

#include "maxrep_construction.hpp"
#include "runs_and_ms_algorithms.hpp"

using namespace std;

namespace fdms{
    Interval fill_ms_slice(const string& t, StreeOhleb<>& st, sdsl::bit_vector& ms, sdsl::bit_vector& runs, sdsl::bit_vector& maxrep,
                           const size_type from, const size_type to,
                           const bool lazy, const bool rank_fail, const bool use_maxrep){
        if(lazy){
            if(rank_fail){
                cerr << " *** computing with a lazy, double_rank_and_fail strategy ... " << endl;
                return fill_ms_slice_lazy_fail(t, st, ms, runs, from, to);
            } else {
                cerr << " *** computing with a lazy, double_rank_no_fail strategy ... " << endl;
                return fill_ms_slice_lazy_nofail(t, st, ms, runs, from, to);
            }
        } else { // nonlazy
            if(rank_fail){
                if(use_maxrep){
                    cerr << " *** computing with a non-lazy, double_rank_and_fail strategy, using maxrep ... " << endl;
                    return fill_ms_slice_nonlazy_fail(t, st, ms, runs, maxrep, from, to);
                } else {
                    cerr << " *** computing with a non-lazy, double_rank_and_fail strategy ... " << endl;
                    return fill_ms_slice_nonlazy_fail(t, st, ms, runs, from, to);
                }
            } else {
                cerr << " *** computing with a non-lazy, double_rank_no_fail strategy ... " << endl;
                return fill_ms_slice_nonlazy_nofail(t, st, ms, runs, from, to);
            }
        }
    }

    runs_rt fill_runs_slice(const string& t, StreeOhleb<>& st, sdsl::bit_vector& runs,
                            node_type v, const Interval slice,
                            const bool rank_fail){
        if(rank_fail){
            cerr << " *** computing with a double_rank_and_fail strategy ... " << endl;
            return fill_runs_slice_fail(t, st, runs, v, slice.first, slice.second);
        }
        cerr << " *** computing with a double_rank_no_fail strategy ... " << endl;
        return fill_runs_slice_nofail(t, st, runs, v, slice.first, slice.second);
    }
}

#endif /* fd_ms_h */
