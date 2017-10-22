//
//  matching_statistics.h
//  fast_ms
//
//  Created by denas on 10/21/17.
//  Copyright © 2017 denas. All rights reserved.
//

#ifndef matching_statistics_h
#define matching_statistics_h

#include <vector>
#include <string>
#include <iostream>
#include <future>
#include <thread>

#include <sdsl/int_vector.hpp>



class MStats{
    typedef StreeOhleb<> cst_t;
    typedef typename cst_t::node_type node_type;
    typedef sdsl::bit_vector vector_t;

private:
    const size_t nthreads;
    
public:
    MStats(const cst_t& st, sdsl::bit_vector runs, ){
        ;
    }
};


#endif /* matching_statistics_h */
