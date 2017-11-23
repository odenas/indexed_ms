//
//  runs_ms.hpp
//  fast_ms
//
//  Created by denas on 10/27/17.
//  Copyright © 2017 denas. All rights reserved.
//

#ifndef runs_ms_h
#define runs_ms_h

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

//#include "stree_sct3.hpp"
#include "slices.hpp"

using namespace std;
using namespace fdms;

namespace fdms {
    template<typename cst_t, typename bitvec_type>
    class MsVectors{
        typedef typename cst_t::size_type       size_type;
        typedef typename cst_t::node_type       node_type;
        typedef std::pair<size_type, size_type> interval_t;
        typedef Slices<size_type>               slice_t;

    private:
        static void _resize_ms(bitvec_type &ms_, size_type min_size, size_type max_size){
            if(min_size > max_size)
                min_size = max_size;
            while(min_size > ms_.size()){
                size_type new_size = ms_.size() * 1.5;
                if(new_size > max_size)
                    new_size = max_size;
                ms_.resize(new_size);
            }
            
        }
        
        static void _set_next_ms_values1(bitvec_type& ms, size_type& ms_idx,
                                         const size_type h, const size_type h_star, const size_type max_ms_size){
            _resize_ms(ms, ms_idx + (h_star - h) + 2, max_ms_size);
            //ms_idx += (h_star -  h + 1);
            for(size_type i = 0; i < (h_star -  h + 1); i++)
                ms[ms_idx++] = 0; // adding 0s
            if(h_star - h + 1 > 0)
                ms[ms_idx++] = 1; // ... and a 1
        }
        
        static size_type _set_next_ms_values2(bitvec_type& ms, bitvec_type& runs, size_type& ms_idx,
                                              const size_type k, const size_type to, const size_type max_ms_size){
            // k_prim: index of the first zero to the right of k in runs
            size_type k_prim = find_k_prim_(k, runs.size(), runs);
            _resize_ms(ms, ms_idx + (k_prim - 1 - k) + 1, max_ms_size);
            for(size_type i = k + 1; i <= k_prim - 1 && i < to; i++)
                ms[ms_idx++] = 1;
            return k_prim;
        }

        
        /* find k': index of the first zero to the right of k in runs */
        static size_type find_k_prim_(size_type __k, size_type max__k, bitvec_type& __runs){
            while(++__k < max__k && __runs[__k] != 0)
                ;
            return __k;
        }


    public:
        //typedef tuple<size_type, size_type, node_type> runs_rt;

        bitvec_type runs;
        vector<bitvec_type> mses; // the ms vector for each thread
        size_t nthreads;
        slice_t slices;

        MsVectors(){
            runs = bitvec_type(0);
            mses = vector<bitvec_type>(0);
            slices = slice_t();
            nthreads = 0;
        }

        MsVectors(const size_type query_size,  size_type const nthr){
            nthreads = nthr;
            runs = bitvec_type(query_size);
            mses = vector<bitvec_type>(nthreads);
            slices = slice_t(query_size, nthreads);
            for(int i=0; i<nthreads; i++){
                mses[i].resize(query_size / nthreads);
                sdsl::util::set_to_value(mses[i], 0);
            }
        }

        MsVectors(const MsVectors &mv) {
            nthreads = mv.nthreads;
            runs = bitvec_type(mv.runs.size());
            mses = vector<bitvec_type>(mv.mses.size());
            slices = mv.slices;

            for(size_type i=0; i<runs.size(); i++)
                runs[i] = mv.runs[i];
            
            for(size_type vi=0; vi < mses.size(); vi++){
                mses[vi].resize(mv.mses[vi].size());
                for(size_type i=0; i<runs.size(); i++)
                    mses[vi][i] = mv.mses[vi][i];
            }
        }
        
        void set_next_ms_values1(const size_type mses_idx, size_type& ms_idx,
                                 const size_type h, const size_type h_star, const size_type max_ms_size){
            _set_next_ms_values1(mses[mses_idx], ms_idx, h, h_star, max_ms_size);
        }
        
        size_type set_next_ms_values2(const size_type mses_idx, size_type& ms_idx,
                                      const size_type k, const size_type to, const size_type max_ms_size){
            return _set_next_ms_values2(mses[mses_idx], runs, ms_idx, k, to, max_ms_size);
        }
        
        size_type ms_size() const {
            size_type total_ms_length = 0;
            for(auto v: mses)
                total_ms_length += v.size();
            return total_ms_length;
        }
        
        void show_runs(std::ostream& out){
            for(size_type i = 0; i < runs.size(); i++)
                out << runs[i] << " ";
            out << endl;
        }

        void show_MS(std::ostream& out){
            for(auto ms : mses){
                size_type k = 0;
                for (size_type i = 0; i < ms.size(); i++){
                    if(ms[i] == 1){
                        out << i - (2*k) << " ";
                        k += 1;
                    }
                }
            }
        }
    };
}

#endif /* runs_ms_h */
