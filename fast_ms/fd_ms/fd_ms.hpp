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
#include <vector>

#include <sdsl/select_support.hpp>
#include "basic.hpp"
#include "stree_sct3.hpp"

using namespace std;

namespace fdms{
    typedef typename StreeOhleb<>::node_type node_type;

    Interval bstep_interval(StreeOhleb<>& st_, Interval& cur_i, char c){
        int cc = st_.csa.char2comp[c];
        return std::make_pair(st_.csa.C[cc] + st_.csa.bwt.rank(cur_i.first, c),
                              st_.csa.C[cc] + st_.csa.bwt.rank(cur_i.second + 1, c) - 1);
    }

    void resize_ms(bvector &ms_, float factor, size_type max_size){
        assert(factor > 1);
        size_type new_size = ms_.size() * factor;
        if(new_size > max_size)
            new_size = max_size;
        ms_.resize(new_size);
    }

    /* find k': index of the first zero to the right of k in runs */
    size_type find_k_prim_(size_type __k, size_type max__k, bvector& __runs){
        while(++__k < max__k && __runs[__k] != 0)
            ;
        return __k;
    }

    Interval fill_ms_slice_lazy(const string &t, StreeOhleb<> &st, bvector &ms, bvector &runs,
                                std::map<size_type, size_type> &consecutive_wl_calls,
                                const size_type from, const size_type to){
        size_type k = from, h_star = k + 1, h = h_star, h_star_prev = h_star, k_prim, ms_idx = 0, ms_size = t.size();
        uint8_t c = t[k];
        node_type v = st.wl(st.root(), c), u = v;

        while(k < to){
            h = h_star;
            h_star_prev = h_star;
            while((!st.is_root(u)) && h_star < ms_size){
                c = t[h_star];
                u = st.lazy_wl(v, c);
                if(!st.is_root(u)){
                    v = u;
                    h_star += 1;
                }
            }

            if(h_star > h_star_prev) // we must have called lazy_wl(). complete the node
                st.lazy_wl_followup(v);

            consecutive_wl_calls[(int) (h_star - h_star_prev)] += 1;

            while(ms_idx + (h_star - h) + 2 > ms.size()){
                resize_ms(ms, 1.5, t.size() * 2);
            }
            //ms_idx += (h_star -  h + 1);
            for(size_type i = 0; i < (h_star -  h + 1); i++)
                ms[ms_idx++] = 0; // adding 0s
            if(h_star - h + 1 > 0)
                ms[ms_idx++] = 1; // ... and a 1

            if(h_star < ms_size){ // remove prefixes of t[k..h*] until you can extend by 'c'
                do{ // remove suffixes of t[k..] until you can extend by 'c'
                    v = st.parent(v);
                    u = st.wl(v, t[h_star]);
                } while(st.is_root(u));

                h_star += 1;
            }
            // k_prim: index of the first zero to the right of k in runs
            k_prim = find_k_prim_(k, ms_size, runs);

            if(ms_idx + (k_prim - 1 - k) >= ms.size()){
                resize_ms(ms, 1.5, t.size() * 2);
            }
            for(size_type i = k + 1; i <= k_prim - 1 && i < to; i++)
                ms[ms_idx++] = 1;
            
            v = st.wl(v, c);
            k = k_prim;
        }

        pair<size_type, size_type> result(ms.size(), ms_idx);
        ms.resize(ms_idx);
        return result;
    }

    void check(const Interval I, const node_type u, const StreeOhleb<> &st){
        if(I.first > I.second)
            assert(st.is_root(u));
        else
            assert (I.first == u.i && I.second == u.j);
    }

    Interval fill_ms_slice_nonlazy(const string &t, StreeOhleb<> &st, bvector &ms, bvector &runs,
                                   std::map<size_type, size_type> &consecutive_wl_calls,
                                   const size_type from, const size_type to){
        size_type k = from, h_star = k + 1, h = h_star, h_star_prev = h_star, k_prim, ms_idx = 0, ms_size = t.size();
        uint8_t c = t[k];
        node_type v = st.wl(st.root(), c), u = v;

        while(k < to){
            h = h_star;
            h_star_prev = h_star;
            u = v;

            while((!st.is_root(u)) && h_star < ms_size){
                c = t[h_star];
                u = st.wl(v, c);
                if(!st.is_root(u)){
                    v = u;
                    h_star += 1;
                }
            }
            consecutive_wl_calls[(int) (h_star - h_star_prev)] += 1;

            while(ms_idx + (h_star - h) + 2 > ms.size()){
                resize_ms(ms, 1.5, t.size() * 2);
            }
            //ms_idx += (h_star -  h + 1);
            for(size_type i = 0; i < (h_star -  h + 1); i++)
                ms[ms_idx++] = 0; // adding 0s
            if(h_star - h + 1 > 0)
                ms[ms_idx++] = 1; // ... and a 1


            if(h_star < ms_size){ // remove prefixes of t[k..h*] until you can extend by 'c'
                do{ // remove suffixes of t[k..] until you can extend by 'c'
                    v = st.parent(v);
                    u = st.wl(v, t[h_star]);
                } while(st.is_root(u));
                h_star += 1;
            }
            v = st.wl(v, c);

            // k_prim: index of the first zero to the right of k in runs
            k_prim = find_k_prim_(k, ms_size, runs);

            if(ms_idx + (k_prim - 1 - k) >= ms.size()){
                resize_ms(ms, 1.5, t.size() * 2);
            }
            for(size_type i = k + 1; i <= k_prim - 1 && i < to; i++)
                ms[ms_idx++] = 1;

            k = k_prim;
        }

        pair<size_type, size_type> result(ms.size(), ms_idx);
        ms.resize(ms_idx);
        return result;
    }

}

#endif /* fd_ms_h */
