//
//  runs_and_ms_algorithms.h
//  fast_ms
//
//  Created by denas on 4/23/17.
//  Copyright © 2017 denas. All rights reserved.
//

#ifndef runs_and_ms_algorithms_h
#define runs_and_ms_algorithms_h

#include <iostream>
#include "utils.hpp"

using namespace fdms;

namespace fdms {
    typedef typename StreeOhleb<>::node_type node_type;
    typedef tuple<size_type, size_type, node_type> runs_rt;

    void output_partial_vec(const sdsl::bit_vector& v, size_type idx, const char name[], bool verbose){
        if (!verbose)
            return ;
        cout << name << ": ";
        for(size_type i = 0; i < (size_type)v.size(); i++)
            cout << v[i];
        cout << endl;
        for(size_type i = 0; i < (size_type)strlen(name) + 2; i++)
            cout << " ";
        for(size_type i = 0; i < (size_type)v.size(); i++)
            cout << (i == idx ? "*" : " ");
        cout << endl;
    }

    /* find k': index of the first zero to the right of k in runs */
    size_type find_k_prim_(size_type __k, size_type max__k, sdsl::bit_vector& __runs){
        while(++__k < max__k && __runs[__k] != 0)
            ;
        return __k;
    }

    void resize_ms__(sdsl::bit_vector &ms_, float factor, size_type max_size){
        assert(factor > 1);
        size_type new_size = ms_.size() * factor;
        if(new_size > max_size)
            new_size = max_size;
        ms_.resize(new_size);
    }

    void _resize_ms(sdsl::bit_vector &ms_, size_type min_size, size_type max_size){
        if(min_size > max_size)
            min_size = max_size;
        while(min_size > ms_.size()){
            size_type new_size = ms_.size() * 1.5;
            if(new_size > max_size)
                new_size = max_size;
            ms_.resize(new_size);
        }

    }

    runs_rt fill_runs_slice_fail(const string &t, StreeOhleb<> &st, sdsl::bit_vector &runs, node_type v,
                                 const size_type from, const size_type to){
        size_type first_fail = 0, last_fail = 0;
        node_type last_fail_node = v;

        size_type k = to, c = t[k - 1];
        bool idx_set = false;

        while(--k > from){
            c = t[k-1];
            if(st.is_root(st.double_rank_fail_wl(v, c))){ // empty
                if(!idx_set){ // first failing wl()
                    first_fail = k;
                    idx_set = true;
                }
                runs[k] = 0;
                do{ // remove suffixes of t[k..] until you can extend by 'c'
                    v = st.parent(v);
                } while(st.is_root(st.double_rank_fail_wl(v, c)));
                // idx of last 0 in runs - 1 (within this block) and corresponding wl(node)
                last_fail_node = st.double_rank_fail_wl(v, c);// given, parent_sequence() above, this has a wl()
                last_fail = k;
            } else {
                runs[k] = 1;
            }
            v = st.double_rank_fail_wl(v, c); // update v
        }
        if(!idx_set){
            first_fail = last_fail = from + 1;
            last_fail_node = v;
        }
        return make_tuple(first_fail, last_fail, last_fail_node);
    }

    runs_rt fill_runs_slice_nofail(const string &t, StreeOhleb<> &st, sdsl::bit_vector &runs, node_type v,
                                   const size_type from, const size_type to){
        size_type first_fail = 0, last_fail = 0;
        node_type last_fail_node = v;

        size_type k = to, c = t[k - 1];
        bool idx_set = false;

        while(--k > from){
            c = t[k-1];
            if(st.is_root(st.double_rank_nofail_wl(v, c))){ // empty
                if(!idx_set){ // first failing wl()
                    first_fail = k;
                    idx_set = true;
                }
                runs[k] = 0;
                do{ // remove suffixes of t[k..] until you can extend by 'c'
                    v = st.parent(v);
                } while(st.is_root(st.double_rank_nofail_wl(v, c)));
                // idx of last 0 in runs - 1 (within this block) and corresponding wl(node)
                last_fail_node = st.double_rank_nofail_wl(v, c);// given, parent_sequence() above, this has a wl()
                last_fail = k;
            } else {
                runs[k] = 1;
            }
            v = st.double_rank_nofail_wl(v, c); // update v
        }
        if(!idx_set){
            first_fail = last_fail = from + 1;
            last_fail_node = v;
        }
        return make_tuple(first_fail, last_fail, last_fail_node);
    }



    Interval fill_ms_slice_nonlazy_fail(const string &t, StreeOhleb<> &st, sdsl::bit_vector &ms, sdsl::bit_vector &runs,
                                          const size_type from, const size_type to){
        size_type k = from, h_star = k + 1, h = h_star, h_star_prev = h_star, k_prim, ms_idx = 0, ms_size = t.size();
        uint8_t c = t[k];
        node_type v = st.double_rank_fail_wl(st.root(), c), u = v;

        while(k < to){
            h = h_star;
            h_star_prev = h_star;
            u = v;

            while((!st.is_root(u)) && h_star < ms_size){
                c = t[h_star];
                u = st.double_rank_fail_wl(v, c);
                if(!st.is_root(u)){
                    v = u;
                    h_star += 1;
                }
            }

            _resize_ms(ms, ms_idx + (h_star - h) + 2, t.size() * 2);
            //ms_idx += (h_star -  h + 1);
            for(size_type i = 0; i < (h_star -  h + 1); i++)
                ms[ms_idx++] = 0; // adding 0s
            if(h_star - h + 1 > 0)
                ms[ms_idx++] = 1; // ... and a 1


            if(h_star < ms_size){ // remove prefixes of t[k..h*] until you can extend by 'c'
                do{ // remove suffixes of t[k..] until you can extend by 'c'
                    v = st.parent(v);
                    u = st.double_rank_fail_wl(v, c);
                } while(st.is_root(u));
                h_star += 1;
            }

            // k_prim: index of the first zero to the right of k in runs
            k_prim = find_k_prim_(k, ms_size, runs);

            _resize_ms(ms, ms_idx + (k_prim - 1 - k) + 1, t.size() * 2);
            for(size_type i = k + 1; i <= k_prim - 1 && i < to; i++)
                ms[ms_idx++] = 1;

            v = st.double_rank_fail_wl(v, c);
            k = k_prim;
        }
        pair<size_type, size_type> result(ms.size(), ms_idx);
        ms.resize(ms_idx);
        return result;
    }

    Interval fill_ms_slice_nonlazy_nofail(const string &t, StreeOhleb<> &st, sdsl::bit_vector &ms, sdsl::bit_vector &runs,
                                          const size_type from, const size_type to){
        size_type k = from, h_star = k + 1, h = h_star, h_star_prev = h_star, k_prim, ms_idx = 0, ms_size = t.size();
        uint8_t c = t[k];
        node_type v = st.double_rank_nofail_wl(st.root(), c), u = v;

        while(k < to){
            h = h_star;
            h_star_prev = h_star;
            u = v;

            while((!st.is_root(u)) && h_star < ms_size){
                c = t[h_star];
                u = st.double_rank_nofail_wl(v, c);
                if(!st.is_root(u)){
                    v = u;
                    h_star += 1;
                }
            }

            _resize_ms(ms, ms_idx + (h_star - h) + 2, t.size() * 2);
            //ms_idx += (h_star -  h + 1);
            for(size_type i = 0; i < (h_star -  h + 1); i++)
                ms[ms_idx++] = 0; // adding 0s
            if(h_star - h + 1 > 0)
                ms[ms_idx++] = 1; // ... and a 1


            if(h_star < ms_size){ // remove prefixes of t[k..h*] until you can extend by 'c'
                do{ // remove suffixes of t[k..] until you can extend by 'c'
                    v = st.parent(v);
                    u = st.double_rank_nofail_wl(v, c);
                } while(st.is_root(u));
                h_star += 1;
            }

            // k_prim: index of the first zero to the right of k in runs
            k_prim = find_k_prim_(k, ms_size, runs);

            _resize_ms(ms, ms_idx + (k_prim - 1 - k), t.size() * 2);
            for(size_type i = k + 1; i <= k_prim - 1 && i < to; i++)
                ms[ms_idx++] = 1;

            v = st.double_rank_nofail_wl(v, c);
            k = k_prim;
        }
        pair<size_type, size_type> result(ms.size(), ms_idx);
        ms.resize(ms_idx);
        return result;
    }

    Interval fill_ms_slice_lazy_nofail(const string &t, StreeOhleb<> &st, sdsl::bit_vector &ms, sdsl::bit_vector &runs,
                                       const size_type from, const size_type to){
        size_type k = from, h_star = k + 1, h = h_star, h_star_prev = h_star, k_prim, ms_idx = 0, ms_size = t.size();
        uint8_t c = t[k];
        node_type v = st.lazy_double_rank_wl(st.root(), c); st.lazy_wl_followup(v);
        node_type u = v;

        while(k < to){
            h = h_star;
            h_star_prev = h_star;
            while((!st.is_root(u)) && h_star < ms_size){
                c = t[h_star];
                u = st.lazy_double_rank_wl(v, c);
                if(!st.is_root(u)){
                    v = u;
                    h_star += 1;
                }
            }
            if(h_star > h_star_prev) // we must have called lazy_wl(). complete the node
                st.lazy_wl_followup(v);

            _resize_ms(ms, ms_idx + (h_star - h) + 2, t.size() * 2);
            //ms_idx += (h_star -  h + 1);
            for(size_type i = 0; i < (h_star -  h + 1); i++)
                ms[ms_idx++] = 0; // adding 0s
            if(h_star - h + 1 > 0)
                ms[ms_idx++] = 1; // ... and a 1

            if(h_star < ms_size){ // remove prefixes of t[k..h*] until you can extend by 'c'
                do{ // remove suffixes of t[k..] until you can extend by 'c'
                    v = st.parent(v);
                    u = st.lazy_double_rank_wl(v, c);
                } while(st.is_root(u));

                h_star += 1;
            }
            // k_prim: index of the first zero to the right of k in runs
            k_prim = find_k_prim_(k, ms_size, runs);

            _resize_ms(ms, ms_idx + (k_prim - 1 - k) + 1, t.size() * 2);
            for(size_type i = k + 1; i <= k_prim - 1 && i < to; i++)
                ms[ms_idx++] = 1;

            v = st.lazy_double_rank_wl(v, c); st.lazy_wl_followup(v);
            k = k_prim;
        }

        pair<size_type, size_type> result(ms.size(), ms_idx);
        ms.resize(ms_idx);
        return result;
    }

    Interval fill_ms_slice_lazy_fail(const string &t, StreeOhleb<> &st, sdsl::bit_vector &ms, sdsl::bit_vector &runs,
                                          const size_type from, const size_type to){
        size_type k = from, h_star = k + 1, h = h_star, h_star_prev = h_star, k_prim, ms_idx = 0, ms_size = t.size();
        uint8_t c = t[k];
        node_type v = st.lazy_double_rank_fail_wl(st.root(), c); //st.lazy_wl_followup(v);
        node_type u = v;

        while(k < to){
            h = h_star;
            h_star_prev = h_star;
            u = v;

            while((!st.is_root(u)) && h_star < ms_size){
                c = t[h_star];
                u = st.lazy_double_rank_fail_wl(v, c);
                if(!st.is_root(u)){
                    v = u;
                    h_star += 1;
                }
            }
            if(h_star > h_star_prev) // we must have called lazy_wl(). complete the node
                st.lazy_wl_followup(v);

            _resize_ms(ms, ms_idx + (h_star - h) + 2, t.size() * 2);
            //ms_idx += (h_star -  h + 1);
            for(size_type i = 0; i < (h_star -  h + 1); i++)
                ms[ms_idx++] = 0; // adding 0s
            if(h_star - h + 1 > 0)
                ms[ms_idx++] = 1; // ... and a 1


            if(h_star < ms_size){ // remove prefixes of t[k..h*] until you can extend by 'c'
                do{ // remove suffixes of t[k..] until you can extend by 'c'
                    v = st.parent(v);
                    u = st.lazy_double_rank_fail_wl(v, c);
                } while(st.is_root(u));
                h_star += 1;
            }
            
            // k_prim: index of the first zero to the right of k in runs
            k_prim = find_k_prim_(k, ms_size, runs);

            _resize_ms(ms, ms_idx + (k_prim - 1 - k) + 1, t.size() * 2);
            for(size_type i = k + 1; i <= k_prim - 1 && i < to; i++)
                ms[ms_idx++] = 1;

            v = st.lazy_double_rank_wl(v, c); st.lazy_wl_followup(v);
            k = k_prim;
        }
        
        pair<size_type, size_type> result(ms.size(), ms_idx);
        ms.resize(ms_idx);
        return result;
    }
}

#endif /* runs_and_ms_algorithms_h */
