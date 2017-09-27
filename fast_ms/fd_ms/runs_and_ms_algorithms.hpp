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

#define IS_INTNODE_MAXIMAL(node)                    ((maxrep[(node).i] == 1) && (maxrep[(node).j] == 1) )
#define IS_MAXIMAL(node) ( ((node).i != (node).j) && (maxrep[(node).i] == 1) && (maxrep[(node).j] == 1) )


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

    Interval bstep_on_interval(StreeOhleb<>& st_, Interval I, const int cc){
        I.first += st_.csa.C[cc];
        I.second += (st_.csa.C[cc] - 1);
        return I;
    }

    Interval init_interval(StreeOhleb<>& st_, const char c){
        int cc = st_.csa.char2comp[c];
        return std::make_pair(st_.csa.C[cc], st_.csa.C[cc + 1] - 1);
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


    runs_rt fill_runs_slice(const string &t, StreeOhleb<> &st,
                            wl_method_t1 wl_f_ptr, double_rank_method dr_f_ptr, parent_seq_method pseq_f_ptr,
                            sdsl::bit_vector &runs,
                            node_type v, const size_type from, const size_type to){
        size_type first_fail = 0, last_fail = 0;
        node_type last_fail_node = v;

        size_type k = to;
        char_type c = t[k - 1];
        bool idx_set = false;

        while(--k > from){
            c = t[k-1];

            if(!st.has_wl(v, c)){
                if(!idx_set){ // first failing wl()
                    first_fail = k;
                    idx_set = true;
                }
                runs[k] = 0;

                // remove suffixes of t[k..] until you can extend by 'c'
                v = CALL_MEMBER_FN(st, pseq_f_ptr)(v, c);

                // idx of last 0 in runs - 1 (within this block) and corresponding wl(node)
                last_fail_node = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);// given, parent_sequence() above, this has a wl()
                last_fail = k;
            } else {
                runs[k] = 1;
            }
            v = CALL_MEMBER_FN(st, wl_f_ptr)(v, c); // update v
        }
        if(!idx_set){
            first_fail = last_fail = from + 1;
            last_fail_node = v;
        }
        return make_tuple(first_fail, last_fail, last_fail_node);
    }

    void _set_next_ms_values1(sdsl::bit_vector& ms, size_type& ms_idx, const size_type h, const size_type h_star, const size_type max_ms_size){
        _resize_ms(ms, ms_idx + (h_star - h) + 2, max_ms_size);
        //ms_idx += (h_star -  h + 1);
        for(size_type i = 0; i < (h_star -  h + 1); i++)
            ms[ms_idx++] = 0; // adding 0s
        if(h_star - h + 1 > 0)
            ms[ms_idx++] = 1; // ... and a 1
    }

    size_type _set_next_ms_values2(sdsl::bit_vector& ms, sdsl::bit_vector& runs, size_type& ms_idx,
                                   const size_type k, const size_type to, const size_type max_ms_size){
        // k_prim: index of the first zero to the right of k in runs
        size_type k_prim = find_k_prim_(k, runs.size(), runs);
        _resize_ms(ms, ms_idx + (k_prim - 1 - k) + 1, max_ms_size);
        for(size_type i = k + 1; i <= k_prim - 1 && i < to; i++)
            ms[ms_idx++] = 1;
        return k_prim;
    }

    Interval fill_ms_slice_maxrep(const string& t, StreeOhleb<>& st,
                                  double_rank_method dr_f_ptr, parent_seq_method pseq_f_ptr,
                                  sdsl::bit_vector& ms, sdsl::bit_vector& runs, sdsl::bit_vector& maxrep,
                                  const size_type from, const size_type to){
        size_type k = from, h_star = k + 1, h = h_star, ms_idx = 0, ms_size = t.size();
        uint8_t c = t[k];
        node_type v = st.double_rank_fail_wl_mrep(st.root(), c, true), u = v;
        bool is_maximal = true;

        while(k < to){
            h = h_star;
            
            while(h_star < ms_size){
                c = t[h_star];
                is_maximal = IS_INTNODE_MAXIMAL(v);
                u = st.double_rank_fail_wl_mrep(v, c, is_maximal);
                if(!st.is_root(u)){
                    v = u;
                    h_star += 1;
                } else
                    break;
            }
            _set_next_ms_values1(ms, ms_idx, h, h_star, t.size() * 2);
            
            if(h_star < ms_size){ // remove prefixes of t[k..h*] until you can extend by 'c'
                is_maximal = false;
                bool has_wl = false;
                u = st.root();
                do{ // remove suffixes of t[k..] until you can extend by 'c'
                    v = st.parent(v);
                    if(!is_maximal)
                        is_maximal = IS_INTNODE_MAXIMAL(v); //since parent of a maximal is a maximal
                    if(is_maximal){
                        u = st.double_rank_fail_wl_mrep(v, c, is_maximal);
                        has_wl = !st.is_root(u);
                    }
                } while(!has_wl && !st.is_root(v)); // since !maximal => no wl
                h_star += 1;
            }
            k = _set_next_ms_values2(ms, runs, ms_idx, k, to, t.size() * 2);
            v = u;
        }
        pair<size_type, size_type> result(ms.size(), ms_idx);
        ms.resize(ms_idx);
        return result;
    }

    Interval fill_ms_slice(const string& t, StreeOhleb<>& st,
                           wl_method_t1 wl_f_ptr, double_rank_method dr_f_ptr, parent_seq_method pseq_f_ptr,
                           sdsl::bit_vector& ms, sdsl::bit_vector& runs,
                           const size_type from, const size_type to){
        size_type k = from, h_star = k + 1, h = h_star, ms_idx = 0, ms_size = t.size();
        uint8_t c = t[k];
        node_type v = st.double_rank_fail_wl(st.root(), c), u = v;

        while(k < to){
            h = h_star;
            
            while(h_star < ms_size){
                c = t[h_star];
                u = st.double_rank_fail_wl(v, c);
                if(!st.is_root(u)){
                    v = u;
                    h_star += 1;
                } else
                    break;
            }
            _set_next_ms_values1(ms, ms_idx, h, h_star, t.size() * 2);

            if(h_star < ms_size){ // remove prefixes of t[k..h*] until you can extend by 'c'
                bool has_wl = false;
                u = st.root();
                do{ // remove suffixes of t[k..] until you can extend by 'c'
                    v = st.parent(v);
                    u = st.double_rank_fail_wl(v, c);
                    has_wl = !st.is_root(u);
                } while(!has_wl && !st.is_root(v));
                h_star += 1;
            }
            k = _set_next_ms_values2(ms, runs, ms_idx, k, to, t.size() * 2);
            v = u;
        }
        pair<size_type, size_type> result(ms.size(), ms_idx);
        ms.resize(ms_idx);
        return result;
    }
}


/*
 this uses the st.has_wl() methods
 */
Interval fill_ms_slice_maxrep_(const string& t, StreeOhleb<>& st,
                               double_rank_method dr_f_ptr, parent_seq_method pseq_f_ptr,
                               sdsl::bit_vector& ms, sdsl::bit_vector& runs, sdsl::bit_vector& maxrep,
                               const size_type from, const size_type to){
    size_type k = from, h_star = k + 1, h = h_star, ms_idx = 0, ms_size = t.size();
    uint8_t c = t[k];
    node_type v = st.double_rank_fail_wl(st.root(), c);
    bool is_maximal;
    
    while(k < to){
        h = h_star;
        
        while(h_star < ms_size){
            c = t[h_star];
            is_maximal = IS_MAXIMAL(v);
            if(st.has_wl_mrep(v, c, is_maximal)){
                v = st.double_rank_fail_wl_mrep(v, c, is_maximal);
                h_star += 1;
            } else
                break;
        }
        _set_next_ms_values1(ms, ms_idx, h, h_star, t.size() * 2);
        
        
        if(h_star < ms_size){ // remove prefixes of t[k..h*] until you can extend by 'c'
            is_maximal = IS_MAXIMAL(v);
            do{ // remove suffixes of t[k..] until you can extend by 'c'
                v = st.parent(v);
                if(!is_maximal)
                    is_maximal = IS_MAXIMAL(v); //since parent of a maximal is a maximal
            } while((!is_maximal) || (!st.has_wl(v, c))); // since !maximal => no wl
            h_star += 1;
        }
        k = _set_next_ms_values2(ms, runs, ms_idx, k, to, t.size() * 2);
        v = st.double_rank_fail_wl_mrep(v, c, IS_MAXIMAL(v));
    }
    pair<size_type, size_type> result(ms.size(), ms_idx);
    ms.resize(ms_idx);
    return result;
}

/*
 this uses the st.has_wl() methods
 */
Interval fill_ms_slice_(const string& t, StreeOhleb<>& st,
                        wl_method_t1 wl_f_ptr, double_rank_method dr_f_ptr, parent_seq_method pseq_f_ptr,
                        sdsl::bit_vector& ms, sdsl::bit_vector& runs,
                        const size_type from, const size_type to){
    
    size_type k = from, h_star = k + 1, h = h_star, ms_idx = 0, ms_size = t.size();
    uint8_t c = t[k];
    node_type v = CALL_MEMBER_FN(st, wl_f_ptr)(st.root(), c);
    
    while(k < to){
        h = h_star;
        
        while(h_star < ms_size) {
            c = t[h_star];
            if(st.has_wl(v, c)){
                //v = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
                v = st.double_rank_fail_wl(v, c);
                h_star += 1;
            } else
                break;
        }
        
        _set_next_ms_values1(ms, ms_idx, h, h_star, t.size() * 2);
        
        if(h_star < ms_size){ // remove prefixes of t[k..h*] until you can extend by 'c'
            //v = CALL_MEMBER_FN(st, pseq_f_ptr)(v, c);
            // TODO: assuming not-lazy
            do{ // remove suffixes of t[k..] until you can extend by 'c'
                v = st.parent(v);
            } while(!st.has_wl(v, c) && !st.is_root(v)); // since !maximal => no wl
            h_star += 1;
        }
        k = _set_next_ms_values2(ms, runs, ms_idx, k, to, t.size() * 2);
        if(st.has_wl(v, c))
            v = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
    }
    pair<size_type, size_type> result(ms.size(), ms_idx);
    ms.resize(ms_idx);
    return result;
}


#endif /* runs_and_ms_algorithms_h */
