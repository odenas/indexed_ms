//
//  runs_and_ms_algorithms.h
//  fast_ms
//
//  Created by denas on 4/23/17.
//  Copyright © 2017 denas. All rights reserved.
//

#ifndef runs_and_ms_algorithms_h
#define runs_and_ms_algorithms_h

#include "maxrep_vector.hpp"
#include "parallel_runs_ms.hpp"


using namespace fdms;


typedef typename StreeOhleb<>::node_type node_type;
typedef typename StreeOhleb<>::size_type size_type;
typedef uint8_t char_type;
typedef map<std::string, size_type> Counter;
typedef std::pair<size_type, size_type> Interval;
typedef tuple<size_type, size_type, node_type> runs_rt;
typedef pair<size_type, node_type> ms_rt;

//Declare various wl strategies
typedef StreeOhleb<>::node_type (StreeOhleb<>::*wl_method_t1) (const StreeOhleb<>::node_type& v, const StreeOhleb<>::char_type c) const;
typedef StreeOhleb<>::node_type (StreeOhleb<>::*wl_method_t2) (const StreeOhleb<>::node_type& v, const StreeOhleb<>::char_type c, const bool is_max) const;
//typedef std::pair<StreeOhleb<>::size_type, StreeOhleb<>::size_type> (sdsl::bwt_of_csa_wt<sdsl::csa_wt<>>::*double_rank_method)(const StreeOhleb<>::size_type i, const StreeOhleb<>::size_type j, const StreeOhleb<>::char_type c)const;
//typedef StreeOhleb<>::node_type (StreeOhleb<>::*parent_seq_method) (const StreeOhleb<>::node_type& v, const StreeOhleb<>::char_type c) const;
#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))


template<typename msvec_t>
runs_rt fill_runs_slice(const string &t, StreeOhleb<> &st,
                        wl_method_t1 wl_f_ptr, const bool lca_parents,
                        msvec_t& runs_ms, node_type v, const Interval slice){
    
    size_type first_fail = 0, last_fail = 0, from = slice.first, to = slice.second;
    node_type last_fail_node = v, u = v;

    size_type k = to;
    char_type c = t[k - 1];
    bool idx_set = false;
    
    while(--k > from){
        c = t[k-1];

        u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
        if(st.is_root(u)){
            if(!idx_set){ // first failing wl()
                first_fail = k;
                idx_set = true;
            }
            runs_ms.runs[k] = 0;
            
            // remove suffixes of t[k..] until you can extend by 'c'
            if(lca_parents){
#ifdef VERBOSE
                cerr << "runs_lca" << endl;
#endif
                v = st.maxrep_ancestor(v, c);
                u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
            } else {
#ifdef VERBOSE
                cerr << "runs_pseq" << endl;
#endif
                if(!st.has_complete_info(v))
                    st.lazy_wl_followup(v);
                bool has_wl = false;
                u = st.root();
                do{ // remove suffixes of t[k..] until you can extend by 'c'
                    v = st.parent(v);
                    u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
                    has_wl = !st.is_root(u);
                } while(!has_wl && !st.is_root(v));
            }
            
            // idx of last 0 in runs - 1 (within this block) and corresponding wl(node)
            last_fail_node = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);// given, parent_sequence() above, this has a wl()
            last_fail = k;
        } else {
            runs_ms.runs[k] = 1;
        }
        v = CALL_MEMBER_FN(st, wl_f_ptr)(v, c); // update v
    }
    if(!idx_set){
        first_fail = last_fail = from + 1;
        last_fail_node = v;
    }
    return make_tuple(first_fail, last_fail, last_fail_node);
}

template<typename msvec_t>
Interval fill_ms_slice(const string& t, StreeOhleb<>& st,
                       wl_method_t1 wl_f_ptr, 
                       const bool lca_parents,
                       msvec_t& runs_ms, const size_type thread_id,
                       const node_type start_node, const Interval slice){
    size_type last_fail = 0;
    size_type from = slice.first, to = slice.second;
    size_type k = from, h_star = k + 1, h = h_star, ms_idx = 0, ms_size = t.size();
    uint8_t c = t[k];
    node_type v = CALL_MEMBER_FN(st, wl_f_ptr)(start_node, c), u = v, last_fail_node = v;

    while(k < to){
        h = h_star;

        while(h_star < ms_size){
            c = t[h_star];
            u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
            if(!st.is_root(u)){
                v = u;
                h_star += 1;
            } else{
                last_fail = h_star;
                break;
            }
        }
        runs_ms.set_next_ms_values1(thread_id, ms_idx, h, h_star, t.size() * 2);

        if(h_star < ms_size){ // remove prefixes of t[k..h*] until you can extend by 'c'
            if(lca_parents){
#ifdef VERBOSE
                cerr << "ms_lca" << endl;
#endif
                v = st.maxrep_ancestor(v, c);
                u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
            } else {
#ifdef VERBOSE
                cerr << "ms_pseq" << endl;
#endif
                if(!st.has_complete_info(v))
                    st.lazy_wl_followup(v);
                bool has_wl = false;
                u = st.root();
                do{ // remove suffixes of t[k..] until you can extend by 'c'
                    v = st.parent(v);
                    u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
                    has_wl = !st.is_root(u);
                } while(!has_wl && !st.is_root(v));

                last_fail = k + 1;
                last_fail_node = u;
            }
            h_star += 1;
        }
        k = runs_ms.set_next_ms_values2(thread_id, ms_idx, k, to, t.size() * 2);
        v = u;
    }
    pair<size_type, size_type> result(runs_ms.mses[thread_id].size(), ms_idx);
    runs_ms.mses[thread_id].resize(ms_idx);
    return result;
}

template<typename maxrep_t, typename msvec_t>
Interval fill_ms_slice_maxrep(const string& t, StreeOhleb<>& st,
                              wl_method_t2 wl_f_ptr,
                              msvec_t& runs_ms, maxrep_t& maxrep, const size_type thread_id,
                              const node_type start_node, const Interval slice){
    size_type from = slice.first, to = slice.second;
    size_type k = from, h_star = k + 1, h = h_star, ms_idx = 0, ms_size = t.size();
    uint8_t c = t[k];
    bool is_maximal = true;
    node_type v = CALL_MEMBER_FN(st, wl_f_ptr)(st.root(), c, is_maximal), u = v;
    
    while(k < to){
        h = h_star;
        
        while(h_star < ms_size){
            c = t[h_star];
            is_maximal = maxrep.is_intnode_maximal(v);
            u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c, is_maximal);
            if(!st.is_root(u)){
                v = u;
                h_star += 1;
            } else
                break;
        }
        runs_ms.set_next_ms_values1(thread_id, ms_idx, h, h_star, t.size() * 2);
        
        if(h_star < ms_size){ // remove prefixes of t[k..h*] until you can extend by 'c'
            is_maximal = false;
            bool has_wl = false;
            u = st.root();
            do{ // remove suffixes of t[k..] until you can extend by 'c'
                v = st.parent(v);
                if(!is_maximal){
                    is_maximal = maxrep.is_intnode_maximal(v);
                }
                if(is_maximal){
                    u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c, is_maximal);
                    has_wl = !st.is_root(u);
                }
            } while(!has_wl && !st.is_root(v)); // since !maximal => no wl
            h_star += 1;
        }
        k = runs_ms.set_next_ms_values2(thread_id, ms_idx, k, to, t.size() * 2);
        v = u;
    }
    pair<size_type, size_type> result(runs_ms.mses[thread_id].size(), ms_idx);
    runs_ms.mses[thread_id].resize(ms_idx);
    return result;
}

#endif /* runs_and_ms_algorithms_h */
