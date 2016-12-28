//
//  fd_ms_algorithms.h
//  fast_ms
//
//  Created by denas on 12/12/16.
//  Copyright © 2016 denas. All rights reserved.
//

#ifndef fd_ms_algorithms_h
#define fd_ms_algorithms_h

#include "basic.hpp"
#include "Bwt.hpp"

#define LAZY_STREE
//#define VERBOSE

namespace fdms {


    void output_partial_vec(const bvector& v, size_type idx, const char name[], bool verbose){
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


    template<class t_cst, class bwt_t>
    size_type build_runs_from_st_and_bwt(const t_cst& st, bwt_t& bwt, const string& t, bvector& runs){
        size_type k = t.size(), c = t[k - 1];
        Interval I{bwt, static_cast<char>(c)};
        typedef typename t_cst::node_type node_type;

        node_type v = st.wl(st.root(), c); // stree node
        while(--k > 0){
            c = t[k-1];
            I.bstep(c);
            if(I.is_empty()){
                runs[k] = 0;
                // update I to the parent of the proper locus of w until we can extend by 'c'
                do{
                    v = st.parent(v);
                    I.set(st.lb(v), st.rb(v));
                    I.bstep(c);
                } while(I.is_empty());
            } else {
                runs[k] = 1;
            }
            v = st.wl(v, c); // update v
#ifdef VERBOSE
            output_partial_vec(runs, k, "runs", true);
#endif
        }
        return 0;
    }


    template<class t_cst, class bwt_t>
    size_type build_ms_from_st_and_bwt(t_cst& st, bwt_t& bwt, const string& t, const string& prefix, bvector& runs, bvector& ms, bool lazy_stree){

        auto get_ms = [] (sdsl::select_support_mcl<1,1>& __ms_select1, size_type __k) -> size_type {
            if(__k == -1)
                return (size_type) 1;
            return __ms_select1(__k + 1) - (2 * __k);
        };

        auto find_k_prim = []  (size_type __k, bvector& __runs, sdsl::rank_support_v<0>& __runs_rank0, sdsl::select_support_mcl<0, 1> __runs_sel0) -> size_type {
            size_t zeros = __runs_rank0(__k + 1);
            //return (__runs_rank0(__runs.size()) == zeros ? __runs.size() : __runs_sel0(zeros + 1));

            if(__runs_rank0(__runs.size()) == zeros)
                return __runs.size();
            return __runs_sel0(zeros + 1);
        };

        auto find_k_prim_ = [] (size_type __k, size_type max__k, bvector& __runs){
            while(++__k < max__k && __runs[__k] != 0)
                ;
            return __k;
        };

        typedef typename t_cst::node_type node_type;

        //TODO: these support data structures might not be needed
        //sdsl::rank_support_v<0> runs_rank0(&runs);
        //sdsl::select_support_mcl<0, 1> runs_select0(&runs);
        size_type size_in_bytes_ms_select1 = 0;


        size_type k = 0, h_star = k + 1, k_prim, ms_idx = 0, ms_size = t.size() ;
        uint8_t c = t[k];
        Interval I{bwt, static_cast<char>(c)};

        node_type v = st.wl(st.root(), c); // stree node
        while(k < ms_size){
#ifdef VERBOSE
            output_partial_vec(ms, ms_idx, "ms", true);
#endif

            sdsl::select_support_mcl<1,1> ms_select1(&ms);
            size_in_bytes_ms_select1 = (size_in_bytes_ms_select1 < sdsl::size_in_bytes(ms_select1) ? sdsl::size_in_bytes(ms_select1) : size_in_bytes_ms_select1);


            if(lazy_stree){
                bool followed_wl = false;

                for(; !I.is_empty() && h_star < ms_size; ){
                    c = t[h_star];
                    I.bstep(c);
                    if(!I.is_empty()){
                        v = st.lazy_wl(v, c);
                        followed_wl = true;
                        h_star++;
                    }
                }

                if (followed_wl) // finish completing the new node
                    st.lazy_wl_followup(v);
            } else { // non-lazy weiner links
                for(; !I.is_empty() && h_star < ms_size; ){
                    c = t[h_star];
                    I.bstep(c);
                    if(!I.is_empty()){
                        v = st.wl(v, c);
                        h_star++;
                    }
                }
            }

            for(int i = 0; i < h_star - k - get_ms(ms_select1, k - 1) + 1; i++)
                ms[ms_idx++] = 0;
            if(h_star - k - get_ms(ms_select1, k - 1) + 1 > 0)
                ms[ms_idx++] = 1;


            if(h_star < ms_size){
                do {
                    v = st.parent(v);
                    I.set(st.lb(v), st.rb(v));
                    I.bstep(t[h_star]);
                } while(I.is_empty());
                h_star = h_star + 1;
            }
            // k_prim: index of the first zero to the right of k in runs
            //k_prim = find_k_prim(k, runs, runs_rank0, runs_select0);
            //size_type __i = find_k_prim(k, runs, runs_rank0, runs_select0);
            k_prim = find_k_prim_(k, ms_size, runs);

            for(size_type i = k + 1; i <= k_prim - 1; i++)
                ms[ms_idx++] = 1;

            // update v
            v = st.wl(v, c);
            k = k_prim;
        }
        //cout << "bstep_counter = " << bstep_counter << ", wl_counter = " << wl_counter << endl;
        //return (sdsl::size_in_bytes(runs_rank0) + sdsl::size_in_bytes(runs_select0) + size_in_bytes_ms_select1);
        return size_in_bytes_ms_select1;
    }

}
#endif /* fd_ms_algorithms_h */
