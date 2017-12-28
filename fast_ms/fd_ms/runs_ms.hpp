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

#include "stats.hpp"
#include "slices.hpp"
#include "maxrep_vector.hpp"

#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))


using namespace std;
using namespace fdms;


namespace fdms {
    template<typename cst_t, typename bitvec_type>
    class MsVectors{
        typedef typename cst_t::size_type       size_type;
        typedef typename cst_t::char_type       char_type;
        typedef typename cst_t::node_type       node_type;
        typedef Maxrep<cst_t, bitvec_type>      maxrep_t;
        typedef std::pair<size_type, size_type> interval_t;
        typedef Slices<size_type>               slice_t;
		//Declare various wl strategies
		typedef node_type (cst_t::*wl_method_t1) (const node_type& v, const char_type c) const;
		typedef node_type (cst_t::*wl_method_t2) (const node_type& v, const char_type c, const bool is_max) const;
		// Declare a parent-sequence strategy type
		typedef node_type (MsVectors<cst_t, bitvec_type>::*pseq_method_t)(const cst_t& st, wl_method_t1 wl_f_ptr, const node_type& v, const char_type c) const;


    private:
        /* find k': index of the first zero to the right of k in runs */
        static size_type find_k_prim_(size_type __k, size_type max__k, bitvec_type& __runs){
            while(++__k < max__k && __runs[__k] != 0)
                ;
            return __k;
        }

        void set_next_ms_values1(size_type& ms_idx, const size_type h, const size_type h_star){
            //ms_idx += (h_star -  h + 1);
            for(size_type i = 0; i < (h_star -  h + 1); i++)
                ms[ms_idx++] = 0; // adding 0s
            if(h_star - h + 1 > 0)
                ms[ms_idx++] = 1; // ... and a 1
        }
        
        size_type set_next_ms_values2(size_type& ms_idx, const size_type k){
            size_type k_prim = find_k_prim_(k, runs.size(), runs);
            for(size_type i = k + 1; i <= k_prim - 1 && i < ms.size(); i++)
                ms[ms_idx++] = 1;
            return k_prim;
        }

        
    public:
        bitvec_type runs, ms;
        size_type query_size;

        MsVectors(){
            query_size = 0;
            runs = bitvec_type(query_size);
            ms = bitvec_type(query_size * 2);
        }

        MsVectors(const size_type query_size){
            runs = bitvec_type(query_size);
            sdsl::util::set_to_value(runs, 0);
            ms = bitvec_type(query_size * 2);
            sdsl::util::set_to_value(ms, 0);
        }

        MsVectors(const MsVectors &mv) {
            runs = bitvec_type(mv.runs.size());
            for(size_type i=0; i<runs.size(); i++)
                runs[i] = mv.runs[i];

            ms = bitvec_type(mv.ms.size());
            for(size_type i=0; i<runs.size(); i++)
                ms[i] = mv.ms[i];
        }

        /*
         * call parent(v) in sequece until reaching a node u for which wl(u, c) exists
         */
        node_type parent_sequence(const cst_t& st, wl_method_t1 wl_f_ptr, const node_type& v, const char_type c) const {
            size_type cc = st.m_csa.char2comp[c];
            node_type vv = v, u = st.root();

            if(!st.has_complete_info(vv))
                st.lazy_wl_followup(vv);

            bool has_wl = false;
            do{ // remove suffixes of t[k..] until you can extend by 'c'
                vv = st.parent(vv);
				u = CALL_MEMBER_FN(st, wl_f_ptr)(vv, c);
                has_wl = !st.is_root(u);
            } while(!has_wl && !st.is_root(vv));
            return vv;
        }

        /*
         * find the ancestor u  of v s.t., wl(u, c) is not the root
         */
        node_type lca_parent(const cst_t& st, wl_method_t1 wl_f_ptr, const node_type& v, const char_type c) const {
            size_type cc = st.m_csa.char2comp[c];
            size_type cnt_c = st.m_csa.C[cc + 1] - st.m_csa.C[cc];

            size_type r = (st.m_csa.bwt.rank(v.j + 1, c) < cnt_c ? st.m_csa.bwt.select_at_dist(c, v.j, 1) : st.size());
            node_type p = r < st.size() ? st.lca(v, st.select_leaf(r + 1)) : st.root();

            if(p.i == v.i)
                return p;

            //index of last occurrence of c before position v.i
            size_type l = (st.m_csa.bwt.rank(v.i, c) > 0 ? st.m_csa.bwt.select_at_dist(c, v.i, 0) : 0);
            if(p.i > l)
                return p;
            node_type q = l > 0 ? st.lca(st.select_leaf(l + (l < st.size())), v) : st.root();
            
            // computing lca(p, q)
            node_type res = (q.j - q.i <= p.j - p.i ? q : p);
            //node_type exp_res = _maxrep_ancestor(v, c);
            //if(res.i != exp_res.i or res.j != exp_res.j)
            //    assert(0);
            //assert (res == _maxrep_ancestor(v, c));
            return res;
		}

		/* runs and ms algorithms for input stats */
        void fill_runs(Stats<cst_t, maxrep_t>& stats, const string& t, const cst_t& st,
						wl_method_t1 wl_f_ptr, pseq_method_t pseq_f_ptr){

			bitvec_type v_{0};
			NodeProperty<cst_t, maxrep_t>NP{st, maxrep_t(v_)};
			size_type k = t.size();
			char_type c = t[k - 1];
			node_type v = CALL_MEMBER_FN(st, wl_f_ptr)(st.root(), c),
					  u = v;

		    while(--k > 0){
		        c = t[k-1];

				u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
				stats.runs_wl_calls[NP.runs_node_label(v, c)] += 1;
		        if(st.is_root(u)){
		            runs[k] = 0;
					u = v; // needed for register_runs_pseq
		        	v = CALL_MEMBER_FN(*this, pseq_f_ptr)(st, wl_f_ptr, v, c);
		        	// register consecutive parent calls and wl() calls therein
		        	stats.register_runs_pseq(st, NP, u, v, c);

					v = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
					stats.runs_wl_calls[NP.runs_node_label(v, c)] += 1;
		        } else {
		            runs[k] = 1;
					v = u;
		        }
		    }
		}

		void fill_ms(Stats<cst_t, maxrep_t>& stats, const string& t, const cst_t& st,
					 wl_method_t1 wl_f_ptr, pseq_method_t pseq_f_ptr,
					 const maxrep_t& maxrep){

			NodeProperty<cst_t, maxrep_t>NP{st, maxrep};
		    size_type k = 0, h_star = k + 1, h = h_star, h_star_prev = h_star, ms_idx = 0;
		    char_type c = t[k];
		    node_type v = st.double_rank_nofail_wl(st.root(), c), u = v;

		    while(k < t.size()){
		        h = h_star;
        
		        h_star_prev = h_star;
		        while(h_star < ms.size()){
		            c = t[h_star];
					u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
					stats.ms_wl_calls[NP.ms_node_label(v,  c)] += 1;
		            if(!st.is_root(u)){
		                v = u;
		                h_star += 1;
		            } else
		                break;
		        }
				stats.ms_wlcalls_seq[(size_type) (h_star - h_star_prev)] += 1;  // record
		        set_next_ms_values1(ms_idx, h, h_star);

		        if(h_star < t.size()){ // remove prefixes of t[k..h*] until you can extend by 'c'
					u = v;
		        	v = CALL_MEMBER_FN(*this, pseq_f_ptr)(st, wl_f_ptr, v, c);
		        	stats.register_ms_pseq(st, NP, u, v, c);
		            h_star += 1;
		        }
		        k = set_next_ms_values2(ms_idx, k);
				v = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
				stats.ms_wl_calls[NP.ms_node_label(v,  c)] += 1;
		    }
		}

		/* runs and ms algorithms */
        void fill_runs(const string& t, const cst_t& st, wl_method_t1 wl_f_ptr, pseq_method_t pseq_f_ptr){

			size_type k = t.size();
			char_type c = t[k - 1];
			node_type v = CALL_MEMBER_FN(st, wl_f_ptr)(st.root(), c),
					  u = v;

		    while(--k > 0){
		        c = t[k-1];

				u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
		        if(st.is_root(u)){
		            runs[k] = 0;
		        	v = CALL_MEMBER_FN(*this, pseq_f_ptr)(st, wl_f_ptr, v, c);
					v = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
		        } else {
		            runs[k] = 1;
					v = u;
		        }
		    }
		}

		void fill_ms(const string& t, const cst_t& st, wl_method_t1 wl_f_ptr, pseq_method_t pseq_f_ptr){

		    size_type k = 0, h_star = k + 1, h = h_star, ms_idx = 0;
		    char_type c = t[k];
		    node_type v = CALL_MEMBER_FN(st, wl_f_ptr)(st.root(), c), u = v;
		    bool is_maximal = true;

		    while(k < t.size()){
		        h = h_star;
        
		        while(h_star < runs.size()){
		            c = t[h_star];
					u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
		            if(!st.is_root(u)){
		                v = u;
		                h_star += 1;
		            } else
		                break;
		        }
		        set_next_ms_values1(ms_idx, h, h_star);

		        if(h_star < t.size()){ // remove prefixes of t[k..h*] until you can extend by 'c'
		        	v = CALL_MEMBER_FN(*this, pseq_f_ptr)(st, wl_f_ptr, v, c);
		            h_star += 1;
		        }
		        k = set_next_ms_values2(ms_idx, k);
				v = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
		    }
		}

		void fill_ms(const string& t, const cst_t& st, wl_method_t2 wl_f_ptr, maxrep_t& maxrep){
            cerr << " ** using maxrep " << endl;

		    size_type k = 0, h_star = k + 1, h = h_star, h_star_prev = h_star, ms_idx = 0;
		    char_type c = t[k];
		    node_type v = st.double_rank_nofail_wl(st.root(), c), u = v;
		    bool is_maximal = true;

		    while(k < t.size()){
		        h = h_star;
        
		        h_star_prev = h_star;
		        while(h_star < ms.size()){
		            c = t[h_star];
					is_maximal = maxrep.is_maximal(v);
					u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c, is_maximal);
		            if(!st.is_root(u)){
		                v = u;
		                h_star += 1;
		            } else
		                break;
		        }
		        set_next_ms_values1(ms_idx, h, h_star);

		        if(h_star < runs.size()){ // remove prefixes of t[k..h*] until you can extend by 'c'
		            //is_maximal = false;
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
		        k = set_next_ms_values2(ms_idx, k);
				v = u;
		    }
		}

        void show_runs(std::ostream& out){
            for(size_type i = 0; i < runs.size(); i++)
                out << runs[i] << " ";
            out << endl;
        }

        void show_MS(std::ostream& out){
        	size_type k = 0;
            for (size_type i = 0; i < ms.size(); i++){
                if(ms[i] == 1){
                    out << i - (2*k) << " ";
                    k += 1;
                }
			}
        }

        pair<size_type, size_type> ms_composition() const {
            size_type ones = 0;
            for (size_type i = 0; i < ms.size(); i++){
                if(ms[i] == 1)
                    ones += 1;
            }
            return std::make_pair<size_type, size_type>(ms.size() - ones, ones);
        }

        pair<size_type, size_type> runs_composition() const {
            size_type ones = 0;
            for (size_type i = 0; i < runs.size(); i++){
                if(runs[i] == 1)
                    ones += 1;
            }
            return std::make_pair<size_type, size_type>(runs.size() - ones, ones);
        }
    };
}

#endif /* runs_ms_h */
