/*
 * Represents the runs vector.
 */


#ifndef RUNS_VECTOR_HPP
#define RUNS_VECTOR_HPP


#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <sdsl/vectors.hpp>

#include "input_spec.hpp"
#include "stats.hpp"
#include "query.hpp"
#include "maxrep_vector.hpp"
#include "stats.hpp"

#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))

namespace fdms {

    template<typename cst_t>
    class runs_vector {
        typedef typename cst_t::size_type size_type;
        typedef typename cst_t::char_type char_type;
        typedef typename cst_t::node_type node_type;

        typedef Maxrep<cst_t, sdsl::bit_vector> maxrep_t;
        typedef sdsl::int_vector_buffer<1> buff_vec_t;

    public:
        // strategies for rank operations
        typedef node_type(cst_t::*wl_method_t1) (const node_type& v, const char_type c) const;
        typedef node_type(cst_t::*wl_method_t2) (const node_type& v, const char_type c, const bool is_max) const;

        // strategies for sequences of parent operations
        typedef node_type(*pseq_method_t) (const cst_t& st, wl_method_t1 wl_f_ptr, const node_type& v, const char_type c);


        /*
         * call parent(v) in sequece until reaching a node u for which wl(u, c) exists
         */
        static node_type parent_sequence(const cst_t& st, wl_method_t1 wl_f_ptr, const node_type& v, const char_type c) {
            node_type vv = v, u = st.root();

            if (!st.has_complete_info(vv))
                st.lazy_wl_followup(vv);

            bool has_wl = false;
            do { // remove suffixes of t[k..] until you can extend by 'c'
                vv = st.parent(vv);
                u = CALL_MEMBER_FN(st, wl_f_ptr)(vv, c);
                has_wl = !st.is_root(u);
            } while (!has_wl && !st.is_root(vv));
            return vv;
        }

        /*
         * find the ancestor u  of v s.t., wl(u, c) is not the root
         */
        static node_type lca_parent(const cst_t& st, wl_method_t1 wl_f_ptr, const node_type& v, const char_type c) {
            size_type cc = st.m_csa.char2comp[c];
            size_type cnt_c = st.m_csa.C[cc + 1] - st.m_csa.C[cc];

            size_type r = (st.m_csa.bwt.rank(v.j + 1, c) < cnt_c ? st.m_csa.bwt.select_at_dist(c, v.j, 1) : st.size());
            node_type p = r < st.size() ? st.lca(v, st.select_leaf(r + 1)) : st.root();

            if (p.i == v.i)
                return p;

            //index of last occurrence of c before position v.i
            size_type l = (st.m_csa.bwt.rank(v.i, c) > 0 ? st.m_csa.bwt.select_at_dist(c, v.i, 0) : 0);
            if (p.i > l)
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

        static void dump(const InputSpec ispec, const cst_t& st,
                wl_method_t1 wl_f_ptr, pseq_method_t pseq_f_ptr, const size_t buffer_size) {


            Query_rev t{ispec.t_fname, buffer_size};
            buff_vec_t runs(ispec.runs_fname, std::ios::out, (uint64_t) buffer_size);
            cerr << " ** dumping (buffersize: " << runs.buffersize() << ") to " << ispec.runs_fname << " ... ";

            size_type k = t.size();
            char_type c = t[k - 1];
            node_type v = CALL_MEMBER_FN(st, wl_f_ptr)(st.root(), c), u = v;
            while (--k > 0) {
                c = t[k - 1];

                u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
                if (st.is_root(u)) {
                    runs[k] = 0;
                    //v = CALL_MEMBER_FN(*this, pseq_f_ptr)(st, wl_f_ptr, v, c);
                    //v = parent_sequence(st, wl_f_ptr, v, c);
                    v = pseq_f_ptr(st, wl_f_ptr, v, c);
                    v = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
                } else {
                    runs[k] = 1;
                    v = u;
                }
            }
        }

        static void dump(const InputSpec ispec, const cst_t& st,
                wl_method_t2 wl_f_ptr, maxrep_t& maxrep, const size_t buffer_size) {

            Query_rev t{ispec.t_fname, buffer_size};
            cerr << " ** using maxrep (runs) " << endl;
            buff_vec_t runs(ispec.runs_fname, std::ios::out, buffer_size);
            cerr << " ** dumping (buffersize: " << runs.buffersize() << ") to " << ispec.runs_fname << " ... ";

            size_type k = t.size();
            char_type c = t[k - 1];
            node_type v = st.double_rank_nofail_wl(st.root(), c), u = v;
            bool is_maximal = true;

            while (--k > 0) {
                c = t[k - 1];
                is_maximal = maxrep.is_maximal(v);

                u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c, is_maximal);
                if (st.is_root(u)) {
                    runs[k] = 0;
                    bool has_wl = false;
                    do {
                        v = st.parent(v);
                        if (!is_maximal) {
                            is_maximal = maxrep.is_intnode_maximal(v);
                        }
                        if (is_maximal) {
                            u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c, is_maximal);
                            has_wl = !st.is_root(u);
                        }
                    } while (!has_wl && !st.is_root(v)); // since !maximal => no wl
                } else {
                    runs[k] = 1;
                }
                v = u;
            }
        }

        /* runs algorithm for input stats */
        static void fill_runs(Stats<cst_t, maxrep_t>& stats, const InputSpec ispec, const cst_t& st,
                wl_method_t1 wl_f_ptr, pseq_method_t pseq_f_ptr, const size_t buffer_size) {

            Query_rev t{ispec.t_fname, buffer_size};
            sdsl::bit_vector v_{0};
            NodeProperty<cst_t, maxrep_t>NP{st, maxrep_t(v_)};
            buff_vec_t runs(ispec.runs_fname, std::ios::out, buffer_size);
            size_type k = t.size();
            char_type c = t[k - 1];
            node_type v = CALL_MEMBER_FN(st, wl_f_ptr)(st.root(), c), u = v;

            while (--k > 0) {
                c = t[k - 1];

                u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
                stats.runs_wl_calls[NP.runs_node_label(v, c)] += 1;
                if (st.is_root(u)) {
                    runs[k] = 0;
                    u = v; // needed for register_runs_pseq
                    v = pseq_f_ptr(st, wl_f_ptr, v, c);
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


        static void show(const string runs_fname, std::ostream& out) {
            buff_vec_t runs(runs_fname, std::ios::in);
            for (size_type i = 0; i < runs.size(); i++)
                out << runs[i] << " ";
            out << endl;
        }

        static std::pair<size_type, size_type> composition(const string runs_fname) {
            buff_vec_t runs(runs_fname, std::ios::in);

            size_type ones = 0;
            for (size_type i = 0; i < runs.size(); i++) {
                if (runs[i] == 1)
                    ones += 1;
            }
            return std::make_pair<size_type, size_type>(runs.size() - ones, ones);
        }
        
        static size_type size(const InputSpec ispec) {
            buff_vec_t ms(ispec.runs_fname, std::ios::in);
            return (size_type) ms.size();
        }

    };
}

#endif /* RUNS_VECTOR_HPP */

