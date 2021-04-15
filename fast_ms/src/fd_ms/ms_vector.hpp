/*
 * Represents the ms vector.
 */

#ifndef MS_VECTOR_HPP
#define MS_VECTOR_HPP


#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <sdsl/vectors.hpp>

#include "input_spec.hpp"
#include "stats.hpp"
#include "query.hpp"
#include "maxrep_vector.hpp"


#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))

namespace fdms {

    template<typename cst_t>
    class ms_vector {
        typedef typename cst_t::size_type size_type;
        typedef typename cst_t::char_type char_type;
        typedef typename cst_t::node_type node_type;

        typedef sdsl::int_vector_buffer<1> buff_vec_t;
        typedef Maxrep<cst_t, sdsl::bit_vector> maxrep_t;

    public:
        // strategies for rank operations
        typedef node_type(cst_t::*wl_method_t1) (const node_type& v, const char_type c) const;
        typedef node_type(cst_t::*wl_method_t2) (const node_type& v, const char_type c, const bool is_max) const;
        // strategies for sequences of parent operations
        typedef node_type(*pseq_method_t) (const cst_t& st, wl_method_t1 wl_f_ptr, const node_type& v, const char_type c);

    private:

        /* find k': index of the first zero to the right of k in runs */
        static size_type find_k_prim_(size_type __k, const size_type max__k, buff_vec_t& __runs) {
            while (++__k < max__k && __runs[__k] != 0)
                ;
            return __k;
        }

        static void set_next_ms_values1(buff_vec_t& ms, size_type& ms_idx, const size_type h, const size_type h_star) {
            //ms_idx += (h_star -  h + 1);
            for (size_type i = 0; i < (h_star - h + 1); i++)
                ms[ms_idx++] = 0; // adding 0s
            if (h_star - h + 1 > 0)
                ms[ms_idx++] = 1; // ... and a 1
        }

        static size_type set_next_ms_values2(buff_vec_t& runs, buff_vec_t& ms, size_type& ms_idx, const size_type k) {
            size_type k_prim = find_k_prim_(k, runs.size(), runs);
            for (size_type i = k + 1; i <= k_prim - 1; i++)
                ms[ms_idx++] = 1;
            return k_prim;
        }

    public:

        /*
         * call parent(v) in sequece until reaching a node u for which wl(u, c) exists
         */
        static node_type parent_sequence(const cst_t& st, wl_method_t1 wl_f_ptr, const node_type& v, const char_type c) {
            size_type cc = st.m_csa.char2comp[c];
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

            Query_fwd t{ispec.t_fname, buffer_size};
            //sdsl::bit_vector runs(t.size());
            //sdsl::load_from_file(runs, ispec.runs_fname);
            //sdsl::bit_vector ms(2*t.size());

            //buff_vec_t frequency(ispec.runs_fname + ".freq", std::ios::out, (uint64_t) buffer_size);

            buff_vec_t runs(ispec.runs_fname, std::ios::in, buffer_size);
            buff_vec_t ms(ispec.ms_fname, std::ios::out, buffer_size);

            // assuming t[0] is in the index
            size_type k = 0, h_star = k + 1, h = h_star, ms_idx = 0;
            char_type c = t[k];
            node_type v = CALL_MEMBER_FN(st, wl_f_ptr)(st.root(), c), u = v;

            while (k < t.size()) {
                h = h_star;

                while (h_star < runs.size()) {
                    c = t[h_star];
                    u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
                    if (!st.is_root(u)) {
                        v = u;
                        h_star += 1;
                    } else
                        break;
                }
                //set_next_ms_values1(ms, ms_idx, h, h_star);
                {// add h* - h + 1 zeros and a 1 (if h*-h+1 > 0)
                    for (size_type i = 0; i < (h_star - h + 1); i++)
                        ms[ms_idx++] = 0; // adding 0s
                    if (h_star - h + 1 > 0)
                        ms[ms_idx++] = 1; // ... and a 1
                }
                {
                        if (!st.has_complete_info(v))
                            st.lazy_wl_followup(v);

                        bool has_wl = false;
                        size_type fk = k;
                        size_type k_prim = 0;
                        while(true){ // remove suffixes of t[k..] until you can extend by 'c'
                            {
                                k_prim = find_k_prim_(fk, runs.size(), runs);
                                assert(k_prim == runs.size() or runs[k_prim] == 0);
                                for(size_type i = fk; i < k_prim; i++)
                                    cout << i << "," << st.size(v) << endl;
                                fk = k_prim;
                            }

                            v = st.parent(v);
                            // remove prefixes of t[k..h*] until you can extend by 'c'
                            has_wl = (!st.is_root(CALL_MEMBER_FN(st, wl_f_ptr)(v, c)) and
                                      h_star < t.size());
                            if(has_wl or st.is_root(v)){
                                for (size_type i = k + 1; i <= k_prim - 1; i++)
                                    ms[ms_idx++] = 1;
                                k = k_prim;
                                break;
                            }
                        }
                    h_star += 1;
                }
                //k = set_next_ms_values2(runs, ms, ms_idx, k);
                v = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
            }
        }

        static void dump(const InputSpec ispec, const cst_t& st, wl_method_t2 wl_f_ptr, maxrep_t& maxrep,
                const size_t buffer_size) {
            cerr << " ** using maxrep (ms) " << endl;

            Query_fwd t{ispec.t_fname, buffer_size};
            //sdsl::bit_vector runs(t.size());
            //sdsl::load_from_file(runs, ispec.runs_fname);
            //sdsl::bit_vector ms(2*t.size());
            buff_vec_t runs(ispec.runs_fname, std::ios::in, buffer_size);
            buff_vec_t ms(ispec.ms_fname, std::ios::out, buffer_size);

            size_type k = 0, h_star = k + 1, h = h_star, h_star_prev = h_star, ms_idx = 0;
            char_type c = t[k];
            node_type v = st.double_rank_nofail_wl(st.root(), c), u = v;
            bool is_maximal = true;

            while (k < t.size()) {
                h = h_star;

                h_star_prev = h_star;
                while (h_star < runs.size()) {
                    c = t[h_star];
                    is_maximal = maxrep.is_maximal(v);
                    u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c, is_maximal);
                    if (!st.is_root(u)) {
                        v = u;
                        h_star += 1;
                    } else
                        break;
                }
                set_next_ms_values1(ms, ms_idx, h, h_star);

                if (h_star < runs.size()) { // remove prefixes of t[k..h*] until you can extend by 'c'
                    //is_maximal = false;
                    bool has_wl = false;
                    u = st.root();
                    do { // remove suffixes of t[k..] until you can extend by 'c'
                        v = st.parent(v);
                        if (!is_maximal) {
                            is_maximal = maxrep.is_intnode_maximal(v);
                        }
                        if (is_maximal) {
                            u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c, is_maximal);
                            has_wl = !st.is_root(u);
                        }
                    } while (!has_wl && !st.is_root(v)); // since !maximal => no wl
                    h_star += 1;
                }
                k = set_next_ms_values2(runs, ms, ms_idx, k);
                v = u;
            }
        }

        static void fill_ms(Stats<cst_t, maxrep_t>& stats, const InputSpec ispec, const cst_t& st,
                wl_method_t1 wl_f_ptr, pseq_method_t pseq_f_ptr,
                const maxrep_t& maxrep, const size_t buffer_size) {

            Query_fwd t{ispec.t_fname, buffer_size};
            buff_vec_t runs(ispec.runs_fname, std::ios::in, buffer_size);
            buff_vec_t ms(ispec.ms_fname, std::ios::out, buffer_size);

            NodeProperty<cst_t, maxrep_t>NP{st, maxrep};
            size_type k = 0, h_star = k + 1, h = h_star, h_star_prev = h_star, ms_idx = 0;
            char_type c = t[k];
            node_type v = st.double_rank_nofail_wl(st.root(), c), u = v;

            while (k < t.size()) {
                h = h_star;

                h_star_prev = h_star;
                while (h_star < ms.size()) {
                    c = t[h_star];
                    u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
                    stats.ms_wl_calls[NP.ms_node_label(v, c)] += 1;
                    if (!st.is_root(u)) {
                        v = u;
                        h_star += 1;
                    } else
                        break;
                }
                stats.ms_wlcalls_seq[(size_type) (h_star - h_star_prev)] += 1; // record
                set_next_ms_values1(ms, ms_idx, h, h_star);

                if (h_star < t.size()) { // remove prefixes of t[k..h*] until you can extend by 'c'
                    u = v;
                    v = pseq_f_ptr(st, wl_f_ptr, v, c);
                    stats.register_ms_pseq(st, NP, u, v, c);
                    h_star += 1;
                }
                k = set_next_ms_values2(runs, ms, ms_idx, k);
                v = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
                stats.ms_wl_calls[NP.ms_node_label(v, c)] += 1;
            }
        }

        static void show_MS(const InputSpec ispec, std::ostream& out) {
            buff_vec_t ms(ispec.ms_fname, std::ios::in);
            //sdsl::bit_vector ms;
            //sdsl::load_from_file(ms, ispec.ms_fname);

            size_type k = 0, max_k = ms.size() / 2;
            for (size_type i = 0; i < ms.size(); i++) {
                if (ms[i] == 1) {
                    out << i - (2 * k) << (++k < max_k ? " " : "");
                }
            }
        }

        static double avg_matching_statistics(const InputSpec ispec) {
            buff_vec_t ms(ispec.ms_fname, std::ios::in);
            //sdsl::bit_vector ms;
            //sdsl::load_from_file(ms, ispec.ms_fname);
            double ans = 0.0;

            size_type k = 0;
            for (size_type i = 0; i < ms.size(); i++) {
                if (ms[i] == 1) {
                    ans += i - (2 * k);
                    k += 1;
                }
            }
            return ans / k;
        }

        static pair<size_type, size_type> ms_composition(const InputSpec ispec) {
            buff_vec_t ms(ispec.ms_fname, std::ios::in);
            //sdsl::bit_vector ms;
            //sdsl::load_from_file(ms, ispec.ms_fname);
            size_type ones = 0;
            for (size_type i = 0; i < ms.size(); i++) {
                if (ms[i] == 1)
                    ones += 1;
            }
            return std::make_pair<size_type, size_type>(ms.size() - ones, ones);
        }

        static size_type size(const InputSpec ispec) {
            buff_vec_t ms(ispec.ms_fname, std::ios::in);
            //sdsl::bit_vector ms;
            //sdsl::load_from_file(ms, ispec.ms_fname);
            return (size_type) ms.size();
        }
    };
}


#endif /* MS_VECTOR_HPP */
