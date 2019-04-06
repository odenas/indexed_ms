/*
 * Represents the ms vector. Implements multithreaded methods.
 */


#ifndef P_MS_VECTOR_HPP
#define P_MS_VECTOR_HPP


#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <sdsl/vectors.hpp>

#include "input_spec.hpp"
#include "stats.hpp"
#include "slices.hpp"
#include "query.hpp"
#include "maxrep_vector.hpp"
#include "stats.hpp"

#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))

namespace fdms {

    class ms_compression {
    public:
        enum compression_types {
            none, rrr, hyb, delta, succint, nibble, rle
        };

    private:
        static std::map<compression_types, string> c2s() {
            std::map<compression_types, string> c2s = {
                {compression_types::none, ""},
                {compression_types::rrr, ".rrr"},
                {compression_types::hyb, ".hyb"},
                {compression_types::delta, ".delta"},
                {compression_types::succint, ".succint"},
                {compression_types::nibble, ".nibble"},
                {compression_types::rle, ".rle"},
            };
            return c2s;
        }

    public:
        ms_compression() = default;

        static string to_str(const compression_types ct){
            return c2s()[ct];
        }

        static compression_types parse_compression(const string c_str){
            for(auto item: c2s()){
                if(item.second == ("." + c_str))
                    return item.first;
            }
            if(c_str == "none")
                return compression_types::none;
            throw string{"bad compression string: " + c_str};
        }
    };

    template<typename cst_t>
    class p_ms_vector {
    public:
        typedef typename cst_t::size_type size_type;
        typedef typename cst_t::char_type char_type;
        typedef typename cst_t::node_type node_type;

        typedef Maxrep<cst_t, sdsl::bit_vector> maxrep_t;
        typedef sdsl::int_vector_buffer<1> buff_vec_t;

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

        static size_type select11(buff_vec_t& v){
            for(size_type i=0; i < v.size(); i++)
                if(v[i] == 1)
                    return i;
            return -1;
        }


    public:
        // strategies for rank operations
        typedef node_type(cst_t::*wl_method_t1) (const node_type& v, const char_type c) const;
        typedef node_type(cst_t::*wl_method_t2) (const node_type& v, const char_type c, const bool is_max) const;

        // strategies for sequences of parent operations
        typedef node_type(*pseq_method_t) (const cst_t& st, wl_method_t1 wl_f_ptr, const node_type& v, const char_type c);

        typedef pair<size_type, size_type> pair_t;
        typedef tuple<size_type, size_type, node_type> alg_state_t;

        size_t m_nthreads, m_buffer_size;
        wl_method_t1 m_wl_f_ptr;
        pseq_method_t m_pseq_f_ptr;

        p_ms_vector() = default;

        p_ms_vector(const size_t nthr, const size_type buffer_size,
                wl_method_t1 wl_f_ptr, pseq_method_t pseq_f_ptr){
            m_nthreads = nthr;
            m_buffer_size = buffer_size;
            m_wl_f_ptr = wl_f_ptr;
            m_pseq_f_ptr = pseq_f_ptr;
        }

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

        static string buff_fname(const string base_fname, size_type thr_id) {
            return base_fname + "." + std::to_string(thr_id);
        }

        size_type fill_slice(const InputSpec& ispec, cst_t& st,
                const pair_t slice, size_type thread_id){
            
            pseq_method_t pseq_f_ptr = m_pseq_f_ptr;
            wl_method_t1 wl_f_ptr = m_wl_f_ptr;
            Query_fwd t{ispec.t_fname, m_buffer_size};

            buff_vec_t runs(ispec.runs_fname, std::ios::in, m_buffer_size);
            buff_vec_t ms(buff_fname(ispec.ms_fname, thread_id), std::ios::out, m_buffer_size);

            size_type from = slice.first, to = slice.second;
            size_type k = from, h_star = k + 1, h = h_star, ms_idx = 0, ms_size = t.size();
            char_type c = t[k];
            node_type v = CALL_MEMBER_FN(st, wl_f_ptr)(st.root(), c), u = v;

            while (k < to) {
                h = h_star;

                while (h_star < ms_size) {
                    c = t[h_star];
                    u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
                    if (!st.is_root(u)) {
                        v = u;
                        h_star += 1;
                    } else
                        break;
                }
                set_next_ms_values1(ms, ms_idx, h, h_star);

                if (h_star < ms_size) { // remove prefixes of t[k..h*] until you can extend by 'c'
                    v = pseq_f_ptr(st, wl_f_ptr, v, c);
                    h_star += 1;
                }
                k = set_next_ms_values2(runs, ms, ms_idx, k);
                v = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
            }
            return ms_idx;
        }

        /*
         * Given slices, find the corresponding MS buffer vector files
         * and merge them into a single buffer vector file. This is not
         * just a concatenation. There is a correction step the details
         * of which I don't remember. 
         */
        void merge(const InputSpec ispec, const Slices<size_type>& slices){

            buff_vec_t out_ms(ispec.ms_fname, std::ios::out, (uint64_t) m_buffer_size);
            size_type ms_idx = 0;
            size_type correction = 0;

            for(int slice_idx = 0; slice_idx < slices.nslices; slice_idx++){
                buff_vec_t in_ms(buff_fname(ispec.ms_fname, slice_idx), std::ios::in, (uint64_t) m_buffer_size);
                size_type in_ms_size = in_ms.size();
                if(in_ms[in_ms_size - 1] != 1)
                    throw string{"expecting 1 at the last position of " +
                                 buff_fname(ispec.ms_fname, slice_idx) +
                                 " found " + to_string(in_ms[in_ms_size - 1])};
                (cerr << " ** adding " << slices.repr(slice_idx) << " from " <<
                        buff_fname(ispec.ms_fname, slice_idx) << endl);

                size_type i = correction;
                bool corrected = true;
                while(i < in_ms_size){
                    if(corrected){
                        out_ms[ms_idx++] = in_ms[i++];
                    } else {
                        i = select11(in_ms);
                        cerr << "i: " << i << ", correction: " << correction << endl;
                        assert(i >= correction);
                        for(size_type j = 0; j < i - correction; j++)
                            out_ms[ms_idx++] = in_ms[j];
                        corrected = true;
                    }
                }
                //cerr << "in_ms_size: " << in_ms_size << ", 2 * slice_len: " << 2 * slices.slice_length(slice_idx) << endl;
                assert(in_ms_size >= (2 * slices.slice_length(slice_idx)));
                correction = in_ms_size - (2 * slices.slice_length(slice_idx));
                //correction = 0;
            }
        }

    };
}

#endif /* P_MS_VECTOR_HPP */



