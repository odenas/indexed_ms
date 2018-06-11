/*
 * Represents the runs vector. Implements multithreaded methods. 
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
#include "slices.hpp"
#include "query.hpp"
#include "maxrep_vector.hpp"
#include "stats.hpp"

#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))

namespace fdms {

    template<typename cst_t>
    class p_runs_vector {
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

        typedef pair<size_type, size_type> pair_t;
        typedef tuple<size_type, size_type, node_type> alg_state_t;

        size_t nthreads, thread_id;
        pair_t slice;

        p_runs_vector(){
            slice = std::make_pair(0, 0);
            nthreads = 0;
            thread_id = 0;
        }

        p_runs_vector(const size_t nthr, const size_t thr_id, const pair_t slice) :
            nthreads{nthr}, thread_id{thr_id}, slice{slice} {}


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
        
        static void merge(const InputSpec ispec, const Slices<size_type>& slices,
                const size_type buffer_size){

            buff_vec_t in_runs;
            buff_vec_t out_runs(ispec.runs_fname, std::ios::out, (uint64_t) buffer_size);
            for(int slice_idx = 0; slice_idx < slices.nslices; slice_idx++){
                (cerr << " ** adding " << slices.repr(slice_idx) << " from " << 
                        buff_fname(ispec.runs_fname, slice_idx) << " ... " 
                        << endl);
                in_runs = buff_vec_t(buff_fname(ispec.runs_fname, slice_idx), std::ios::in, (uint64_t) buffer_size);
                for(size_type i = slices[slice_idx].first; i < slices[slice_idx].second; i++)
                    out_runs[i] = in_runs[i];
            }
        }

        static int fill_inter_slice(const InputSpec ispec, const cst_t& st,
                wl_method_t1 wl_f_ptr, pseq_method_t pseq_f_ptr, 
                node_type v, size_type slice_idx, const Slices<size_type>& slices,
                const size_type buffer_size){

            Query_rev t{ispec.t_fname, buffer_size};
            pair_t slice = slices[slice_idx];
            size_type from = slice.first, to = slice.second;
            node_type u = v;
            size_type k = to;
            char_type c = t[k - 1];

            buff_vec_t runs(buff_fname(ispec.runs_fname, slice_idx), std::ios::out, (uint64_t) buffer_size);
            while (--k > from) {
                if(k < slices[slice_idx].first){
                    slice_idx -= 1;
                    runs  = buff_vec_t(buff_fname(ispec.runs_fname, slice_idx), std::ios::out, (uint64_t) buffer_size);
                }
                assert(k > from && k < to);
                c = t[k - 1];

                u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
                if (st.is_root(u)) {
                    runs[k] = 0;
                    // remove suffixes of t[k..] until you can extend by 'c'
                    v = pseq_f_ptr(st, wl_f_ptr, v, c);
                    u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
                } else {
                    runs[k] = 1;
                }
                v = CALL_MEMBER_FN(st, wl_f_ptr)(v, c); // update v
            }
            return 0;
        }

        alg_state_t fill_slice(const InputSpec ispec, const cst_t& st,
                wl_method_t1 wl_f_ptr, pseq_method_t pseq_f_ptr, node_type v,
                const size_type buffer_size){
            const string runs_fname = buff_fname(ispec.runs_fname, thread_id);

            Query_rev t{ispec.t_fname, buffer_size};
            buff_vec_t runs(runs_fname, std::ios::out, (uint64_t) buffer_size);
            //cerr << " ** dumping (buffersize: " << runs.buffersize() << ") to " << runs_fname << " ... ";

            size_type first_fail = 0, last_fail = 0, from = slice.first, to = slice.second;
            node_type last_fail_node = v, u = v;

            size_type k = to;
            char_type c = t[k - 1];
            bool idx_set = false;

            while (--k > from) {
                assert(k > from && k < to);
                c = t[k - 1];

                u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
                if (st.is_root(u)) {
                    if (!idx_set) { // first failing wl()
                        first_fail = k;
                        idx_set = true;
                    }
                    runs[k] = 0;

                    // remove suffixes of t[k..] until you can extend by 'c'
                    v = pseq_f_ptr(st, wl_f_ptr, v, c);
                    u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);

                    // idx of last 0 in runs - 1 (within this block) and corresponding wl(node)
                    last_fail_node = CALL_MEMBER_FN(st, wl_f_ptr)(v, c); // given, parent_sequence() above, this has a wl()
                    last_fail = k;
                } else {
                    runs[k] = 1;
                }
                v = CALL_MEMBER_FN(st, wl_f_ptr)(v, c); // update v
            }
            if (!idx_set) {
                first_fail = last_fail = from + 1;
                last_fail_node = v;
            }
            return make_tuple(first_fail, last_fail, last_fail_node);
        }

        // deprecated
        void dump(const InputSpec ispec, const cst_t& st,
                wl_method_t1 wl_f_ptr, pseq_method_t pseq_f_ptr, const size_t buffer_size) {


            Query_rev t{ispec.t_fname, buffer_size};
            //buff_vec_t runs(ispec.runs_fname, std::ios::out, (uint64_t) buffer_size);
            sdsl::bit_vector runs(t.size());
            //cerr << " ** dumping (buffersize: " << runs.buffersize() << ") to " << ispec.runs_fname << " ... ";

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
            cerr << "DONE" << endl;
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
    };
}

#endif /* RUNS_VECTOR_HPP */


