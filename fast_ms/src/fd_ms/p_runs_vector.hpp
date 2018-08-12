/*
 * Represents the runs vector. Implements multithreaded methods. 
 */


#ifndef P_RUNS_VECTOR_HPP
#define P_RUNS_VECTOR_HPP


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

        class p_runs_state {
        public:
            size_type ff_index, lf_index;
            node_type lf_node;

            p_runs_state(const size_type ff_idx = 0, const size_type lf_idx = 0, const node_type lfn = node_type()) :
            ff_index{ff_idx}, lf_index{lf_idx}, lf_node{lfn} {}

            p_runs_state(const p_runs_state& other) :
            ff_index{other.ff_index}, lf_index{other.lf_index}, lf_node{other.lf_node} {}

            p_runs_state& operator=(const p_runs_state& other) {
                ff_index = other.ff_index;
                lf_index = other.lf_index;
                lf_node = other.lf_node;
                return *this;
            }

            string buff_fname(const string base_fname, const Slices<size_type>& slices) const {
                return (base_fname + "." +
                        std::to_string(slices.slice_idx((int) ff_index)) + "_" +
                        std::to_string(slices.slice_idx((int) lf_index)));
            }

            string repr() const {
                return ("[[" + std::to_string(ff_index) + "," + 
                        std::to_string(lf_index) + "), node(" + 
                        std::to_string(lf_node.i) + "," + 
                        std::to_string(lf_node.j) + ")]");
            }

        };
        // strategies for wl operations
        typedef node_type(cst_t::*wl_method_t1) (const node_type& v, const char_type c) const;
        typedef node_type(cst_t::*wl_method_t2) (const node_type& v, const char_type c, const bool is_max) const;

        // strategies for sequences of parent operations
        typedef node_type(*pseq_method_t) (const cst_t& st, wl_method_t1 wl_f_ptr, const node_type& v, const char_type c);

        typedef pair<size_type, size_type> pair_t;
        typedef tuple<size_type, size_type, node_type> alg_state_t;

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
                const vector<p_runs_state> states,
                const size_type buffer_size){

            buff_vec_t out_runs(ispec.runs_fname, std::ios::out, (uint64_t) buffer_size);
            size_type ms_idx = 0;

            for(int slice_idx = 0; slice_idx < slices.nslices; slice_idx++){
                buff_vec_t in_runs(buff_fname(ispec.runs_fname, slice_idx), std::ios::in, (uint64_t) buffer_size);
                size_type in_runs_size = in_runs.size();
                (cerr << " ** adding " << slices.repr(slice_idx) << " from " << 
                        buff_fname(ispec.runs_fname, slice_idx) <<
                        " (" << in_runs_size << " values) ... " << endl);
                // notice that values in in_runs.i are already offset
                for(size_type i = slices[slice_idx].first; i < slices[slice_idx].second; i++)
                    out_runs[i] = in_runs[i];
            }
            // border corrections
            for(int state_idx = 0; state_idx < states.size(); state_idx++){
                p_runs_state state = states[state_idx];

                buff_vec_t in_runs(state.buff_fname(ispec.runs_fname, slices), std::ios::in, (uint64_t) buffer_size);
                size_type in_runs_size = in_runs.size();
                (cerr << " ** adding " << state.repr() << " from " << 
                        state.buff_fname(ispec.runs_fname, slices) <<
                        " (" << in_runs_size << " values) ... " << endl);
                // notice that values in in_runs.i are already offset
                for(size_type i = state.ff_index - 1; i < state.lf_index; i++)
                    out_runs[i] = in_runs[i];
            }
            
        }

        static int fill_inter_slice(const InputSpec ispec, const cst_t& st,
                wl_method_t1 wl_f_ptr, pseq_method_t pseq_f_ptr, 
                const Slices<size_type>& slices, const p_runs_state& state,
                const size_type buffer_size){

            Query_rev t{ispec.t_fname, buffer_size};
            size_type from = state.ff_index - 1, to = state.lf_index, k = to;
            size_type slice_idx = slices.slice_idx(k);
            node_type v = state.lf_node, u = v;
            char_type c = t[k - 1];

            buff_vec_t runs_out(state.buff_fname(ispec.runs_fname, slices), std::ios::out, (uint64_t) buffer_size);
            buff_vec_t runs_in(buff_fname(ispec.runs_fname, slice_idx), std::ios::in, (uint64_t) buffer_size);
            while (--k > from) {
                if(k < slices[slice_idx].first){
                    assert(slice_idx > 0);
                    slice_idx -= 1;
                    runs_in  = buff_vec_t(buff_fname(ispec.runs_fname, slice_idx), std::ios::in, (uint64_t) buffer_size);
                }
                assert(k > from && k < to);
                c = t[k - 1];

                u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
                if (st.is_root(u)) {
                    runs_out[k] = 0;
                    // remove suffixes of t[k..] until you can extend by 'c'
                    v = pseq_f_ptr(st, wl_f_ptr, v, c);
                    u = CALL_MEMBER_FN(st, wl_f_ptr)(v, c);
                } else {
                    runs_out[k] = 1;
                }
                v = CALL_MEMBER_FN(st, wl_f_ptr)(v, c); // update v
            }
            return 0;
        }

        static p_runs_state fill_slice(const InputSpec ispec, const cst_t& st,
                wl_method_t1 wl_f_ptr, pseq_method_t pseq_f_ptr, node_type v,
                const Slices<size_type>& slices, const size_type thread_id, 
                const size_type buffer_size){
            const pair_t slice = slices[thread_id];
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
                    
                    // TODO: last_fail_node and v might be calling wl() unnecessarily

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
            return p_runs_state(first_fail, last_fail, last_fail_node);
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

#endif /* P_RUNS_VECTOR_HPP */


