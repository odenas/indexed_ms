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

        typedef pair<size_type, size_type> pair_t;

        /*
         * The state of the method that fills a slice (during runs construction).
         * Composed of 3 elements:
         *
         * (1) first fail index
         * (2) last fail index
         * (3) the node in which the last fail index hapapened
         * 
         * Slices are filled backwards, hence normally (1) > (2). However, during
         * correction, new states are constructed (see reduce() method) from two
         * different states and the condition is flipped. E.g.,
         *                1 1     2 2     2
         *  0     6 7     3 4     0 1     8
         * |.......|.......|.......|.......|
         *    l f    l        l  f   l  f
         *           f
         * during construction will result in [(4, 2, v4), (8, 8, v8), (19, 16, v19), (25, 22, v22)]
         * but during correction [(4, 16, v16), (19, 22, v22)] (second slice skipped)
         */
        class p_runs_state {
        public:
            size_type ff_index, lf_index;
            node_type lf_node;

            p_runs_state(const size_type ff_idx = 0, const size_type lf_idx = 0, 
                         const node_type lfn = node_type()) :
            ff_index{ff_idx}, lf_index{lf_idx}, lf_node{lfn} 
            {}

            p_runs_state(const p_runs_state& other) :
            ff_index{other.ff_index}, lf_index{other.lf_index}, lf_node{other.lf_node} {}

            p_runs_state& operator=(const p_runs_state& other) {
                ff_index = other.ff_index;
                lf_index = other.lf_index;
                lf_node = other.lf_node;
                return *this;
            }

            string buff_fname(const string base_fname, const Slices<size_type>& slices) const {
                size_type idx1, idx2;
                try{
                    idx1 = slices.slice_idx(ff_index);
                } catch (string s){
                    throw string{"Failed to compute first slice idx for runs buffer filename with array idx: " +
                                 std::to_string(ff_index) + ".\nOriginal message is: " + s};
                }

                try{
                    idx2 = slices.slice_idx((lf_index > 0 ? lf_index - 1 : 0));
                } catch (string s){
                    throw string{"failed to compute second slice idx for runs buffer filename with array idx: " +
                                 std::to_string((lf_index > 0 ? lf_index - 1 : 0)) +
                                 ".\nOriginal message is: " + s};
                }
                return (base_fname + "." + std::to_string(idx1) + "_" + std::to_string(idx2));
            }

            string repr() const {
                return ("[[" + std::to_string(lf_index) + "," + 
                        std::to_string(ff_index) + "), node(" + 
                        std::to_string(lf_node.i) + ", " + 
                        std::to_string(lf_node.j) + ", " +
                        std::to_string(lf_node.ipos) + ", " +
                        std::to_string(lf_node.cipos) + ", " +
                        std::to_string(lf_node.jp1pos) + ")]");
            }

            static p_runs_state slice_covering(const pair_t slice, const node_type v){
                return p_runs_state(slice.second + 1, slice.second + 1, v);
            }
        };
        // strategies for wl operations
        typedef node_type(cst_t::*wl_method_t1) (const node_type& v, const char_type c) const;
        typedef node_type(cst_t::*wl_method_t2) (const node_type& v, const char_type c, const bool is_max) const;

        // strategies for sequences of parent operations
        typedef node_type(*pseq_method_t) (const cst_t& st, wl_method_t1 wl_f_ptr, const node_type& v, const char_type c);

        size_type m_buffer_size;
        Slices<size_type> m_slices;
        wl_method_t1 m_wl_f_ptr;
        pseq_method_t m_pseq_f_ptr;

        p_runs_vector(const size_type buffer_size, const Slices<size_type>& slices, 
                wl_method_t1 wl_f_ptr, pseq_method_t pseq_f_ptr) {
            m_slices = slices;
            m_buffer_size = buffer_size;
            m_wl_f_ptr = wl_f_ptr;
            m_pseq_f_ptr = pseq_f_ptr;
        }
        
        p_runs_vector() = default;
        
        Slices<size_type> slices() const { 
            Slices<size_type> tmp = m_slices;
            return tmp;
        }

        /*
         * Given ending states of the fill_slice method, generate new
         * new states on which to start the fill_inter_slice method.
         *
         * Input states, indicate the first & last position of a failed
         * wl() within a slice(first <= last). 
         *
         * Output states span slices they are of the form
         *  - ff_index=prev_state.ff_index (prev_state is some previous state)
         *  - lf_index=state.lf_index
         *  - lf_node=state.lf_node
         * hence ff_index < lf_index. 
         *
         * If v has size 1, returns an empty vector.
         */
        vector<p_runs_state> reduce(const vector<p_runs_state> v) {
            vector<p_runs_state> u;
            u.reserve(v.size());

            int n = v.size();
            int i = n - 1, j = i - 1;
            while (j >= 0) {
                assert(i == j + 1);

                p_runs_state state = v[i], prev_state = v[j];
                pair_t prev_slice = m_slices[j];
                while(prev_state.ff_index == prev_slice.first + 1 && prev_state.lf_index == prev_slice.first + 1 && j > 0) { // j-th slice had a full match
                    j -= 1;
                }
                assert(j >= 0); // because of the for-loop condition
                prev_state = v[j];
                u.push_back(p_runs_state(prev_state.ff_index, state.lf_index, state.lf_node));
                i = j;
                j = i - 1;
            }
            return u;
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
        
        void merge(const InputSpec ispec, const vector<p_runs_state> states){
            buff_vec_t out_runs(ispec.runs_fname, std::ios::out, (uint64_t) m_buffer_size);

            for(int slice_idx = 0; slice_idx < m_slices.nslices; slice_idx++){
                buff_vec_t in_runs(buff_fname(ispec.runs_fname, slice_idx), std::ios::in, (uint64_t) m_buffer_size);
                //(cerr << " ** adding " << m_slices.repr(slice_idx) << " from " << 
                //        buff_fname(ispec.runs_fname, slice_idx) << endl);
                // notice that values in in_runs.i are already offset
                for(size_type i = m_slices[slice_idx].first; i < m_slices[slice_idx].second; i++){
                    //cerr << "*** out_runs[" << i << "] = " << in_runs[i] << endl;
                    out_runs[i] = in_runs[i];
                }
            }
            // border corrections
            for(int state_idx = 0; state_idx < states.size(); state_idx++){
                p_runs_state state = states[state_idx];

                buff_vec_t in_runs(state.buff_fname(ispec.runs_fname, m_slices), std::ios::in, (uint64_t) m_buffer_size);
                //(cerr << " ** adding " << state.repr() << " from " << 
                //        state.buff_fname(ispec.runs_fname, m_slices) << endl);
                // notice that values in in_runs.i are already offset
                for(size_type i = state.ff_index; i < state.lf_index; i++){
                    //cerr << "*** out_runs[" << i << "] <- " << in_runs[i] << endl;
                    out_runs[i] = in_runs[i];
                }
            }
        }

        /*
        * Generate interleaved run_states from the given sequence of run_states,
        * In the process, remove run_states that spann a full slice.
        */
        size_type fill_inter_slice(const InputSpec ispec, const cst_t& st, 
                const p_runs_state& state){
            
            if(state.ff_index >= state.lf_index)
                throw string("ff_index(" + to_string(state.ff_index) +") < " +
                             "lf_index(" + to_string(state.lf_index) +")");
            if(state.ff_index == 0)
                throw string("expecting ff_index > 0, found: " +
                             to_string(state.ff_index));
            assert(state.ff_index < state.lf_index);

            Query_rev t{ispec.t_fname, m_buffer_size};
            assert(state.ff_index > 0);
            assert(state.lf_index > 0);
            assert(state.ff_index <= state.lf_index);

            size_type from = state.ff_index - 1, k = state.lf_index;
            node_type v = state.lf_node;

            size_type slice_idx = 0;
            try{
                m_slices.slice_idx(k - 1); 
            } catch (string s) {
                throw string{"failed to get slice_idx with message: \n" +
                    s + "\n state is: " + state.repr()};
            }

            buff_vec_t runs_out(state.buff_fname(ispec.runs_fname, m_slices), std::ios::out, (uint64_t) m_buffer_size);
            buff_vec_t runs_in(buff_fname(ispec.runs_fname, slice_idx), std::ios::in, (uint64_t) m_buffer_size);
            while (--k > from) {
                if(k < m_slices[slice_idx].first){
                    assert(slice_idx > 0);
                    slice_idx -= 1;
                    runs_in.close(false);
                    runs_in  = buff_vec_t(buff_fname(ispec.runs_fname, slice_idx), std::ios::in, (uint64_t) m_buffer_size);
                }
                assert(k >= from && k < state.lf_index);
                char_type c = t[k - 1];

                node_type u = CALL_MEMBER_FN(st, m_wl_f_ptr)(v, c);
                if (st.is_root(u)) {
                    runs_out[k] = 0;
                    // remove suffixes of t[k..] until you can extend by 'c'
                    v = m_pseq_f_ptr(st, m_wl_f_ptr, v, c);
                    u = CALL_MEMBER_FN(st, m_wl_f_ptr)(v, c);
                } else {
                    runs_out[k] = 1;
                }
                v = CALL_MEMBER_FN(st, m_wl_f_ptr)(v, c); // update v
            }
            return (size_type) (state.lf_index - state.ff_index);
        }

        p_runs_state fill_slice(const InputSpec ispec, const cst_t& st,
                const size_type thread_id){

            const pair_t slice = m_slices[thread_id];
            const string runs_fname = buff_fname(ispec.runs_fname, thread_id);
            Query_rev t{ispec.t_fname, m_buffer_size};

            if(slice.first >= slice.second)
                throw string{"empty slice: [" + to_string(slice.first) + ", " +
                    to_string(slice.second) + ")"};

            if(slice.second > t.size())
                throw string{"slice endpoint " + to_string(slice.second) + 
                    " > string length " + to_string(t.size())};

            buff_vec_t runs(runs_fname, std::ios::out, (uint64_t) m_buffer_size);

            size_type first_fail = 0, last_fail = 0, from = slice.first, to = slice.second;
            assert(from <= to);

            size_type k = to;
            node_type v = st.double_rank_nofail_wl(st.root(), t[k - 1]);
            node_type u;
            node_type last_fail_node = st.root();
            bool idx_set = false;

            while (--k > from) {
                assert(k > from && k < to);
                char_type c = t[k - 1]; // prev char
                u = CALL_MEMBER_FN(st, m_wl_f_ptr)(v, c);

                if (st.is_root(u)) {
                    if (!idx_set) { // first failing wl()
                        first_fail = k;
                        idx_set = true;
                    }
                    runs[k] = 0;

                    // remove suffixes of t[k..] until you can extend by 'c'
                    v = m_pseq_f_ptr(st, m_wl_f_ptr, v, c);
                    u = CALL_MEMBER_FN(st, m_wl_f_ptr)(v, c);
                    assert(!st.is_root(u));

                    // idx of last 0 in runs - 1 (within this block) and corresponding wl(node)
                    last_fail_node = u; 
                    last_fail = k;
                } else {
                    runs[k] = 1;
                }
                v = u; // update v
            }
            assert(k == from);
            if (!idx_set) {
                /*
                 * (cerr << "**** full match: " << endl
                 *       << "\tfrom = " << from << ", k = " << k << endl
                 *       << "\tu: <" << u.i << ", " << u.j << ">" << endl
                 *       << "\tv: <" << v.i << ", " << v.j << ">" << endl
                 *       << endl);
                 */ 
                first_fail = k + 1;
                last_fail = k + 1;
                last_fail_node = u;
            }
            p_runs_state res(first_fail, last_fail, last_fail_node);
            assert(res.lf_index > 0); // because last_fail > from
            assert(!st.is_root(res.lf_node));
            return res;
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


