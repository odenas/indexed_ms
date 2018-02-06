/* sdsl - succinct data structures library
    Copyright (C) 2013 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file wt_pc.hpp
    \brief wt_pc.hpp contains a class for the wavelet tree of byte sequences.
           The wavelet tree shape is parametrized by a prefix code.
    \author Simon Gog, Timo Beller
*/
#ifndef INCLUDED_SDSL_WT_PC
#define INCLUDED_SDSL_WT_PC

#include "bit_vectors.hpp"
#include "rank_support.hpp"
#include "select_support.hpp"
#include "wt_helper.hpp"
#include <vector>
#include <utility>
#include <tuple>

//! Namespace for the succinct data structure library.
namespace sdsl
{

//! A prefix code-shaped wavelet.
/*!
 * \tparam t_shape       Shape of the tree ().
 * \tparam t_bitvector   Underlying bitvector structure.
 * \tparam t_rank        Rank support for pattern `1` on the bitvector.
 * \tparam t_select      Select support for pattern `1` on the bitvector.
 * \tparam t_select_zero Select support for pattern `0` on the bitvector.
 * \tparam t_tree_strat  Tree strategy determines alphabet and the tree
 *                       class used to navigate the WT.
 *
 *  @ingroup wt
 */
template<class t_shape,
         class t_bitvector   = bit_vector,
         class t_rank        = typename t_bitvector::rank_1_type,
         class t_select      = typename t_bitvector::select_1_type,
         class t_select_zero = typename t_bitvector::select_0_type,
         class t_tree_strat  = byte_tree<>
         >
class wt_pc
{
    public:
        typedef typename
        t_tree_strat::template type<wt_pc>            tree_strat_type;
        typedef int_vector<>::size_type               size_type;
        typedef typename
        tree_strat_type::value_type                   value_type;
        typedef typename t_bitvector::difference_type difference_type;
        typedef random_access_const_iterator<wt_pc>   const_iterator;
        typedef const_iterator                        iterator;
        typedef t_bitvector                           bit_vector_type;
        typedef t_rank                                rank_1_type;
        typedef t_select                              select_1_type;
        typedef t_select_zero                         select_0_type;
        typedef wt_tag                                index_category;
        typedef typename
        tree_strat_type::alphabet_category            alphabet_category;
        typedef typename
        t_shape::template type<wt_pc>                 shape_type;
        enum { lex_ordered=shape_type::lex_ordered };
        using node_type = typename tree_strat_type::node_type;

    private:

#ifdef WT_PC_CACHE
        mutable value_type m_last_access_answer;
        mutable size_type  m_last_access_i;
        mutable size_type  m_last_access_rl;
#endif

        size_type        m_size  = 0;    // original text size
        size_type        m_sigma = 0;    // alphabet size
        bit_vector_type  m_bv;           // bit vector to store the wavelet tree
        rank_1_type      m_bv_rank;      // rank support for the wavelet tree bit vector
        select_1_type    m_bv_select1;   // select support for the wavelet tree bit vector
        select_0_type    m_bv_select0;
        tree_strat_type  m_tree;

        void copy(const wt_pc& wt)
        {
            m_size            = wt.m_size;
            m_sigma           = wt.m_sigma;
            m_bv              = wt.m_bv;
            m_bv_rank         = wt.m_bv_rank;
            m_bv_rank.set_vector(&m_bv);
            m_bv_select1    = wt.m_bv_select1;
            m_bv_select1.set_vector(&m_bv);
            m_bv_select0    = wt.m_bv_select0;
            m_bv_select0.set_vector(&m_bv);
            m_tree          = wt.m_tree;
        }

        // insert a character into the wavelet tree, see construct method
        void insert_char(value_type old_chr, std::vector<uint64_t>& bv_node_pos,
                         size_type times, bit_vector& bv)
        {
            uint64_t p = m_tree.bit_path(old_chr);
            uint32_t path_len = p>>56;
            node_type v = m_tree.root();
            for (uint32_t l=0; l<path_len; ++l, p >>= 1) {
                if (p&1) {
                    bv.set_int(bv_node_pos[v], 0xFFFFFFFFFFFFFFFFULL,times);
                }
                bv_node_pos[v] += times;
                v = m_tree.child(v, p&1);
            }
        }



        // calculates the tree shape returns the size of the WT bit vector
        size_type construct_tree_shape(const std::vector<size_type>& C)
        {
            // vector  for node of the tree
            std::vector<pc_node> temp_nodes; //(2*m_sigma-1);
            shape_type::construct_tree(C, temp_nodes);
            // Convert code tree into BFS order in memory and
            // calculate bv_pos values
            size_type bv_size = 0;
            tree_strat_type temp_tree(temp_nodes, bv_size, this);
            m_tree.swap(temp_tree);
            return bv_size;
        }

        void construct_init_rank_select()
        {
            util::init_support(m_bv_rank, &m_bv);
            util::init_support(m_bv_select0, &m_bv);
            util::init_support(m_bv_select1, &m_bv);
        }

        // recursive internal version of the method interval_symbols
        void
        _interval_symbols(size_type i, size_type j, size_type& k,
                          std::vector<value_type>& cs,
                          std::vector<size_type>& rank_c_i,
                          std::vector<size_type>& rank_c_j, node_type v) const
        {
            // invariant: j>i
            size_type i_new = (m_bv_rank(m_tree.bv_pos(v) + i)
                               - m_tree.bv_pos_rank(v));
            size_type j_new = (m_bv_rank(m_tree.bv_pos(v) + j)
                               - m_tree.bv_pos_rank(v));
            // goto left child
            i -= i_new; j -= j_new;
            if (i != j) {
                node_type v_new = m_tree.child(v, 0);
                if (!m_tree.is_leaf(v_new)) {
                    _interval_symbols(i, j, k, cs, rank_c_i, rank_c_j, v_new);
                } else {
                    rank_c_i[k] = i;
                    rank_c_j[k] = j;
                    cs[k++] = m_tree.bv_pos_rank(v_new);
                }
            }
            // goto right child
            if (i_new!=j_new) {
                node_type v_new = m_tree.child(v, 1);
                if (!m_tree.is_leaf(v_new)) {
                    _interval_symbols(i_new, j_new, k, cs, rank_c_i, rank_c_j,
                                      v_new);
                } else {
                    rank_c_i[k] = i_new;
                    rank_c_j[k] = j_new;
                    cs[k++] = m_tree.bv_pos_rank(v_new);
                }
            }
        }

    public:

        const size_type&       sigma = m_sigma;
        const bit_vector_type& bv  = m_bv;

        // some accessor methods
        rank_1_type get_bv_rank() const { return m_bv_rank; }
        tree_strat_type get_tree() const { return m_tree; }
        select_1_type get_bv_select1() const { return m_bv_select1; }
        select_0_type get_bv_select0() const { return m_bv_select0; }
        size_type get_sigma() const {return m_sigma; }
        bit_vector_type get_bv() const {return m_bv; }


        // Default constructor
        wt_pc() {};

        //! Construct the wavelet tree from a file_buffer
        /*!
         * \param input_buf    File buffer of the input.
         * \param size         The length of the prefix.
         * \par Time complexity
         *      \f$ \Order{n\log|\Sigma|}\f$, where \f$n=size\f$
         */
        wt_pc(int_vector_buffer<tree_strat_type::int_width>& input_buf,
              size_type size):m_size(size)
        {
            if (0 == m_size)
                return;
            // O(n + |\Sigma|\log|\Sigma|) algorithm for calculating node sizes
            // TODO: C should also depend on the tree_strategy. C is just a mapping
            // from a symbol to its frequency. So a map<uint64_t,uint64_t> could be
            // used for integer alphabets...
            std::vector<size_type> C;
            // 1. Count occurrences of characters
            calculate_character_occurences(input_buf, m_size, C);
            // 2. Calculate effective alphabet size
            calculate_effective_alphabet_size(C, m_sigma);
            // 3. Generate tree shape
            size_type tree_size = construct_tree_shape(C);
            // 4. Generate wavelet tree bit sequence m_bv
            bit_vector temp_bv(tree_size, 0);

            // Initializing starting position of wavelet tree nodes
            std::vector<uint64_t> bv_node_pos(m_tree.size(), 0);
            for (size_type v=0; v < m_tree.size(); ++v) {
                bv_node_pos[v] = m_tree.bv_pos(v);
            }
            if (input_buf.size() < size) {
                throw std::logic_error("Stream size is smaller than size!");
                return;
            }
            value_type old_chr = input_buf[0];
            uint32_t times = 0;
            for (size_type i=0; i < m_size; ++i) {
                value_type chr = input_buf[i];
                if (chr != old_chr) {
                    insert_char(old_chr, bv_node_pos, times, temp_bv);
                    times = 1;
                    old_chr = chr;
                } else { // chr == old_chr
                    ++times;
                    if (times == 64) {
                        insert_char(old_chr, bv_node_pos, times, temp_bv);
                        times = 0;
                    }
                }
            }
            if (times > 0) {
                insert_char(old_chr, bv_node_pos, times, temp_bv);
            }
            m_bv = bit_vector_type(std::move(temp_bv));
            // 5. Initialize rank and select data structures for m_bv
            construct_init_rank_select();
            // 6. Finish inner nodes by precalculating the bv_pos_rank values
            m_tree.init_node_ranks(m_bv_rank);
        }


        //! Copy constructor
        wt_pc(const wt_pc& wt) { copy(wt); }

        wt_pc(wt_pc&& wt)
        {
            *this = std::move(wt);
        }

        //! Assignment operator
        wt_pc& operator=(const wt_pc& wt)
        {
            if (this != &wt) {
                copy(wt);
            }
            return *this;
        }

        //! Assignment operator
        wt_pc& operator=(wt_pc&& wt)
        {
            if (this != &wt) {
                m_size            = wt.m_size;
                m_sigma           = wt.m_sigma;
                m_bv              = std::move(wt.m_bv);
                m_bv_rank         = std::move(wt.m_bv_rank);
                m_bv_rank.set_vector(&m_bv);
                m_bv_select1    = std::move(wt.m_bv_select1);
                m_bv_select1.set_vector(&m_bv);
                m_bv_select0    = std::move(wt.m_bv_select0);
                m_bv_select0.set_vector(&m_bv);
                m_tree          = std::move(wt.m_tree);
            }
            return *this;
        }


        //! Swap operator
        void swap(wt_pc& wt)
        {
            if (this != &wt) {
                std::swap(m_size, wt.m_size);
                std::swap(m_sigma,  wt.m_sigma);
                m_bv.swap(wt.m_bv);
                util::swap_support(m_bv_rank, wt.m_bv_rank,
                                   &m_bv, &(wt.m_bv));

                util::swap_support(m_bv_select1, wt.m_bv_select1,
                                   &m_bv, &(wt.m_bv));
                util::swap_support(m_bv_select0, wt.m_bv_select0,
                                   &m_bv, &(wt.m_bv));
                m_tree.swap(wt.m_tree);
            }
        }

        //! Returns the size of the original vector.
        size_type size()const { return m_size; }

        //! Returns whether the wavelet tree contains no data.
        bool empty()const { return m_size == 0; }

        //! Recovers the i-th symbol of the original vector.
        /*!
         * \param i Index in the original vector.
         * \return The i-th symbol of the original vector.
         * \par Time complexity
         *      \f$ \Order{H_0} \f$ on average, where \f$ H_0 \f$ is the
         *      zero order entropy of the sequence
         *
         * \par Precondition
         *      \f$ i < size() \f$
         */
        value_type operator[](size_type i)const
        {
            assert(i < size());
            // which stores how many of the next symbols are equal
            // with the current char
            node_type v = m_tree.root(); // start at root node
            while (!m_tree.is_leaf(v)) {   // while  not a leaf
                if (m_bv[ m_tree.bv_pos(v) + i]) {  // goto right child
                    i = m_bv_rank(m_tree.bv_pos(v) + i)
                        - m_tree.bv_pos_rank(v);
                    v = m_tree.child(v,1);
                } else { // goto the left child
                    i -= (m_bv_rank(m_tree.bv_pos(v) + i)
                          - m_tree.bv_pos_rank(v));
                    v = m_tree.child(v,0);
                }
            }
            // if v is a leaf bv_pos_rank returns symbol itself
            return m_tree.bv_pos_rank(v);
        };

        //! Calculates how many symbols c are in the prefix [0..i-1].
        /*!
         * \param i Exclusive right bound of the range.
         * \param c Symbol c.
         * \return Number of occurrences of symbol c in the prefix [0..i-1].
         * \par Time complexity
         *      \f$ \Order{H_0} \f$ on average, where \f$ H_0 \f$ is the
         *      zero order entropy of the sequence
         *
         * \par Precondition
         *      \f$ i \leq size() \f$
         */
        size_type rank(size_type i, value_type c)const
        {
            assert(i <= size());
            if (!m_tree.is_valid(m_tree.c_to_leaf(c))) {
                return 0;  // if `c` was not in the text
            }
            if (m_sigma == 1) {
                return i; // if m_sigma == 1 answer is trivial
            }
            uint64_t p = m_tree.bit_path(c);
            uint32_t path_len = (p>>56);
            size_type result = i;
            node_type v = m_tree.root();
            for (uint32_t l=0; l<path_len; ++l, p >>= 1) {
                if (p&1) {
                    result  = (m_bv_rank(m_tree.bv_pos(v)+result) -  m_tree.bv_pos_rank(v));
                } else {
                    result -= (m_bv_rank(m_tree.bv_pos(v)+result) -  m_tree.bv_pos_rank(v));
                }
                if(l == (path_len - 1))
                    return result;
                v = m_tree.child(v, p&1); // goto child
            }
            return result;
        };

        /**
         * this is bein called on end of a non-maximal repeat interval
         */
        size_type rank_and_check(size_type i, value_type c)const
        {
            size_type j;

            assert(i <= size());
            assert(i > 0);
            if (!m_tree.is_valid(m_tree.c_to_leaf(c))) {
                return 0;  // if `c` was not in the text
            }
            if (m_sigma == 1) {
                return i; // if m_sigma == 1 answer is trivial
            }
            // access wt[i-1]
            {
                j = i - 1;
                uint64_t p = m_tree.bit_path(c);
                bool c_in_right_subtree = (p&1);
                node_type v = m_tree.root(); // start at root node
                while (!m_tree.is_leaf(v)) {
                    if (m_bv[m_tree.bv_pos(v) + j]) {  // goto right child
                        if(!c_in_right_subtree)
                            return 0;
                        j = m_bv_rank(m_tree.bv_pos(v) + j) - m_tree.bv_pos_rank(v);
                        v = m_tree.child(v,1);
                    } else { // goto the left child
                        if(c_in_right_subtree)
                            return 0;
                        j -= (m_bv_rank(m_tree.bv_pos(v) + j) - m_tree.bv_pos_rank(v));
                        v = m_tree.child(v,0);
                    }
                    p >>= 1;
                    c_in_right_subtree = (p&1);
                }
                // if v is a leaf bv_pos_rank returns symbol itself
                if (m_tree.bv_pos_rank(v) != c)
                    return 0; // c is not in the interval
                // we have already navigated down to leaf
            }
            return j + 1;
        };


        //! Calculates how many symbols c are in the prefix [0..i-1] and [0..j-1].
        /*!
         *  \par This is equivalent to calling rank(i, c) and rank(j, c). However, this
         *  is more efficient since it goes down the WT only once. Furthermore, it will
         *  return early if the two ranks are similar.
         *
         *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in [0..size()]\f$.
         *  \param j The exclusive index of the prefix range [0..j-1], so \f$i\in [0..size()]\f$.
         *  \param c The symbol to count the occurrences in the prefix.
         *    \returns A std::pair with the number of occurrences of symbol c in the prefix [0..i-1] and [0..j-1].
         *  \par Time complexity
         *        \f$ \Order{\log |\Sigma|} \f$
         */
        std::pair<size_type, size_type>
        double_rank_and_fail(size_type i, size_type j, value_type c) const{
            assert(i <= j);
            assert(j <= size());
            if (!m_tree.is_valid(m_tree.c_to_leaf(c))) {
                return std::make_pair(0, 0); // if `c` was not in the text
            }
            if (m_sigma == 1) {
                return std::make_pair(i, j); // if m_sigma == 1 answer is trivial
            }
            uint64_t p = m_tree.bit_path(c);
            uint32_t path_len = (p>>56);
            size_type result_i = i;
            size_type result_j = j;

            node_type v = m_tree.root();
            uint32_t l = 0;
            for (;result_i; ++l, p >>= 1) {
                std::pair<size_type, size_type> a = m_bv_rank.double_rank(m_tree.bv_pos(v) + result_i,
                                                                          m_tree.bv_pos(v) + result_j);
                size_type b = m_tree.bv_pos_rank(v);

                if(p&1){
                    result_i = (a.first - b);
                    result_j = (a.second - b);
                } else {
                    result_i -= (a.first -  b);
                    result_j -= (a.second -  b);
                }
                // If l==path_len-1 we have reached a leaf, so we can return result. 
                // i and j have the same rank. the Wl call will fail
                if(l == (path_len - 1) or result_i == result_j)
                    return std::make_pair(result_i, result_j); // i and j have the same rank. the Wl call will fail
                v = m_tree.child(v, p&1); // goto child
            }
            for (; l<path_len and result_j; ++l, p >>= 1) {
                if(p&1){
                    result_j = (m_bv_rank(m_tree.bv_pos(v)+result_j) -  m_tree.bv_pos_rank(v));
                } else {
                    result_j -= (m_bv_rank(m_tree.bv_pos(v)+result_j) -  m_tree.bv_pos_rank(v));
                }
                v = m_tree.child(v, p&1); // goto child
            }
            return std::make_pair(result_i, result_j);
        }

        std::pair<size_type, size_type>
        double_rank(size_type i, size_type j, value_type c) const{
            assert(i <= size());
            assert(j <= size());
            if (!m_tree.is_valid(m_tree.c_to_leaf(c))) {
                return std::make_pair(0, 0); // if `c` was not in the text
            }
            if (m_sigma == 1) {
                return std::make_pair(i, j); // if m_sigma == 1 answer is trivial
            }
            uint64_t p = m_tree.bit_path(c);
            uint32_t path_len = (p>>56);
            size_type result_i = i;
            size_type result_j = j;

            node_type v = m_tree.root();
            uint32_t l = 0;
            for (;result_i; ++l, p >>= 1) {
                std::pair<size_type, size_type> a = m_bv_rank.double_rank(m_tree.bv_pos(v) + result_i,
                                                                          m_tree.bv_pos(v) + result_j);
                size_type b = m_tree.bv_pos_rank(v);

                if(p&1){
                    result_i = (a.first - b);
                    result_j = (a.second - b);
                } else {
                    result_i -= (a.first -  b);
                    result_j -= (a.second -  b);
                }
                if(l == path_len-1)
                    return std::make_pair(result_i, result_j);	
                v = m_tree.child(v, p&1); // goto child
            }
            for (; (l < path_len) and result_j; ++l, p >>= 1) {
                if(p&1){
                    result_j = (m_bv_rank(m_tree.bv_pos(v)+result_j) -  m_tree.bv_pos_rank(v));
                } else {
                    result_j -= (m_bv_rank(m_tree.bv_pos(v)+result_j) -  m_tree.bv_pos_rank(v));
                }
                v = m_tree.child(v, p&1); // goto child
            }
            return std::make_pair(result_i, result_j);
        }

        //! Calculates how many times symbol wt[i] occurs in the prefix [0..i-1].
        /*!
         * \param i The index of the symbol.
         * \return  Pair (rank(wt[i],i),wt[i])
         * \par Time complexity
         *      \f$ \Order{H_0} \f$
         *
         * \par Precondition
         *      \f$ i < size() \f$
         */
        std::pair<size_type, value_type>
        inverse_select(size_type i)const
        {
            assert(i < size());
            node_type v = m_tree.root();
            while (!m_tree.is_leaf(v)) {   // while not a leaf
                if (m_bv[m_tree.bv_pos(v) + i]) {   //  goto right child
                    i = (m_bv_rank(m_tree.bv_pos(v) + i)
                         - m_tree.bv_pos_rank(v));
                    v = m_tree.child(v, 1);
                } else { // goto left child
                    i -= (m_bv_rank(m_tree.bv_pos(v) + i)
                          - m_tree.bv_pos_rank(v));
                    v = m_tree.child(v,0);
                }
            }
            // if v is a leaf bv_pos_rank returns symbol itself
            return std::make_pair(i, (value_type)m_tree.bv_pos_rank(v));
        }

        //! Calculates the ith occurrence of the symbol c in the supported vector.
        /*!
         * \param i The ith occurrence.
         * \param c The symbol c.
         * \par Time complexity
         *      \f$ \Order{H_0} \f$ on average, where \f$ H_0 \f$ is the zero order
         *       entropy of the sequence
         *
         * \par Precondition
         *      \f$ 1 \leq i \leq rank(size(), c) \f$
         */
        size_type select(size_type i, value_type c)const
        {
            assert(1 <= i and i <= rank(size(), c));
            node_type v = m_tree.c_to_leaf(c);
            if (!m_tree.is_valid(v)) {   // if c was not in the text
                return m_size;         // -> return a position right to the end
            }
            if (m_sigma == 1) {
                return std::min(i-1,m_size);
            }
            size_type result = i-1;    // otherwise
            uint64_t p = m_tree.bit_path(c);
            uint32_t path_len = (p>>56);
            // path_len > 0, since we have handled m_sigma = 1.
            p <<= (64-path_len);
            for (uint32_t l=0; l<path_len; ++l, p <<= 1) {
                if ((p & 0x8000000000000000ULL)==0) { // node was a left child
                    v  = m_tree.parent(v);
                    result = m_bv_select0(m_tree.bv_pos(v)
                                          - m_tree.bv_pos_rank(v) + result + 1)
                             - m_tree.bv_pos(v);
                } else { // node was a right child
                    v   = m_tree.parent(v);
                    result = m_bv_select1(m_tree.bv_pos_rank(v) + result + 1)
                             - m_tree.bv_pos(v);
                }
            }
            return result;
        };


        /* functions below are to be used as aliases */
        inline size_type bit_rank1(const node_type v, const size_type i) const {
            return m_bv_rank(m_tree.bv_pos(v) + i) - m_tree.bv_pos_rank(v);
        }

        inline size_type bit_rank0(const node_type v, const size_type i) const {
            return i - bit_rank1(v, i);
        }

        inline size_type bit_select1(const node_type v, const size_type i) const {
            return m_bv_select1(m_tree.bv_pos_rank(v) + i) - m_tree.bv_pos(v);
        }

        inline size_type bit_select0(const node_type v, const size_type i) const {
            return m_bv_select0(m_tree.bv_pos(v) - m_tree.bv_pos_rank(v) + i) - m_tree.bv_pos(v);
        }

        inline size_type bit_select_at_dist0(const node_type v, const size_type i, const size_type cnt) const {
            if(i == 0xFFFFFFFFFFFFFFFF) // i was supposed to be -1
                return bit_select0(v, cnt);

            if(cnt == 0)
                return bit_select0(v, bit_rank0(v, i + 1) + cnt);

            uint32_t idx = (uint32_t)(i + m_tree.bv_pos(v));
            uint64_t word = *(m_bv.data() + (idx >> 6));
            size_type word_offset = (idx % 64) + 1;
            word = ~word;
            word = (word >> word_offset) << word_offset; // remove ones up to index i

            if(word_offset == 64 || word == 0 || sdsl::bits::cnt(word) < cnt)
                return bit_select0(v, bit_rank0(v, i + 1) + cnt);
            uint32_t word_ans = sdsl::bits::sel(word, (uint32_t)cnt);
            uint32_t ans = word_ans + ((idx >> 6) << 6) - (uint32_t)m_tree.bv_pos(v);
            assert(ans > 0);
            assert(ans == bit_select0(v, bit_rank0(v, i + 1) + cnt));
            return ans;
        }

        inline size_type bit_select_at_dist1(const node_type v, const size_type i, const size_type cnt) const {
            /*
            load the current word and check if there are cnt 1s in it
             if yes, compute the index and return it
             if not, call select(next(i+1), cnt) on the bit vector
            */

            if(i == 0xFFFFFFFFFFFFFFFF) // i was supposed to be -1
                return bit_select1(v, cnt);

            if(cnt == 0)
                return bit_select1(v, bit_rank1(v, i + 1) + cnt);

            uint32_t idx = (uint32_t) (i + m_tree.bv_pos(v));
            uint64_t word = *(m_bv.data() + (idx >> 6));
            size_type word_offset = (idx % 64) + 1;
            word = (word >> word_offset) << word_offset; // remove ones up to index i
            if(word_offset == 64 || word == 0 || sdsl::bits::cnt(word) < cnt)
                return bit_select1(v, bit_rank1(v, i + 1) + cnt);
            uint32_t word_ans = sdsl::bits::sel(word, (uint32_t)cnt);
            uint32_t ans = word_ans + ((idx >> 6) << 6) - (uint32_t)m_tree.bv_pos(v);
            assert(ans > 0);
            assert(ans == bit_select1(v, bit_rank1(v, i + 1) + cnt));
            return ans;
        }


        //! Calculate the cnt-th occurrence of the symbol c in positions [i..size() - 1]
        /*!
         * \param c The symbol
         * \param i The index
         * \param cnt The cnt-th occurrence.
         *
         * \par Precondition
         *      \f$ 1 \leq i \leq size() \f$
         * \par Precondition
         *      \f$ 1 \leq d \leq rank(size(), c) \f$
         */
        size_type select_at_dist(const value_type c, const size_type i, const size_type cnt) const
        {
            uint64_t p = m_tree.bit_path(c), pt_i = m_tree.bit_path((*this)[i]);
            uint32_t p_len = (p >> 56);
            std::vector<bool> equal_prefix(p_len);
            node_type v = m_tree.root();
            std::vector<size_type> i_vec(p_len); // place the i-values here
            i_vec[0] = i; equal_prefix[0] = true;

            for(uint32_t i = 1; i < p_len; i++, p >>= 1, pt_i >>= 1) {
                if(((p&1) == (pt_i&1)) && equal_prefix[i - 1])
                    equal_prefix[i] = true;
                if (p&1)
                    i_vec[i] = bit_rank1(v, i_vec[i - 1]);
                else
                    i_vec[i] = bit_rank0(v, i_vec[i - 1]);
                //cout << "i" << i << " = r" << (p&1) << "(" << i_vec[i - 1] << ") = " << i_vec[i] << endl;
                v = m_tree.child(v, p&1); // goto child
            }
            /*
            cout << "i_vec: [";
            for(int aa = 0; aa < p_len - 1; aa++)
                cout << i_vec[aa] << ", ";
            cout <<  i_vec[p_len - 1] << "]" << endl;
            */

            // reset path
            p = m_tree.bit_path(c);
            p <<= (64- (p >> 56));

            size_type jk = 0, j_prev = 0;
            if(!equal_prefix[p_len - 1])
                i_vec[p_len - 1] -= 1; // this migh undeflow, but that's fine
            if((p & 0x8000000000000000ULL) == 0)
                j_prev = bit_select_at_dist0(v, i_vec[p_len - 1], cnt);
            else
                j_prev = bit_select_at_dist1(v, i_vec[p_len - 1], cnt);
            //cout << (p != pt_i ? "*" : "") << "j" << p_len << " = sd" << ((p & 0x8000000000000000ULL)!=0) << "(" << i_vec[p_len - 1] << ", " << cnt << ") = " << j_prev << endl;
            v  = m_tree.parent(v);
            p <<= 1;

            for(uint32_t idx = p_len - 1; idx > 0; idx--, p <<= 1){
                if(!equal_prefix[idx - 1])
                    i_vec[idx - 1] -= 1;// this migh undeflow, but that's fine
                if ((p & 0x8000000000000000ULL)==0)
                    jk = bit_select_at_dist0(v, i_vec[idx - 1], j_prev - i_vec[idx]);
                else
                    jk = bit_select_at_dist1(v, i_vec[idx - 1], j_prev - i_vec[idx]);
                //cout << (p != pt_i ? "*" : "") << "j" << idx << " = sd" << ((p & 0x8000000000000000ULL)!=0) << "(" << i_vec[idx - 1] << ", " << cnt << ") = " << jk << endl;
                v = m_tree.parent(v);
                j_prev = jk;
            }
            assert(j_prev == select(rank(i + 1, c) + cnt, c));
            return j_prev;
        }


        //! For each symbol c in wt[i..j-1] get rank(i,c) and rank(j,c).
        /*!
         * \param i        The start index (inclusive) of the interval.
         * \param j        The end index (exclusive) of the interval.
         * \param k        Reference for number of different symbols in [i..j-1].
         * \param cs       Reference to a vector that will contain in
         *                 cs[0..k-1] all symbols that occur in [i..j-1] in
         *                 arbitrary order (if lex_ordered = false) and ascending
         *                 order (if lex_ordered = true).
         * \param rank_c_i Reference to a vector which equals
         *                 rank_c_i[p] = rank(i,cs[p]), for \f$ 0 \leq p < k \f$.
         * \param rank_c_j Reference to a vector which equals
         *                 rank_c_j[p] = rank(j,cs[p]), for \f$ 0 \leq p < k \f$.
         * \par Time complexity
         *      \f$ \Order{\min{\sigma, k \log \sigma}} \f$
         *
         * \par Precondition
         *      \f$ i \leq j \leq size() \f$
         *      \f$ cs.size() \geq \sigma \f$
         *      \f$ rank_{c_i}.size() \geq \sigma \f$
         *      \f$ rank_{c_j}.size() \geq \sigma \f$
         */
        void interval_symbols(size_type i, size_type j, size_type& k,
                              std::vector<value_type>& cs,
                              std::vector<size_type>& rank_c_i,
                              std::vector<size_type>& rank_c_j) const
        {
            assert(i <= j and j <= size());
            if (i==j) {
                k = 0;
            } else if (1==m_sigma) {
                k = 1;
                cs[0] = m_tree.bv_pos_rank(m_tree.root());
                rank_c_i[0] = std::min(i,m_size);
                rank_c_j[0] = std::min(j,m_size);
            } else if ((j-i)==1) {
                k = 1;
                auto rc = inverse_select(i);
                rank_c_i[0] = rc.first; cs[0] = rc.second;
                rank_c_j[0] = rank_c_i[0]+1;
            } else if ((j-i)==2) {
                auto rc = inverse_select(i);
                rank_c_i[0] = rc.first; cs[0] = rc.second;
                rc = inverse_select(i+1);
                rank_c_i[1] = rc.first; cs[1] = rc.second;

                if (cs[0]==cs[1]) {
                    k = 1;
                    rank_c_j[0] = rank_c_i[0]+2;
                } else {
                    k = 2;
                    if (lex_ordered and cs[0] > cs[1]) {
                        std::swap(cs[0], cs[1]);
                        std::swap(rank_c_i[0], rank_c_i[1]);
                    }
                    rank_c_j[0] = rank_c_i[0]+1;
                    rank_c_j[1] = rank_c_i[1]+1;
                }
            } else {
                k = 0;
                _interval_symbols(i, j, k, cs, rank_c_i, rank_c_j, 0);
            }
        }


        //! How many symbols are lexicographic smaller/greater than c in [i..j-1].
        /*!
         * \param i       Start index (inclusive) of the interval.
         * \param j       End index (exclusive) of the interval.
         * \param c       Symbol c.
         * \return A triple containing:
         *         * rank(i,c)
         *         * #symbols smaller than c in [i..j-1]
         *         * #symbols greater than c in [i..j-1]
         *
         * \par Precondition
         *       \f$ i \leq j \leq size() \f$
         * \note
         * This method is only available if lex_ordered = true
         */
        template<class t_ret_type = std::tuple<size_type, size_type, size_type>>
        typename std::enable_if<shape_type::lex_ordered, t_ret_type>::type
        lex_count(size_type i, size_type j, value_type c) const
        {
            assert(i <= j and j <= size());
            if (1==m_sigma) {
                value_type _c = m_tree.bv_pos_rank(m_tree.root());
                if (c == _c) { // c is the only symbol in the wt
                    return t_ret_type {i,0,0};
                } else if (c < _c) {
                    return t_ret_type {0,0,j-i};
                } else {
                    return t_ret_type {0,j-i,0};
                }
            }
            if (i==j) {
                return t_ret_type {rank(i,c),0,0};
            }
            uint64_t p = m_tree.bit_path(c);
            uint32_t path_len = p>>56;
            if (path_len == 0) {  // path_len=0: => c is not present
                value_type _c = (value_type)p;
                if (c == _c) {    // c is smaller than any symbol in wt
                    return t_ret_type {0, 0, j-i};
                }
                auto res = lex_count(i, j, _c);
                return t_ret_type {0, j-i-std::get<2>(res),std::get<2>(res)};
            }
            size_type smaller = 0, greater = 0;
            node_type v = m_tree.root();
            for (uint32_t l=0; l<path_len; ++l, p >>= 1) {
                size_type r1_1 = (m_bv_rank(m_tree.bv_pos(v)+i)
                                  - m_tree.bv_pos_rank(v));
                size_type r1_2 = (m_bv_rank(m_tree.bv_pos(v)+j)
                                  - m_tree.bv_pos_rank(v));

                if (p&1) {
                    smaller += j - r1_2 - i + r1_1;
                    i = r1_1;
                    j = r1_2;
                } else {
                    greater += r1_2 - r1_1;
                    i -= r1_1;
                    j -= r1_2;
                }
                v = m_tree.child(v, p&1);
            }
            return t_ret_type {i, smaller, greater};
        }

        //! How many symbols are lexicographic smaller than c in [0..i-1].
        /*!
         * \param i Exclusive right bound of the range.
         * \param c Symbol c.
         * \return A tuple containing:
         *         * rank(i,c)
         *         * #symbols smaller than c in [0..i-1]
         * \par Precondition
         *       \f$ i \leq size() \f$
         * \note
         * This method is only available if lex_ordered = true
         */
        template<class t_ret_type = std::tuple<size_type, size_type>>
        typename std::enable_if<shape_type::lex_ordered, t_ret_type>::type
        lex_smaller_count(size_type i, value_type c)const
        {
            assert(i <= size());
            if (1==m_sigma) {
                value_type _c = m_tree.bv_pos_rank(m_tree.root());
                if (c == _c) { // c is the only symbol in the wt
                    return t_ret_type {i,0};
                } else if (c < _c) {
                    return t_ret_type {0,0};
                } else {
                    return t_ret_type {0,i};
                }
            }

            uint64_t p = m_tree.bit_path(c);
            uint32_t path_len = p>>56;
            if (path_len == 0) {  // path_len=0: => c is not present
                value_type _c = (value_type)p;
                if (c == _c) {    // c is smaller than any symbol in wt
                    return t_ret_type {0, 0};
                }
                auto res = lex_smaller_count(i, _c);
                return t_ret_type {0, std::get<0>(res)+std::get<1>(res)};
            }
            size_type result = 0;
            size_type all    = i; // possible occurrences of c
            node_type v = m_tree.root();
            for (uint32_t l=0; l<path_len and all; ++l, p >>= 1) {
                size_type ones = (m_bv_rank(m_tree.bv_pos(v)+all)
                                  - m_tree.bv_pos_rank(v));
                if (p&1) {
                    result += all - ones;
                    all    = ones;
                } else {
                    all    -= ones;
                }
                v = m_tree.child(v, p&1);
            }
            return t_ret_type {all, result};
        }

        //! Returns a const_iterator to the first element.
        const_iterator begin()const
        {
            return const_iterator(this, 0);
        }

        //! Returns a const_iterator to the element after the last element.
        const_iterator end()const
        {
            return const_iterator(this, size());
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr,
                            std::string name="") const
        {
            structure_tree_node* child = structure_tree::add_child(
                                             v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_size,out,child, "size");
            written_bytes += write_member(m_sigma,out,child, "sigma");
            written_bytes += m_bv.serialize(out,child,"bv");
            written_bytes += m_bv_rank.serialize(out,child,"bv_rank");
            written_bytes += m_bv_select1.serialize(out,child,"bv_select_1");
            written_bytes += m_bv_select0.serialize(out,child,"bv_select_0");
            written_bytes += m_tree.serialize(out,child,"tree");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in)
        {
            read_member(m_size, in);
            read_member(m_sigma, in);
            m_bv.load(in);
            m_bv_rank.load(in, &m_bv);
            m_bv_select1.load(in, &m_bv);
            m_bv_select0.load(in, &m_bv);
            m_tree.load(in);
        }

        //! Random access container to bitvector of node v
        auto bit_vec(const node_type& v) const -> node_bv_container<t_bitvector> {
            return node_bv_container<t_bitvector>(begin(v), end(v));
        }

        //! Random access container to sequence of node v
        auto seq(const node_type& v) const -> random_access_container<std::function<value_type(size_type)>> {
            return random_access_container<std::function<value_type(size_type)>>([&v, this](size_type i)
            {
                node_type vv = v;
                while (!is_leaf(vv)) {
                    auto vs = expand(vv);
                    auto rs = expand(vv, {0, i});
                    bool bit = *(begin(vv)+i);
                    i = std::get<1>(rs[bit]);
                    vv = vs[bit];
                }
                return sym(vv);
            }, size(v));
        }

        //! Checks if the node is a leaf node
        bool is_leaf(const node_type& v) const
        {
            return m_tree.is_leaf(v);
        }

        //! Symbol for a leaf
        value_type sym(const node_type& v) const
        {
            return m_tree.bv_pos_rank(v);
        }

        //! Indicates if node v is empty
        bool empty(const node_type& v) const
        {
            return size(v)==0;
        }

        //! Return the size of node v
        auto size(const node_type& v) const -> decltype(m_tree.size(v))
        {
            if (is_leaf(v)) {
                if (v == root())
                    return size();
                else {
                    auto parent = m_tree.parent(v);
                    auto rs = expand(parent, {0, size(parent)-1});
                    if (m_tree.child(parent, 0) == v)
                        return std::get<1>(std::get<0>(rs))-std::get<0>((std::get<0>(rs)))+1;
                    else
                        return std::get<1>(std::get<1>(rs))-std::get<0>((std::get<1>(rs)))+1;
                }
            } else {
                return m_tree.size(v);
            }
        }

        //! Returns the root node
        node_type root() const
        {
            return m_tree.root();
        }

        //! Returns the two child nodes of an inner node
        /*! \param v An inner node of a wavelet tree.
         *  \return Return a pair of nodes (left child, right child).
         *  \pre !is_leaf(v)
         */
        std::array<node_type, 2>
        expand(const node_type& v) const
        {
            return {{m_tree.child(v,0), m_tree.child(v,1)}};
        }

        //! Returns for each range its left and right child ranges
        /*! \param v      An inner node of an wavelet tree.
         *  \param ranges A vector of ranges. Each range [s,e]
         *                has to be contained in v=[v_s,v_e].
         *  \return A vector a range pairs. The first element of each
         *          range pair correspond to the original range
         *          mapped to the left child of v; the second element to the
         *          range mapped to the right child of v.
         *  \pre !is_leaf(v) and s>=v_s and e<=v_e
         */
        std::array<range_vec_type, 2>
        expand(const node_type& v,
               const range_vec_type& ranges) const
        {
            auto ranges_copy = ranges;
            return expand(v, std::move(ranges_copy));
        }

        //! Returns for each range its left and right child ranges
        /*! \param v      An inner node of an wavelet tree.
         *  \param ranges A vector of ranges. Each range [s,e]
         *                has to be contained in v=[v_s,v_e].
         *  \return A vector a range pairs. The first element of each
         *          range pair correspond to the original range
         *          mapped to the left child of v; the second element to the
         *          range mapped to the right child of v.
         *  \pre !is_leaf(v) and s>=v_s and e<=v_e
         */
        std::array<range_vec_type, 2>
        expand(const node_type& v,
               range_vec_type&& ranges) const
        {
            auto v_sp_rank = m_tree.bv_pos_rank(v);
            range_vec_type res(ranges.size());
            size_t i = 0;
            for (auto& r : ranges) {
                auto sp_rank    = m_bv_rank(m_tree.bv_pos(v) + r[0]);
                auto right_size = m_bv_rank(m_tree.bv_pos(v) + r[1] + 1)
                                  - sp_rank;
                auto left_size  = (r[1]-r[0]+1)-right_size;

                auto right_sp = sp_rank - v_sp_rank;
                auto left_sp  = r[0] - right_sp;

                r = {left_sp, left_sp + left_size - 1};
                res[i++] = {right_sp, right_sp + right_size - 1};
            }
            return {ranges, std::move(res)};
        }

        //! Returns for a range its left and right child ranges
        /*! \param v An inner node of an wavelet tree.
         *  \param r A ranges [s,e], such that [s,e] is
         *           contained in v=[v_s,v_e].
         *  \return A range pair. The first element of the
         *          range pair correspond to the original range
         *          mapped to the left child of v; the second element to the
         *          range mapped to the right child of v.
         *  \pre !is_leaf(v) and s>=v_s and e<=v_e
         */
        std::array<range_type, 2>
        expand(const node_type& v, const range_type& r) const
        {
            auto v_sp_rank = m_tree.bv_pos_rank(v);
            auto sp_rank    = m_bv_rank(m_tree.bv_pos(v) + r[0]);
            auto right_size = m_bv_rank(m_tree.bv_pos(v) + r[1] + 1)
                              - sp_rank;
            auto left_size  = (r[1]-r[0]+1)-right_size;

            auto right_sp = sp_rank - v_sp_rank;
            auto left_sp  = r[0] - right_sp;

            return {{{{left_sp, left_sp + left_size - 1}},
                    {{right_sp, right_sp + right_size - 1}}
                }
            };
        }

        //! return the path to the leaf for a given symbol
        std::pair<uint64_t,uint64_t> path(value_type c) const
        {
            uint64_t path = m_tree.bit_path(c);
            uint64_t path_len = path >> 56;
            // reverse the path till we fix the ordering
            path = bits::rev(path);
            path = path >> (64-path_len); // remove the length
            return {path_len,path};
        }

        //! Returns for a symbol c the next larger or equal symbol in the WT.
        /*! \param c the symbol
         *  \return A pair. The first element of the pair consititues if
         *          a valid answer was found (true) or no valid answer (false)
         *          could be found. The second element contains the found symbol.
         */
        std::pair<bool, value_type> symbol_gte(value_type c) const
        {
            return m_tree.symbol_gte(c);
        }

        //! Returns for a symbol c the previous smaller or equal symbol in the WT.
        /*! \param c the symbol
         *  \return A pair. The first element of the pair consititues if
         *          a valid answer was found (true) or no valid answer (false)
         *          could be found. The second element contains the found symbol.
         */
        std::pair<bool, value_type> symbol_lte(value_type c) const
        {
            return m_tree.symbol_lte(c);
        }

    private:

        //! Iterator to the begin of the bitvector of inner node v
        auto begin(const node_type& v) const -> decltype(m_bv.begin() + m_tree.bv_pos(v))
        {
            return m_bv.begin() + m_tree.bv_pos(v);
        }

        //! Iterator to the begin of the bitvector of inner node v
        auto end(const node_type& v) const -> decltype(m_bv.begin() + m_tree.bv_pos(v) + m_tree.size(v))
        {
            return m_bv.begin() + m_tree.bv_pos(v) + m_tree.size(v);
        }
};

}

#endif
