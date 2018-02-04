/* sdsl - succinct data structures library
 Copyright (C) 2010 Simon Gog

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
/*! \file StreeOhleb.hpp
 \brief StreeOhleb.hpp contains an implementation of the interval based CST.
 \author Simon Gog
 */
#ifndef stree_ohleb
#define stree_ohleb

#include <iostream>
#include <cassert>


#include <sdsl/int_vector.hpp>
#include <sdsl/suffix_trees.hpp>



namespace fdms
{

    // Declaration of the CST's node type
    template<class t_int = sdsl::int_vector<>::size_type>
    struct bp_interval;

    //! A class for the Compressed Suffix Tree (CST) proposed by Ohlebusch and Gog.
    /*!
     * \tparam t_csa       Type of a CSA (member of this type is accessible via
     *                     member `csa`, default class is sdsl::csa_sada).
     * \tparam t_lcp       Type of a LCP structure (member is accessible via member
     *                     `lcp`, default class is sdsl::lcp_support_sada),
     * \tparam t_bp_support Type of a BPS structure (member accessible via member
     *                      `bp_support`, default class is sdsl::bp_support_sada),
     * \tparam t_rank       Type of rank structure which supports the bitvector
     *                      which indicates the leftmost child of the nodes.
     *
     * It also contains a sdsl::bit_vector which represents the BP sequence of the
     * Super-Cartesian tree of the LCP array. This bitvector can be accessed via
     * the member `bp`. Another sdsl::bit_vector stores information, if a node is
     * the leftmost child of another node. This bitvector can be access via the
     * member first_child_bv and takes n bits.
     *
     * A node \f$v\f$ of the csa_sct is represented by an sdsl::bp_interval. The
     * size of the sdsl::StreeOhleb is smaller than the size of a sdsl::cst_sada
     * since the tree topology needs only \f$2n+n=3n\f$ bits in contrast to the
     * \f$4n\f$ bits in sdsl::cst_sada.
     *
     * \par Reference
     * Enno Ohlebusch, Johannes Fischer, Simon Gog:
     * CST++.
     * SPIRE 2010: 322-333
     *
     * \par Applications of the CST
     * The compressed suffix tree could be used for string matching and many other
     * application in sequence analysis. 17 applications are in the book
     * "Algorithms on Strings, Trees, and Sequences" of Dan Gusfield.
     *
     * @ingroup cst
     */
    template<class t_csa = sdsl::csa_wt<>,
             class t_lcp = sdsl::lcp_dac<>,
             class t_bp_support = sdsl::bp_support_sada<>,
             class t_bv = sdsl::bit_vector,
             class t_rank = typename std::conditional<std::is_same<t_bv, sdsl::bit_vector>::value, sdsl::rank_support_v5<>, typename t_bv::rank_1_type>::type,
             class t_sel =  typename std::conditional<std::is_same<t_bv, sdsl::bit_vector>::value and std::is_same<typename t_csa::alphabet_category, sdsl::byte_alphabet_tag>::value, sdsl::select_support_scan<>, typename t_bv::select_1_type>::type>
    class StreeOhleb
    {
        static_assert(std::is_same<typename sdsl::index_tag<t_csa>::type, sdsl::csa_tag>::value,
                      "First template argument has to be a compressed suffix array.");
    public:
        typedef typename t_csa::size_type                      size_type;
        typedef ptrdiff_t                                      difference_type;
        typedef t_csa                                          csa_type;
        typedef typename t_lcp::template type<StreeOhleb>      lcp_type;
        typedef t_bp_support                                   bp_support_type;
        typedef typename t_csa::char_type                      char_type;
        typedef typename t_csa::string_type                    string_type;
        typedef bp_interval<size_type>                         node_type; //!< Type for the nodes in the tree
        typedef t_bv                                           bv_type;
        typedef t_rank                                         rank_type;
        typedef t_sel                                          sel_type;

        typedef typename t_csa::alphabet_type::comp_char_type  comp_char_type;
        typedef typename t_csa::alphabet_type::sigma_type      sigma_type;

        typedef typename t_csa::alphabet_category              alphabet_category;
        typedef sdsl::cst_tag                                  index_category;
        typedef sdsl::bit_vector                               bit_vector;
        csa_type        m_csa;

    private:

        // for the Super Cartesian Tree
        bit_vector      m_bp;
        bp_support_type m_bp_support;

        bv_type         m_first_child;
        rank_type       m_first_child_rank;
        sel_type        m_first_child_select;
        size_type       m_nodes;

        void copy(const StreeOhleb& cst)
        {
            m_csa              = cst.m_csa;
            m_bp               = cst.m_bp;
            m_bp_support       = cst.m_bp_support;
            m_bp_support.set_vector(&m_bp);
            m_first_child      = cst.m_first_child;
            m_first_child_rank = cst.m_first_child_rank;
            m_first_child_rank.set_vector(&m_first_child);
            m_first_child_select = cst.m_first_child_select;
            m_first_child_select.set_vector(&m_first_child);
            m_nodes            = cst.m_nodes;
        }

        // Get the next smaller value.
        /*
         * \param i    Position in the original vector.
         * \param ipos Position of the corresponding opening parenthesis in BP.
         * \return Position of the next smaller value in [i+1..n-1], and n when
         *         no such value exists.
         * \par Time complexity
         *      \f$ \Order{1} \f$
         */
        // possible optimization: calculate also position of nsv,
        // i.e. next ( following position cipos
        inline size_type nsv(SDSL_UNUSED size_type i, size_type ipos)const
        {
            size_type cipos = m_bp_support.find_close(ipos);
            size_type result = m_bp_support.rank(cipos);
            return result;
        }

        // Get the previous smaller value.
        /*
         * \param i      Position in the original vector.
         * \param ipos   Corresponding opening parenthesis in m_bp
         * \param cipos  Corresponding closing parenthesis to ipos
         * \par Time complexity
         *    \f$ \Order{\frac{\sigma}{w}} \f$, where w=64 is the word size,
         *    can be implemented in \f$\Order{1}\f$ with rank and select.
         */
        inline size_type psv(SDSL_UNUSED size_type i, size_type ipos,
                             size_type cipos, size_type& psvpos,
                             size_type& psvcpos)const
        {
            // if lcp[i]==0 => psv is the 0-th index by definition
            if ((cipos + (size_type)m_csa.sigma) >= m_bp.size()) {
                psvpos = 0;
                psvcpos = m_bp.size()-1;
                return 0;
            }
            if (m_bp[cipos+1]) {
                psvpos = m_bp_support.enclose(ipos);
                psvcpos = m_bp_support.find_close(psvpos);
                return m_bp_support.rank(psvpos)-1;
            }
            // r0 = index of clothing parenthesis in m_first_child
            size_type r0 = cipos - m_bp_support.rank(cipos);
            size_type next_first_child = 0;
            const uint64_t* p = m_first_child.data() + (r0>>6);
            uint64_t w = (*p) >> (r0&0x3F);
            if (w) { // if w!=0
                next_first_child = cipos + sdsl::bits::lo(w);
                if (cipos == next_first_child and m_bp[next_first_child+1]) {
                    psvpos = m_bp_support.enclose(ipos);
                    psvcpos = m_bp_support.find_close(psvpos);
                    return m_bp_support.rank(psvpos)-1;
                }
            } else {
                size_type delta = 63-(r0&0x3F);
                ++p;
                int steps = 4;
                while (!(w=*p) and steps-- > 0) { // while w==0
                    ++p;
                    delta += 64;
                }
                if (w != 0) {
                    delta += sdsl::bits::lo(w) + 1;
                } else {
                    auto pos = m_first_child_select(m_first_child_rank(r0+1)+1);
                    delta    = pos - r0;
                }
                next_first_child = cipos + delta;
            }
            if (!m_bp[next_first_child+1]) { // if next parenthesis is a closing one
                psvcpos = next_first_child+1;
                psvpos = m_bp_support.find_open(psvcpos);
                return m_bp_support.rank(psvpos)-1;
            } else {
                psvpos = m_bp_support.enclose(m_bp_support.find_open(next_first_child));
                psvcpos = m_bp_support.find_close(psvpos);
                return m_bp_support.rank(psvpos)-1;
            }
        }

        // Range minimum query based on the rr_enclose method.
        /* \par Time complexity
         *   \f$ \Order{\rrenclose} \f$
         */
        inline size_type rmq(size_type l, size_type r)const
        {
            size_type i     = m_bp_support.select(l+1);
            size_type j     = m_bp_support.select(r+1);
            size_type fc_i     = m_bp_support.find_close(i);
            if (j < fc_i) { // i < j < find_close(j) < find_close(i)
                return l;
            } else { // i < find_close(i) < j < find_close(j)
                size_type ec = m_bp_support.rr_enclose(i,j);
                if (ec == m_bp_support.size()) {// no restricted enclosing pair found
                    return r;
                } else { // found range restricted enclosing pair
                    return m_bp_support.rank(ec)-1; // subtract 1, as the index is 0 based
                }
            }
        }

        // Get the first l-index of a node
        // if there exists no ith l-index return node.j+1
        /* \param v Node
         * \param i l-index in [1..degree()]
         * \paran
         */
        size_type select_l_index(const node_type& v, size_type& kpos, size_type& ckpos)const
        {
            if (v.cipos > v.jp1pos) { // corresponds to m_lcp[i] <= m_lcp[j+1]
                ckpos    = v.jp1pos-1;
            } else { // corresponds to m_lcp[i] > m_lcp[j+1]
                ckpos    = v.cipos-1;
            }

            assert(m_bp[ckpos] == 0);   // at least the first l-index should be present, i.e. node is not leaf
            kpos    = m_bp_support.find_open(ckpos);
            return m_bp_support.rank(kpos) - 1;
        }

    public:
        const csa_type&             csa              = m_csa;
        const bit_vector&           bp               = m_bp;
        const bp_support_type&      bp_support       = m_bp_support;

        const bv_type&   first_child_bv     = m_first_child;
        const rank_type& first_child_rank   = m_first_child_rank;
        const sel_type&  first_child_select = m_first_child_select;

        /*! \defgroup StreeOhleb_constructors Constructors of StreeOhleb */
        /* @{ */

        //! Default constructor
        StreeOhleb() {}

        //! Construct CST from cache config
        StreeOhleb(sdsl::cache_config& cache, bool build_only_bps=false);

        //! Copy constructor
        /*!
         *  \param cst The StreeOhleb which should be copied.
         *  \par Time complexity
         *       \f$ \Order{n} \f$, where \f$n=\f$StreeOhleb.size()
         */
        StreeOhleb(const StreeOhleb& cst)
        {
            copy(cst);
        }

        //! Move constructor
        /*!
         *  \param cst The StreeOhleb which should be moved.
         */
        StreeOhleb(StreeOhleb&& cst)
        {
            *this = std::move(cst);
        }

        /* @} */

        //! Number of leaves of the suffix tree.
        /*! Required for the Container Concept of the STL.
         *  \sa max_size, empty
         */
        size_type size()const
        {
            return m_bp.size()>>1;
        }

        //! Returns the largest size that StreeOhleb can ever have.
        /*! Required for the Container Concept of the STL.
         *  \sa size
         */
        static size_type max_size()
        {
            return t_csa::max_size();
        }

        //! Returns if the data structure is empty.
        /*! Required for the Container Concept of the STL.
         * \sa size
         */
        bool empty()const
        {
            return m_csa.empty();
        }

        //! Swap method for StreeOhleb
        /*! The swap method can be defined in terms of assignment.
         This requires three assignments, each of which, for a container type, is linear
         in the container's size. In a sense, then, a.swap(b) is redundant.
         This implementation guaranties a run-time complexity that is constant rather than linear.
         \param cst StreeOhleb to swap.

         Required for the Assignable Conecpt of the STL.
         */
        void swap(StreeOhleb& cst)
        {
            if (this != &cst) {
                m_csa.swap(cst.m_csa);
                m_bp.swap(cst.m_bp);
                sdsl::util::swap_support(m_bp_support, cst.m_bp_support, &m_bp, &(cst.m_bp));
                m_first_child.swap(cst.m_first_child);
                sdsl::util::swap_support(m_first_child_rank, cst.m_first_child_rank, &m_first_child, &(cst.m_first_child));
                sdsl::util::swap_support(m_first_child_select, cst.m_first_child_select, &m_first_child, &(cst.m_first_child));
                std::swap(m_nodes, cst.m_nodes);
            }
        }

        //! Assignment Operator.
        /*!
         *    Required for the Assignable Concept of the STL.
         */
        StreeOhleb& operator=(const StreeOhleb& cst);

        //! Assignment Move Operator.
        /*!
         *    Required for the Assignable Concept of the STL.
         */
        StreeOhleb& operator=(StreeOhleb&& cst);

        //! Serialize to a stream.
        /*! \param out Outstream to write the data structure.
         *  \return The number of written bytes.
         */
        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")const;

        //! Load from a stream.
        /*! \param in Inputstream to load the data structure from.
         */
        void load(std::istream& in);

        /*! \defgroup StreeOhleb_tree_methods Tree methods of StreeOhleb */
        /* @{ */

        //! Return the root of the suffix tree.
        /*!
         * \par Time complexity O(1)
         *      \f$ \Order{1} \f$
         */
        node_type root() const
        {
            return node_type(0, size()-1, 0, m_bp.size()-1, m_bp.size());
        }

        bool is_root(const node_type& v)const {
            return (v.i == 0 &&
                    v.j == size() - 1 &&
                    v.ipos == 0 &&
                    v.cipos == m_bp.size() - 1 &&
                    v.jp1pos == m_bp.size());
        }

        //! Decide if a node is a leaf.
        /*!
         * \param v A valid node.
         * \returns A boolean value indicating if v is a leaf.
         * \par Time complexity
         *      \f$ \Order{1} \f$
         */
        bool is_leaf(const node_type& v)const
        {
            return v.i==v.j;
        }

        //! Return the i-th leaf (1-based from left to right).
        /*!
         * \param i 1-based position of the leaf.
         * \return The i-th leave.
         * \par Time complexity
         *      \f$ \Order{1} \f$
         *       \f$ 1 \leq i \leq size() \f$
         */
        node_type select_leaf(size_type i)const
        {
            assert(i > 0 and i <= size());
            size_type ipos = m_bp_support.select(i);
            size_type jp1pos;
            if (i == size())
                jp1pos = m_bp.size();
            else if (m_bp[ipos+1])
                jp1pos = ipos+1;
            else
                jp1pos = m_bp_support.select(i+1);
            return node_type(i-1, i-1, ipos, m_bp_support.find_close(ipos), jp1pos);
        }

        //! Calculate the number of leaves in the subtree rooted at node v.
        /*! \param v A valid node of the suffix tree.
         *  \return The number of leaves in the subtree rooted at node v.
         *  \par Time complexity
         *    \f$ \Order{1} \f$
         */
        size_type size(const node_type& v)const
        {
            return v.j-v.i+1;
        }

        //! Calculates the leftmost leaf in the subtree rooted at node v.
        /*! \param v A valid node of the suffix tree.
         *  \return The leftmost leaf in the subtree rooted at node v.
         *  \par Time complexity
         *    \f$ \Order{1} \f$
         */
        node_type leftmost_leaf(const node_type& v)const
        {
            return select_leaf(v.i+1);
        }

        //! Calculates the rightmost leaf in the subtree rooted at node v.
        /*! \param v A valid node of the suffix tree.
         *  \return The rightmost leaf in the subtree rooted at node v.
         *  \par Time complexity
         *    \f$ \Order{1} \f$
         */
        node_type rightmost_leaf(const node_type& v)const
        {
            return select_leaf(v.j+1);
        }

        //! Calculates the index of the leftmost leaf in the corresponding suffix array.
        /*! \param v A valid node of the suffix tree.
         *  \return The index of the leftmost leaf in the corresponding suffix array.
         *  \par Time complexity
         *    \f$ \Order{1} \f$
         *  \par Note
         *  lb is an abbreviation for ,,left bound''
         */
        size_type lb(const node_type& v)const
        {
            return v.i;
        }

        //! Calculates the index of the rightmost leaf in the corresponding suffix array.
        /*! \param v A valid node of the suffix tree.
         *     \return The index of the rightmost leaf in the corresponding suffix array.
         *    \par Time complexity
         *        \f$ \Order{1} \f$
         *  \par Note
         *   rb is an abbreviation for ,,right bound''
         */
        size_type rb(const node_type& v)const
        {
            return v.j;
        }

        //! Calculate the LCA of two nodes `v` and `w`
        /*!
         * \param v The first node.
         * \param w The second node.
         * \return The lowest common ancestor of v and w.
         * \par Time complexity
         *   \f$ \Order{\rrenclose}\   \f$
         */

        node_type lca(node_type v, node_type w)const
        {
            if (v.i > w.i or(v.i == w.i and v.j < w.j)) {
                std::swap(v, w);
            }
            if (v.j >= w.j) { // v encloses w or v==w
                return v;
            } else { // v.i < v.j < w.i < w.j
                size_type min_index = rmq(v.i+1, w.j);
                size_type min_index_pos     = m_bp_support.select(min_index+1);
                size_type min_index_cpos     = m_bp_support.find_close(min_index_pos);

                if (min_index_cpos >= (m_bp.size() - m_csa.sigma)) {   // if lcp[min_index]==0 => return root
                    return root();
                }
                size_type new_j = nsv(min_index, min_index_pos)-1;
                size_type new_ipos, new_icpos;
                size_type new_i = psv(min_index, min_index_pos, min_index_cpos, new_ipos, new_icpos);
                size_type jp1pos = m_bp.size();
                if (new_j < size()-1) {
                    jp1pos = m_bp_support.select(new_j+2);
                }
                return node_type(new_i, new_j, new_ipos, new_icpos, jp1pos);
            }
        }
/*
        bool has_wl(const node_type v, const char_type c) const {
            //node_type u = double_rank_fail_wl(v, c);
            std::pair<size_type, size_type> I = m_csa.bwt.double_rank_and_fail(v.i, v.j  + 1, c);
            int cc = m_csa.char2comp[c];
            I.first += m_csa.C[cc];
            I.second += (m_csa.C[cc] - 1);
            return (m_csa.C[cc] > 0) and (I.first <= I.second);
        }
*/
        /*
         * call parent(v) in sequece until reaching a node u for which wl(u, c) exists
         */
        node_type parent_sequence(const node_type& v, const char_type c) const {
            std::pair<size_type, size_type> I = std::make_pair(0, 0);
            size_type cc = m_csa.char2comp[c];
            node_type vv = v, u = root();

            if(!has_complete_info(vv))
                lazy_wl_followup(vv);

            bool has_wl = false;
            do{ // remove suffixes of t[k..] until you can extend by 'c'
                vv = parent(vv);
                u = single_rank_wl(vv, c);
                has_wl = !is_root(u);
            } while(!has_wl && !is_root(vv));

            /*
            do{ // remove suffixes of t[k..] until you can extend by 'c'
                vv = parent(vv);
                I.first = m_csa.C[cc] + m_csa.bwt.rank(vv.i, c);
                I.second = m_csa.C[cc] + m_csa.bwt.rank(vv.j + 1, c) - 1;
            } while(I.first > I.second);
            */
            return vv;
        }

        /*
         * USE select_at_dist()
         * call parent(v) in sequece until reaching a node u for which wl(u, c) exists
         */
        node_type maxrep_ancestor(const node_type& v, const char_type c) const {
            //return parent_sequence(v, c);
            //return _maxrep_ancestor(v, c);
            
            size_type cnt_c = m_csa.C[m_csa.char2comp[c] + 1] - m_csa.C[m_csa.char2comp[c]];
  
            //index of first occurrence of c after position v.j
            size_type r = (m_csa.bwt.rank(v.j + 1, c) < cnt_c ? m_csa.bwt.select_at_dist(c, v.j, 1) : size());
            node_type p = r < size() ? lca(v, select_leaf(r + 1)) : root();

            if(p.i == v.i)
                return p;
            
            //index of last occurrence of c before position v.i
            size_type l = (m_csa.bwt.rank(v.i, c) > 0 ? m_csa.bwt.select_at_dist(c, v.i, 0) : 0);
            if(p.i > l)
                return p;
            node_type q = l > 0 ? lca(select_leaf(l + (l < size())), v) : root();
            
            // computing lca(p, q)
            node_type res = (q.j - q.i <= p.j - p.i ? q : p);
            //node_type exp_res = _maxrep_ancestor(v, c);
            //if(res.i != exp_res.i or res.j != exp_res.j)
            //    assert(0);
            //assert (res == _maxrep_ancestor(v, c));
            return res;
        }
        

        /*
         * call parent(v) in sequece until reaching a node u for which wl(u, c) exists
         * TODO: this fails with c in t but c notin s
         */
        node_type _maxrep_ancestor(const node_type& v_, const char_type c) const {
            node_type u, p, q, v = v_;

            // ** not needed **/
            //if(!has_complete_info(v))
            //    lazy_wl_followup(v);

            std::pair<size_type, size_type> left_right_cnt_c = m_csa.bwt.double_rank(v.i, v.j + 1, c);
            size_type cnt_c = m_csa.C[m_csa.char2comp[c] + 1] - m_csa.C[m_csa.char2comp[c]];

            // computing p
            size_type right_cnt_c = left_right_cnt_c.second;
            if(right_cnt_c < cnt_c){
                // r = select(rank(v.i) + 1, c) = select_next(v.i, c) = select_at_dist(c, v.i, 1)
                size_type r = m_csa.bwt.select(right_cnt_c + 1, c);
                u = select_leaf(r + 1);
                p = lca(v, u);
                if(p.i == v.i)
                    return p;
            } else
                p = root();

            // computing q
            size_type left_cnt_c = left_right_cnt_c.first;
            if(left_cnt_c > 0){
                // l = select(rank(v.j + 1) + 1, c) = select_next(v.j + 1, c) = select_at_dist(c, v.j + 1, 1)

                size_type l = m_csa.bwt.select(left_cnt_c, c);
                if(p.i > l)
                    return p;
                u = select_leaf(l + 1);
                q = lca(u, v);
            } else
                q = root();

            // computing lca(p, q)
            if(q.j - q.i <= p.j - p.i)
                return q;
            return p;
        }

        //! Calculate the parent node of a node v.
        /*! \param v A valid node of the suffix tree.
         *  \return The parent node of v or the root if v==root().
         *  \par Time complexity
         *     \f$ \Order{1}\f$
         */
        node_type parent(const node_type& v) const
        {
            if (v.cipos > v.jp1pos) { // LCP[i] <= LCP[j+1]
                size_type psv_pos, psv_cpos, psv_v, nsv_v, nsv_p1pos;
                psv_v = psv(v.j+1, v.jp1pos, m_bp_support.find_close(v.jp1pos), psv_pos, psv_cpos);
                nsv_v = nsv(v.j+1, v.jp1pos)-1;
                if (nsv_v == size()-1) {
                    nsv_p1pos = m_bp.size();
                } else { // nsv_v < size()-1
                    nsv_p1pos = m_bp_support.select(nsv_v+2);
                }
                return node_type(psv_v, nsv_v, psv_pos, psv_cpos, nsv_p1pos);
            } else { // LCP[i] > LCP[j+1]
                size_type psv_pos, psv_cpos, psv_v;
                psv_v = psv(v.i, v.ipos, v.cipos, psv_pos, psv_cpos);
                return node_type(psv_v, v.j, psv_pos, psv_cpos, v.jp1pos);
            }
        }

        //! Returns the next sibling of node v.
        /*!
         * \param v A valid node v of the suffix tree.
         * \return The next (right) sibling of node v or root() if v has no next (right) sibling.
         * \par Time complexity
         *   \f$ \Order{1} \f$
         */
        node_type sibling(const node_type& v)const
        {
            //Procedure:(1) Determine, if v has a right sibling.
            if (v.cipos < v.jp1pos) { // LCP[i] > LCP[j+1] => v has the same right border as parent(v) => no right sibling
                return root();
            }
            //Procedure:(2) There exists a right sibling, LCP[j+1] >= LCP[i] and j>i
            // Now it holds:  v.cipos > v.jp1pos
            size_type cjp1posm1 = m_bp_support.find_close(v.jp1pos)-1; // v.cipos-2 ???
            // m_bp[cjp1posm1] equals 1 =>  v is the last child
            bool last_child = m_bp[cjp1posm1];
            // otherwise if m_bp[cjp1posm1] equals 0 => we don't know if it is the last child
            if (!last_child) {
                size_type first_child_idx = cjp1posm1 - m_bp_support.rank(cjp1posm1);
                last_child = m_first_child[first_child_idx]; // if first_child indicator is true => the new sibling is the rightmost sibling
            }
            if (last_child) {
                size_type nsv_v = nsv(v.j+1, v.jp1pos)-1, nsv_p1pos;
                if (nsv_v == size()-1) {
                    nsv_p1pos = m_bp.size();
                } else {
                    nsv_p1pos = m_bp_support.select(nsv_v+2);
                }
                return node_type(v.j+1, nsv_v, v.jp1pos, m_bp_support.find_close(v.jp1pos), nsv_p1pos);
            } else {
                size_type new_j = m_bp_support.rank(m_bp_support.find_open(cjp1posm1))-2;
                return node_type(v.j+1, new_j, v.jp1pos, m_bp_support.find_close(v.jp1pos), m_bp_support.select(new_j+2));
            }
        }

        //! Get the first child of a node v.
        /*!
         * \param v A valid tree node of the cst.
         * \return The first child node of v or root() if v has no children.
         * \par Time complexity
         * \f$ \Order{\frac{\sigma}{w}} \f$, where w=64 is the word size,
         * can be implemented in \f$\Order{1}\f$ with rank and select.
         *      \f$ 1 \leq i \leq degree(v) \f$
         */
        node_type first_child(const node_type& v) const
        {
            if (is_leaf(v))  // if v is a leave, v has no child
                return root();

            // v is not a leaf: v has at least two children
            size_type k = 0, kpos = 0, k_find_close = 0;
            k = select_l_index(v, kpos, k_find_close);// get first l-index k and the position of k
            return node_type(v.i, k-1, v.ipos, v.cipos, kpos);
        }

        //! Compute the suffix link of node v.
        /*!
         * \param v A valid node of a StreeOhleb.
         * \return The suffix link of node v.
         * \par Time complexity
         *      \f$ \Order{ \rrenclose } \f$
         */
        node_type sl(const node_type& v)const
        {
            if (v == root())
                return root();
            // get interval with first char deleted
            size_type i     = m_csa.psi[v.i];
            if (is_leaf(v)) {
                if (v.i==0 and v.j==0) // if( v.l==1 )
                    return root();
                else
                    return select_leaf(i+1);
            }
            size_type j     = m_csa.psi[v.j];
            assert(i < j);
            size_type min_index = rmq(i+1, j); // rmq
            size_type min_index_pos     = m_bp_support.select(min_index+1);
            size_type min_index_cpos     = m_bp_support.find_close(min_index_pos);
            if (min_index_cpos >= (m_bp.size() - m_csa.sigma)) {  // if lcp[min_index]==0 => return root
                return root();
            }
            size_type new_j = nsv(min_index, min_index_pos)-1;
            size_type new_ipos, new_icpos;
            size_type new_i = psv(min_index, min_index_pos, min_index_cpos, new_ipos, new_icpos);
            size_type jp1pos = m_bp.size();
            if (new_j < size()-1) {
                jp1pos = m_bp_support.select(new_j+2);
            }
            return node_type(new_i, new_j, new_ipos, new_icpos, jp1pos);
        }

        node_type single_rank_wl(const node_type& v, const char_type c) const
        {
#ifdef VERBOSE
            cerr  << "single_rank_wl" << endl;
#endif
            size_type c_left    = m_csa.bwt.rank(v.i, c);
            size_type c_right   = m_csa.bwt.rank(v.j+1, c);
            return _wl_from_interval(std::make_pair(c_left, c_right), c);
        }

        //! Compute the Weiner link of node v and character c.
        /*!
         * \param v A valid node of a StreeOhleb.
         * \param c The character which should be prepended to the string of the current node.
         * \return  root() if the Weiner link of (v, c) does not exist,
         *          otherwise the Weiner link is returned.
         * \par Time complexity
         *        \f$ \Order{ t_{rank\_bwt} } \f$
         */
        node_type double_rank_nofail_wl(const node_type& v, const char_type c) const
        {
#ifdef VERBOSE
            cerr  << "double_rank_nofail_wl" << endl;
#endif
            // what in single_rank_wl is (c_left, c_right)
            std::pair<size_type, size_type> lr = m_csa.bwt.double_rank(v.i, v.j+1, c);
            return _wl_from_interval(lr, c);
        }

        // as the above, but fail's if early if Weiner Link doesn't exist
        node_type double_rank_fail_wl(const node_type& v, const char_type c) const
        {
#ifdef VERBOSE
            cerr  << "double_rank_fail_wl" << endl;
#endif
            // what in single_rank_wl is (c_left, c_right)
            std::pair<size_type, size_type> lr = m_csa.bwt.double_rank_and_fail(v.i, v.j+1, c);
            return _wl_from_interval(lr, c);
        }

        /*
          if non-maxrep, check and return after a single rank
        */
        node_type double_rank_fail_wl_mrep_vanilla(const node_type& v, const char_type c, const bool is_maximal) const
        {
#ifdef VERBOSE
            cerr  << "double_rank_fail_wl_mrep_vanilla" << endl;
#endif
            if(is_maximal)
                return double_rank_fail_wl(v, c);
            
            if (m_csa.bwt[v.j] != c)
                return root();
            size_type c_right = m_csa.bwt.rank(v.j + 1, c);

            std::pair<size_type, size_type> lr = std::make_pair(c_right - (v.j - v.i + 1), c_right);
            return _wl_from_interval(lr, c);
        }

        /*
          if non-maxrep, rank_and_check and return
        */
        node_type double_rank_fail_wl_mrep_rc(const node_type& v, const char_type c, const bool is_maximal) const
        {
#ifdef VERBOSE
            cerr  << "double_rank_fail_wl_mrep_rc" << endl;
#endif
            if(is_maximal)
                return double_rank_fail_wl(v, c);
            
            size_type c_right = m_csa.bwt.rank_and_check(v.j + 1, c);
            if(c_right == 0)
                return root();
            
            std::pair<size_type, size_type> lr = std::make_pair(c_right - (v.j - v.i + 1), c_right);
            return _wl_from_interval(lr, c);
        }

        //! Weiner link instructions common to wl() all implementations listed above
        /*!
         * \param lr An interval containing the rank information
         * \param c The character that should be prepended to the string of the current node.
         * \return  root() if the Weiner link of (v, c) does not exist,
         *          otherwise the Weiner link is returned.
         */
        node_type _wl_from_interval(const std::pair<size_type, size_type> &lr, const char_type c) const {
            if (lr.first == lr.second)  // there exists no Weiner link
                return root();

            if (lr.first+1 == lr.second)
                return select_leaf(m_csa.C[m_csa.char2comp[c]] + lr.first + 1);

            size_type left     = m_csa.C[m_csa.char2comp[c]] + lr.first;
            size_type right    = m_csa.C[m_csa.char2comp[c]] + lr.second - 1;
            assert(left < right);

            size_type ipos   = m_bp_support.select(left+1);
            size_type jp1pos = m_bp.size();
            if (right < size()-1) {
                jp1pos = m_bp_support.select(right+2);
            }
            return node_type(left, right, ipos,
                             m_bp_support.find_close(ipos), jp1pos);
        }
        
        /*! The Weiner Link of the given node, with only the fields needed to compute it's wl.
         * \param v A valid node of a StreeOhleb.
         * \param c The character which should be prepended to the string of the current node.
         * \return  root() if the Weiner link of (v, c) does not exist,
         *          otherwise the Weiner link is returned.
         *  \par Time complexity
         *        \f$ \Order{ t_{rank\_bwt} } \f$
         */
        node_type lazy_single_rank_wl(const node_type& v, const char_type c) const
        {
#ifdef VERBOSE
            cerr  << "lazy_single_rank_wl" << endl;
#endif
            size_type c_left    = m_csa.bwt.rank(v.i, c);
            size_type c_right   = m_csa.bwt.rank(v.j+1, c);
            return _lazy_wl_from_interval(std::make_pair(c_left, c_right), c);
        }

        /*! The Weiner Link of the given node, with only the fields needed to compute it's wl.
         * \param v A valid node of a StreeOhleb.
         * \param c The character which should be prepended to the string of the current node.
         * \return  root() if the Weiner link of (v, c) does not exist,
         *          otherwise the Weiner link is returned.
         *  \par Time complexity
         *        \f$ \Order{ t_{rank\_bwt} } \f$
         */
        node_type lazy_double_rank_nofail_wl(const node_type& v, const char_type c) const
        {
#ifdef VERBOSE
            cerr  << "lazy_double_rank_nofail_wl" << endl;
#endif
            // what in single_rank_wl is (c_left, c_right)
            std::pair<size_type, size_type> lr = m_csa.bwt.double_rank(v.i, v.j+1, c);
            return _lazy_wl_from_interval(lr, c);
        }

        /*! The Weiner Link of the given node, with only the fields needed to compute it's wl.
         * \param v A valid node of a StreeOhleb.
         * \param c The character which should be prepended to the string of the current node.
         * \return  root() if the Weiner link of (v, c) does not exist,
         *          otherwise the Weiner link is returned.
         *  \par Time complexity
         *        \f$ \Order{ t_{rank\_bwt} } \f$
         */
        node_type lazy_double_rank_fail_wl(const node_type& v, const char_type c) const
        {
#ifdef VERBOSE
            cerr  << "lazy_double_rank_fail_wl" << endl;
#endif
            std::pair<size_type, size_type> lr = m_csa.bwt.double_rank_and_fail(v.i, v.j+1, c);
            return _lazy_wl_from_interval(lr, c);
        }

        bool has_complete_info(const node_type v) const {
            return !((v.ipos == 0) && (v.cipos == 0) && (v.jp1pos == 0));
        }

        //! Weiner link instructions common to lazy_wl() all implementations listed above
        /*!
         * \param lr An interval containing the rank information
         * \param c The character that should be prepended to the string of the current node.
         * \return  root() if the Weiner link of (v, c) does not exist,
         *          otherwise the Weiner link is returned.
         */
        node_type _lazy_wl_from_interval(const std::pair<size_type, size_type> &lr, const char_type c) const {
            if (lr.first == lr.second)  // there exists no Weiner link
                return root();

            size_type left = m_csa.C[m_csa.char2comp[c]] + lr.first;
            if (lr.first + 1 == lr.second)
                return node_type(left, left, 0, 0, 0);

            size_type right    = m_csa.C[m_csa.char2comp[c]] + lr.second - 1;
            assert(left < right);
            return node_type(left, right, 0, 0, 0);
        }

        //! Complete the lazy_wl call on the node
        /*!
         * \param v A valid node returned by a call to lazy_wl()
         */
        void lazy_wl_followup(node_type& v) const {
            //size_type left = v.i;
            //size_type right = v.j;

            size_type ipos = m_bp_support.select(v.i + 1);
            size_type jp1pos = m_bp.size();
            if (v.j < size()-1) {
                jp1pos = m_bp_support.select(v.j + 2);
            }
            v.ipos = ipos;
            v.cipos = m_bp_support.find_close(ipos);
            v.jp1pos = jp1pos;
        }


        //! Computes the suffix number of a leaf node v.
        /*!\param v A valid leaf node of a StreeOhleb.
         * \return The suffix array value corresponding to the leaf node v.
         * \par Time complexity
         *   \f$ \Order{ \saaccess } \f$
         */
        size_type sn(const node_type& v)const
        {
            assert(is_leaf(v));
            return m_csa[v.i];
        }

        //! Computes a unique identification number for a node of the suffx tree in the range [0..nodes()-1]
        /*!
         * \param v A valid node of a StreeOhleb.
         * \return A unique identification number for the node v in the range [0..nodes()-1]
         * \par Time complexity
         *    \f$ \Order{1} \f$
         */
        size_type id(const node_type& v)const
        {
            if (is_leaf(v)) { // return id in the range from 0..csa.size()-1
                return v.i;
            }
            size_type ckpos; // closing parentheses of the l-index
            if (v.cipos > v.jp1pos) { // corresponds to m_lcp[i] <= m_lcp[j+1]
                ckpos     = v.jp1pos-1;
            } else { // corresponds to m_lcp[i] > m_lcp[j+1]
                ckpos    = v.cipos-1;
            }
            assert(m_bp[ckpos]==0);
            size_type r0ckpos = ckpos-m_bp_support.rank(ckpos); // determine the rank of the closing parenthesis
            return size()+m_first_child_rank(r0ckpos);
        }

        //! Computes the node for such that id(v)=id.
        /*!
         * \param id An id in the range [0..nodes()-1].
         * \return A node v of the CST such that id(v)=id.
         * \par Time complexity
         *   \f$ \Order{1} \f$ for leaves and \f$ \Order{\log size()} \f$ for inner nodes
         * \sa id(node_type v)
         */
        node_type inv_id(size_type id)
        {
            if (id < size()) {  // the corresponding node is a leaf
                return select_leaf(id+1);
            } else { // the corresponding node is a inner node
                // (1) get index of the closing parenthesis in m_first_child
                size_type r0ckpos = 0;
                {
                    //binary search for the position of the (id-size()+1)-th set bit in
                    id = id-size()+1;
                    size_type lb = 0, rb = m_bp.size(); // lb inclusive, rb exclusive
                    // invariant: arg(lb) < id, arg(rb) >= id
                    while (rb-lb > 1) {
                        size_type mid = lb + (rb-lb)/2;
                        size_type arg = m_first_child_rank(mid); // ones in the prefix [0..mid-1]
                        if (arg < id) {
                            lb = mid;
                        } else { // arg >= id
                            rb = mid;
                        }
                    }
                    r0ckpos = lb;
                }
                // (2) determine position clpos of the r0clpos-th closing parentheses in the parentheses sequence
                size_type ckpos = 0;
                {
                    // binary search for the position of the (r0ckpos+1)-th closing parenthesis
                    size_type lb = 0, rb = m_bp.size(); // lb inclusive, rb exclusive
                    // invariant: arg(lb) < r0ckpos+1,  arg(rb) >= r0ckpos+1
                    while (rb-lb > 1) {
                        size_type mid = lb + (rb-lb)/2;
                        size_type arg = mid - m_bp_support.rank(mid-1);  // zeros in the prefix [0..mid-1]
                        if (arg < r0ckpos+1) {
                            lb = mid;
                        } else { // arg >= x
                            rb = mid;
                        }
                    }
                    ckpos = lb;
                }
                if (ckpos == m_bp.size()-1) {
                    return root();
                }
                if (m_bp[ckpos+1]) {  // jp1pos < cipos
                    size_type jp1pos= ckpos+1;
                    size_type j     = m_bp_support.rank(jp1pos-1)-1;
                    size_type kpos  = m_bp_support.find_open(ckpos);
                    size_type ipos    = m_bp_support.enclose(kpos);
                    size_type cipos = m_bp_support.find_close(ipos);
                    size_type i        = m_bp_support.rank(ipos-1);
                    return node_type(i, j, ipos, cipos, jp1pos);
                } else { //
                    size_type cipos = ckpos+1;
                    size_type ipos  = m_bp_support.find_open(cipos);
                    size_type i     = m_bp_support.rank(ipos-1);
                    size_type j     = nsv(i, ipos)-1;
                    size_type jp1pos= m_bp.size();
                    if (j != size()-1) {
                        jp1pos = m_bp_support.select(j+2);
                    }
                    return node_type(i, j, ipos, cipos, jp1pos);
                }
            }
        }

        //! Get the number of nodes of the suffix tree.
        size_type nodes()const
        {
            return m_nodes;
        }

        //! Get the node in the suffix tree which corresponds to the lcp-interval [lb..rb]
        /* \param lb Left bound of the lcp-interval [lb..rb] (inclusive).
         * \param rb Right bound of the lcp-interval [lb..rb] (inclusive).
         * \return The node in the suffix tree corresponding lcp-interval [lb..rb]
         * \par Time complexity
         *        \f$ \Order{1} \f$
         */
        node_type node(size_type lb, size_type rb) const
        {
            size_type ipos = m_bp_support.select(lb+1);
            size_type jp1pos;
            if (rb == size()-1) {
                jp1pos = m_bp.size();
            } else {
                jp1pos = m_bp_support.select(rb+2);
            }
            return node_type(lb, rb, ipos, m_bp_support.find_close(ipos), jp1pos);
        }

        //! Maps an index i to the position in TLCP where LCP[i] can be found
        /*!
         * \param i The index in the LCP array
         * \return The corresponding position in the TLCP array
         */
        size_type tlcp_idx(size_type i) const
        {
            size_type ipos     = m_bp_support.select(i+1);
            size_type cipos = m_bp_support.find_close(ipos);
            return m_first_child_rank.rank(((ipos+cipos-1)>>1)-i);
        }
        /* @} */
    };
    

    // == template functions ==


    template<class t_csa, class t_lcp, class t_bp_support, class t_bv, class t_rank, class t_sel>
    StreeOhleb<t_csa, t_lcp, t_bp_support, t_bv, t_rank, t_sel>::StreeOhleb(sdsl::cache_config& config, bool build_only_bps)
    {
        {
            auto event = sdsl::memory_monitor::event("bps-sct");
            sdsl::int_vector_buffer<> lcp_buf(cache_file_name(sdsl::conf::KEY_LCP, config));
            m_nodes = construct_supercartesian_tree_bp_succinct_and_first_child(lcp_buf, m_bp, m_first_child);
            m_nodes += m_bp.size()/2;
            if (m_bp.size() == 2) {  // handle special case, when the tree consists only of the root node
                m_nodes = 1;
            }
        }
        {
            auto event = sdsl::memory_monitor::event("bpss-sct");
            sdsl::util::init_support(m_bp_support, &m_bp);
            sdsl::util::init_support(m_first_child_rank, &m_first_child);
            sdsl::util::init_support(m_first_child_select, &m_first_child);
        }
        if (!build_only_bps) {
            auto event = sdsl::memory_monitor::event("clcp");
            sdsl::cache_config tmp_config(false, config.dir, config.id, config.file_map);
            config.file_map = tmp_config.file_map;
        }
        if (!build_only_bps) {
            auto event = sdsl::memory_monitor::event("load csa");
            load_from_cache(m_csa,std::string(sdsl::conf::KEY_CSA)+"_"+sdsl::util::class_to_hash(m_csa), config);
        }
    }

    template<class t_csa, class t_lcp, class t_bp_support, class t_bv, class t_rank, class t_sel>
    auto StreeOhleb<t_csa, t_lcp, t_bp_support, t_bv, t_rank, t_sel>::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const -> size_type
    {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += m_csa.serialize(out, child, "csa");
        written_bytes += m_bp.serialize(out, child, "bp");
        written_bytes += m_bp_support.serialize(out, child, "bp_support");
        written_bytes += m_first_child.serialize(out, child, "mark_child");
        written_bytes += m_first_child_rank.serialize(out, child, "mark_child_rank");
        written_bytes += m_first_child_select.serialize(out, child, "mark_child_select");
        written_bytes += write_member(m_nodes, out, child, "node_cnt");
        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    template<class t_csa, class t_lcp, class t_bp_support, class t_bv, class t_rank, class t_sel>
    void StreeOhleb<t_csa, t_lcp, t_bp_support, t_bv, t_rank, t_sel>::load(std::istream& in)
    {
        m_csa.load(in);
        m_bp.load(in);
        m_bp_support.load(in, &m_bp);
        m_first_child.load(in);
        m_first_child_rank.load(in,&m_first_child);
        m_first_child_select.load(in,&m_first_child);
        sdsl::read_member(m_nodes, in);
    }

    template<class t_csa, class t_lcp, class t_bp_support, class t_bv, class t_rank, class t_sel>
    StreeOhleb<t_csa, t_lcp, t_bp_support, t_bv, t_rank, t_sel>& StreeOhleb<t_csa, t_lcp, t_bp_support, t_bv, t_rank, t_sel>::operator=(const StreeOhleb& cst)
    {
        if (this != &cst) {
            copy(cst);
        }
        return *this;
    }

    template<class t_csa, class t_lcp, class t_bp_support, class t_bv, class t_rank, class t_sel>
    StreeOhleb<t_csa, t_lcp, t_bp_support, t_bv, t_rank, t_sel>& StreeOhleb<t_csa, t_lcp, t_bp_support, t_bv, t_rank, t_sel>::operator=(StreeOhleb&& cst)
    {
        if (this != &cst) {
            m_csa              = std::move(cst.m_csa);
            m_bp               = std::move(cst.m_bp);
            m_bp_support       = std::move(cst.m_bp_support);
            m_bp_support.set_vector(&m_bp);
            m_first_child      = std::move(cst.m_first_child);
            m_first_child_rank = std::move(cst.m_first_child_rank);
            m_first_child_rank.set_vector(&m_first_child);
            m_first_child_select = std::move(cst.m_first_child_select);
            m_first_child_select.set_vector(&m_first_child);
            m_nodes            = std::move(cst.m_nodes);
        }
        return *this;
    }

    StreeOhleb<>::size_type load_or_build(StreeOhleb<>& st, const std::string &s, const std::string potential_stree_fname, const bool load){
        using timer = std::chrono::high_resolution_clock;
        auto start = timer::now();
        if(load){
            std::cerr << " * loading the CST from " << potential_stree_fname << " ";
            sdsl::load_from_file(st, potential_stree_fname);
        } else {
            std::cerr << " * building the CST of length " << s.size() << " ";
            sdsl::construct_im(st, s, 1);
        }
        auto stop = timer::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    }


    template<class t_int>
    struct bp_interval {
        t_int i;     //!< The left border of the lcp-interval \f$\ell-[left..right]\f$.
        t_int j;     //!< The right border of the lcp-interval \f$\ell-[left..right]\f$.
        t_int ipos;  // position of the i+1th opening parenthesis in the balanced parentheses sequence
        t_int cipos; // position of the matching closing parenthesis of the i+1th opening parenthesis in the balanced parentheses sequence
        t_int jp1pos;// position of the j+2th opening parenthesis in the balanced parentheses sequence

        //! Constructor
        bp_interval(t_int i=0, t_int j=0, t_int ipos=0, t_int cipos=0, t_int jp1pos=0):i(i),j(j),ipos(ipos),cipos(cipos),jp1pos(jp1pos) {};

        //! Copy constructor
        bp_interval(const bp_interval& iv) = default;
        //! Move copy constructor
        bp_interval(bp_interval&& iv) = default;

        bool operator<(const bp_interval& interval)const
        {
            if (i!=interval.i)
                return i<interval.i;
            return j<interval.j;
        }
        
        //! Equality operator.
        /*! Two lcp-intervals are equal if and only if all their corresponding member variables have the same values.
         */
        bool operator==(const bp_interval& interval)const
        {
            return i==interval.i and j==interval.j;
        }
        
        //! Inequality operator.
        /*! Two lcp-intervals are not equal if and only if not all their corresponding member variables have the same values.
         */
        bool operator!=(const bp_interval& interval)const
        {
            return !(*this==interval);
        }
        
        //! Assignment operator.
        bp_interval& operator=(const bp_interval& interval) = default;
        //! Move assignment
        bp_interval& operator=(bp_interval&& interval) = default;
    };
    
    
} // end namespace sdsl
#endif
