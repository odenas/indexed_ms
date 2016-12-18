//
//  stree.h
//  fast_ms
//
//  Created by denas on 11/12/16.
//  Copyright © 2016 denas. All rights reserved.
//

#ifndef stree_sada_h
#define stree_sada_h

#include "bp_geary.hpp"
#include "bp_sada.hpp"

#include "Bwt.hpp"

using namespace std;


namespace fdms {

    template<class t_bp_support>
    class StreeSada{
    public:
        typedef t_bp_support bp_support_type;
        typedef size_type    node_type;

    private:
        bp_support_type&       m_bp_supp;
        Bwt& m_bwt;

    public:

        sdsl::rank_support_v5<10,2>    m_bp_rank10;
        sdsl::select_support_mcl<10,2> m_bp_select10;
        size_type size_in_bytes__rank, size_in_bytes__select;

        StreeSada(bp_support_type& bp_supp, Bwt& bwt) : m_bwt{bwt}, m_bp_supp{bp_supp} {
            sdsl::util::init_support(m_bp_rank10, m_bp_supp.m_bp);
            sdsl::util::init_support(m_bp_select10, m_bp_supp.m_bp);
            size_in_bytes__rank = sdsl::size_in_bytes(m_bp_rank10); //10 rank support
            size_in_bytes__select = sdsl::size_in_bytes(m_bp_select10); //10 select support
        }

        node_type root() const { return 0; }

        bool is_leaf(node_type v) const {
            assert((*m_bp_supp.m_bp)[v] == 1);  // assert v is a valid node
            // if there is a closing parenthesis at position v+1, the node is a leaf
            return !(*m_bp_supp.m_bp)[v + 1];
        }

        node_type parent(node_type v) const {
            assert((*m_bp_supp.m_bp)[v] == 1); //assert valid node
            if(v == root())
                return root();
            return m_bp_supp.enclose(v);
        }

        //!Calculates the index of the leftmost leaf in the corresponding suffix array.
        /*!\param v A valid node of the suffix tree.
         * \return The index of the leftmost leaf in the corresponding suffix array.
         * \par Time complexity
         *  \f$ \Order{1} \f$
         * \par Note
         * lb is an abbreviation for ,,left bound''
         */
        size_type lb(const node_type v) const { return m_bp_rank10(v); }

        //! Calculates the index of the rightmost leaf in the corresponding suffix array.
        /*! \param v A valid node of the suffix tree.
         *  \return The index of the rightmost leaf in the corresponding suffix array.
         *  \par Time complexity
         *    \f$ \Order{1} \f$
         *  \par Note
         *   rb is an abbreviation for ,,right bound''
         */
        size_type rb(const node_type v) const {
            size_type r = m_bp_supp.find_close(v);
            return m_bp_rank10(r + 1) - 1;
        }

        node_type select_leaf(size_type i) const {
            assert(i > 0 and i <= m_bwt.bwt_len);
            // -1 as select(i) returns the postion of the 0 of pattern 10
            return m_bp_select10.select(i) - 1;
        }

        //! Calculate the lowest common ancestor (lca) of two nodes v and w of the suffix tree.
        /*!
         * \param v The first node for which the lca with the second node should be computed.
         * \param w The second node for which the lca with the first node should be computed.
         * \return A node that is the lowest common ancestor of v and w in the suffix tree.
         * \par Time complexity
         *    \f$ \Order{\rrenclose}\   \f$
         */
        node_type lca(node_type v, node_type w) const {
            assert((*m_bp_supp.m_bp)[v] == 1 and (*m_bp_supp.m_bp)[w] == 1);
            if (v > w)
                std::swap(v, w);
            else if (v == w)
                return v;

            if (v == root())
                return root();
            return m_bp_supp.double_enclose(v, w);
        }

        //! Compute the Weiner link of node v and character c.
        /*
         * \param v A valid node of a cst_sada.
         * \param c The character which should be prepended to the string of the current node.
         *   \return root() if the Weiner link of (v, c) does not exist, otherwise the Weiner link is returned.
         * \par Time complexity
         *    \f$ \Order{ t_{rank\_bwt} + t_{lca}}\f$
         */
        node_type wl(node_type v, const char_type c) const {
            // get leftmost leaf in the tree rooted at v
            size_type left        = m_bp_rank10(v);
            // get the rightmost leaf in the tree rooted at v
            size_type right = is_leaf(v) ? left : m_bp_rank10(m_bp_supp.find_close(v)) - 1;

            size_type c_left    = m_bwt.rank(left, c);
            size_type c_right   = m_bwt.rank(right + 1, c);

            if (c_left == c_right)  // there exists no Weiner link
                return root();

            if (c_left + 1 == c_right)
                return select_leaf(m_bwt.C[m_bwt.char2int[c]] + c_left + 1);
            else {
                size_type left    = m_bwt.C[m_bwt.char2int[c]] + c_left;
                size_type right    = m_bwt.C[m_bwt.char2int[c]] + c_right - 1;
                assert(left < right);
                node_type left_leaf = select_leaf(left+1);
                node_type right_leaf= select_leaf(right+1);
                return lca(left_leaf, right_leaf);
            }
        }

        node_type lazy_wl(node_type v, const char_type c) const { return wl(v, c); }

    };

    void stree_test(){
        string Sfwd{"aabbaba"};
        string Sfwdbp{"11011010110100011101001000"};

        bvector bfwd(Sfwdbp.size());
        for(size_type i=0; i<Sfwdbp.size(); i++)
            bfwd[i] = ((unsigned char)Sfwdbp[i] - 48);

        Bwt Bwtfwd(Sfwd);
        fdms::bp_support_sada<> Bpsfwd(&bfwd);
        StreeSada<fdms::bp_support_sada<>> st(Bpsfwd, Bwtfwd);
        
        //{"11011010110100011101001000"};
        //  01 34 6 89 1   567 9  2
        //            1         2
        
        size_type idx[] {0, 1, 3, 4, 6, 8, 9, 11, 15, 16, 17, 19, 22};
        for(size_type i = 0; i < 12; i++){
            cout << "[" << idx[i] << "]" << endl;
        }
        
        for(size_type i = 0; i < 12; i++){
            cout << "[" << idx[i] << "]" << endl;
            cout << "'a'" << st.wl(idx[i], 'a');
            cout << "'b'" << st.wl(idx[i], 'b');
            cout << "'#'" << st.wl(idx[i], '#') << endl;
        }
    }

}

#endif /* stree_h */
