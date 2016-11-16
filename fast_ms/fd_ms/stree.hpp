//
//  stree.h
//  fast_ms
//
//  Created by denas on 11/12/16.
//  Copyright © 2016 denas. All rights reserved.
//

#ifndef stree_h
#define stree_h

#include "bp.hpp"
#include "fd_ms.hpp"

using namespace std;


namespace fdms {

typedef uint8_t       char_type;
typedef unsigned long size_type;
typedef size_type     node_type;


/*
// Gets ISA[SA[idx]+d]
// d = depth of the character 0 = first position
size_type get_char_pos(size_type idx, size_type d, const t_csa& csa) {
    if (d == 0)
        return idx;
    // if we have to apply \f$\LF\f$ or \f$\Phi\f$ more
    // than 2*d times to calc csa(csa[idx]+d), we opt to
    // apply \f$ \Phi \f$ d times
    if (csa.sa_sample_dens + csa.isa_sample_dens > 2*d+2) {
        for (size_type i=0; i < d; ++i)
            idx = csa.psi[idx];
        return idx;
    }
    return csa.isa[csa[idx] + d];
}
*/
class Stree{
private:
    fdms::bp_support_sada<>&       m_bp_supp;
    sdsl::rank_support_v5<10,2>    m_bp_rank10;
    sdsl::select_support_mcl<10,2> m_bp_select10;
    Bwt& m_bwt;

public:
    Stree(fdms::bp_support_sada<>& bp_supp, Bwt& bwt) : m_bwt{bwt}, m_bp_supp{bp_supp}
    {
        sdsl::rank_support_v5<10,2> m_bp_rank10 (m_bp_supp.m_bp);
        sdsl::select_support_mcl<10,2> m_bp_select10(m_bp_supp.m_bp);
    }

    node_type root() const { return 0; }

    bool is_leaf(node_type v) const {
        assert((*m_bp_supp.m_bp)[v] == 1);  // assert v is a valid node
        // if there is a closing parenthesis at position v+1, the node is a leaf
        return !(*m_bp_supp.m_bp)[v + 1];
    }

    size_type parent(node_type v) const {
        assert((*m_bp_supp.m_bp)[v] == 1); //assert valid node
        if(v == root())
            return root();
        return m_bp_supp.enclose(v);
    }

    node_type child(node_type v, const char_type c, size_type& char_pos)const
    {
        if (is_leaf(v))  // if v is a leaf = (), v has no child
            return root();
        // else v = ( (     ))
        size_type cc = m_bwt.char2int[c];
        if (cc == 0 and c != 0) // TODO: aendere char2comp so ab, dass man diesen sonderfall nicht braucht
            return root();
        size_type char_ex_max_pos = m_bwt.C[cc+1], char_inc_min_pos = m_bwt.C[cc];

        //size_type d = depth(v);  // time complexity: \lcpaccess
        size_type res = v+1;
        while (true) {
            //if (is_leaf(res)) {
            //    char_pos = get_char_pos(m_bp_rank10(res), d, m_csa);
            //} else {
            //    char_pos = get_char_pos(inorder(res), d, m_csa);
            //}
            if (char_pos >= char_ex_max_pos)  // if the current char is lex. greater than the searched char: exit
                return root();
            if (char_pos >= char_inc_min_pos)  // if the current char is lex. equal with the
                return res;
            res = m_bp_supp.find_close(res) + 1;
            if (!(*m_bp_supp.m_bp)[res]) // closing parenthesis: there exists no next child
                return root();
        }
    }
};

}

#endif /* stree_h */
