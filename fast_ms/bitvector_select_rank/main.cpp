//
//  main.cpp
//  bitvector_select_rank
//
//  Created by denas on 8/30/17.
//  Copyright © 2017 denas. All rights reserved.
//

#include <iostream>
#include <sdsl/int_vector.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/select_support.hpp>
#include <sdsl/util.hpp>

using namespace std;

class mywt_pc
{
    typedef sdsl::wt_huff<>::node_type node_type;
    typedef sdsl::wt_huff<>::size_type size_type;
    typedef sdsl::wt_huff<>::bit_vector_type bit_vector_type;
    typedef sdsl::wt_huff<>::rank_1_type rank_1_type;
    typedef sdsl::wt_huff<>::select_1_type select_1_type;
    typedef sdsl::wt_huff<>::select_0_type select_0_type;
    typedef sdsl::wt_huff<>::tree_strat_type tree_strat_type;

private:
    sdsl::wt_huff<> m_wt;
    
    size_type        m_size  = 0;    // original text size
    size_type        m_sigma = 0;    // alphabet size
    bit_vector_type  m_bv;           // bit vector to store the wavelet tree
    rank_1_type      m_bv_rank;      // rank support for the wavelet tree bit vector
    select_1_type    m_bv_select1;   // select support for the wavelet tree bit vector
    select_0_type    m_bv_select0;
    tree_strat_type  m_tree;


public:
    mywt_pc(sdsl::wt_huff<>& wt) {
        m_wt = wt;
        m_size = wt.size();
        m_sigma = wt.get_sigma();
        m_bv_rank = wt.get_bv_rank();
        m_tree = wt.get_tree();
        m_bv_select0 = wt.get_bv_select0();
        m_bv_select1 = wt.get_bv_select1();

    }

    size_type bit_rank1(const node_type v, const size_type i) const {
        return m_bv_rank(m_tree.bv_pos(v) + i) - m_tree.bv_pos_rank(v);
    }
    
    size_type bit_rank0(const node_type v, const size_type i) const {
        return i - bit_rank1(v, i);
    }
    
    size_type bit_select1(const node_type v, const size_type i) const {
        return m_bv_select1(m_tree.bv_pos_rank(v) + i + 1) - m_tree.bv_pos(v);
    }
    
    size_type bit_select0(const node_type v, const size_type i) const {
        return m_bv_select0(m_tree.bv_pos(v) - m_tree.bv_pos_rank(v) + i + 1) - m_tree.bv_pos(v);
    }
    
    size_type bit_select_at_dist0(const node_type v, const size_type i, const size_type cnt) const {
        return bit_select0(v, bit_rank0(v, i + 1) + cnt);
    }
    
    size_type bit_select_at_dist1(const node_type v, const size_type i, const size_type cnt) const {
        return bit_select1(v, bit_rank1(v, i + 1) + cnt);
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
    size_type select_at_dist(const sdsl::wt_huff<>::value_type c, const size_type i, const size_type cnt) const
    {
        uint64_t p = m_tree.bit_path(c);
        uint32_t path_len = (p>>56);
        node_type v = m_tree.root();
        std::vector<size_type> i_vec(path_len); // place the i-values here
        std::vector<size_type> j_vec(path_len); // place the j-values here -- TODO: only need 2 values actually
        i_vec[0] = i;
        

        for(uint32_t i = 1; i < path_len; i++, p >>= 1) {
            if (p&1)
                i_vec[i] = bit_rank1(v, i_vec[i - 1]);
            else
                i_vec[i] = bit_rank0(v, i_vec[i - 1]);
            v = m_tree.child(v, p&1); // goto child
        }
        
        // reset path
        p = m_tree.bit_path(c);
        p <<= (64-path_len);
        
        size_type jk = 0, j_prev = 0;
        p <<= 1;
        if((p & 0x8000000000000000ULL) == 0)
            j_prev = bit_select_at_dist0(v, i_vec[path_len - 1], cnt);
        else
            j_prev = bit_select_at_dist1(v, i_vec[path_len - 1], cnt);
        v  = m_tree.parent(v);
        
        for(uint32_t i = path_len - 1; i > 0; i--, p <<= 1){
            if ((p & 0x8000000000000000ULL)==0)
                jk = bit_select_at_dist0(v, i_vec[i - 1], j_prev - i_vec[i]);
            else
                jk = bit_select_at_dist1(v, i_vec[i - 1], j_prev - i_vec[i]);
            v   = m_tree.parent(v);
            j_prev = jk;
        }
        return jk;
    }
};

void aa(){
    sdsl::bit_vector v(10);
    v[0] = 1;
    v[1] = 0;
    v[2] = 0;
    v[3] = 0;
    v[4] = 1;
    v[5] = 1;
    v[6] = 0;
    v[7] = 0;
    v[8] = 1;
    v[9] = 0;

    sdsl::bit_vector::select_1_type bv_select;
    sdsl::bit_vector::select_0_type bv_select0;
    sdsl::util::init_support(bv_select, &v);
    sdsl::util::init_support(bv_select0, &v);
    
    
    cout << "index of 2nd 1 = " << bv_select(2) << endl;
}

void ab(){
    sdsl::wt_huff<> wt;
    sdsl::construct(wt, "/Users/denas/projects/matching_statistics/indexed_ms/tests/wt_txt", 1);
    mywt_pc wwt(wt);
    
    //cout << wt.rank(6, 'c') << endl;
    //cout << wt.rank(7, 'c') << endl;
    cout << wwt.select_at_dist('n', 8, 1) << endl;
}

int main(int argc, const char * argv[]) {
    ab();
    return 0;
}
