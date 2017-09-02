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

public:
    sdsl::wt_huff<> m_wt;
    
    size_type        m_size  = 0;    // original text size
    size_type        m_sigma = 0;    // alphabet size
    bit_vector_type  m_bv;           // bit vector to store the wavelet tree
    rank_1_type      m_bv_rank;      // rank support for the wavelet tree bit vector
    select_1_type    m_bv_select1;   // select support for the wavelet tree bit vector
    select_0_type    m_bv_select0;
    tree_strat_type  m_tree;

    
    mywt_pc(sdsl::wt_huff<>& wt) {
        m_wt = wt;
        m_size = wt.size();
        m_sigma = wt.get_sigma();
        m_bv_rank = wt.get_bv_rank();
        m_tree = wt.get_tree();
        m_bv_select0 = wt.get_bv_select0();
        m_bv_select1 = wt.get_bv_select1();

    }

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
        return bit_select0(v, bit_rank0(v, i + 1) + cnt);
    }
    
    inline size_type bit_select_at_dist1(const node_type v, const size_type i, const size_type cnt) const {
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
        
#define SAME_PREFIX_PATH(p1, p2) (((p1) << 1) == ((p2) << 1))
#define PATH_LENGTH(p) ((p) >> 56)
#define REVERSED_PATH(p) ( (p) << (64 - (PATH_LENGTH(p))) )

        uint64_t p = m_tree.bit_path(c), pt_i = m_tree.bit_path(m_wt[i]);
        uint32_t p_len = (p >> 56), pt_i_len = (pt_i >> 56);
        node_type v = m_tree.root();
        std::vector<size_type> i_vec(p_len); // place the i-values here
        std::vector<size_type> j_vec(p_len); // place the j-values here -- TODO: only need 2 values actually
        i_vec[0] = i;

        for(uint32_t i = 1; i < p_len; i++, p >>= 1) {
            if (p&1)
                i_vec[i] = bit_rank1(v, i_vec[i - 1]);
            else
                i_vec[i] = bit_rank0(v, i_vec[i - 1]);
            cout << "i" << i << " = r" << (p&1) << "(" << i_vec[i - 1] << ") = " << i_vec[i] << endl;
            v = m_tree.child(v, p&1); // goto child
        }
        cout << "i_vec: [";
        for(int aa = 0; aa < p_len - 1; aa++)
            cout << i_vec[aa] << ", ";
        cout <<  i_vec[p_len - 1] << "]" << endl;

        // reset path
        p = m_tree.bit_path(c);
        p <<= (64- (p >> 56));
        pt_i <<= (64 - (pt_i >> 56) - (p_len > pt_i_len ? p_len - pt_i_len : 0));

        size_type jk = 0, j_prev = 0;
        //if(p != pt_i)
        if(!SAME_PREFIX_PATH(p, pt_i))
            i_vec[p_len - 1] -= 1;
        if((p & 0x8000000000000000ULL) == 0)
            j_prev = bit_select_at_dist0(v, i_vec[p_len - 1], cnt);
        else
            j_prev = bit_select_at_dist1(v, i_vec[p_len - 1], cnt);
        cout << (p != pt_i ? "*" : "") << "j" << p_len << " = sd" << ((p & 0x8000000000000000ULL)!=0) << "(" << i_vec[p_len - 1] << ", " << cnt << ") = " << j_prev << endl;
        v  = m_tree.parent(v);
        p <<= 1;
        pt_i <<= 1;
        
        for(uint32_t idx = p_len - 1; idx > 0; idx--, p <<= 1, pt_i <<= 1){
            //if(p != pt_i)
            if(!SAME_PREFIX_PATH(p, pt_i))
                i_vec[idx - 1] -= 1;
            if ((p & 0x8000000000000000ULL)==0)
                jk = bit_select_at_dist0(v, i_vec[idx - 1], j_prev - i_vec[idx]);
            else
                jk = bit_select_at_dist1(v, i_vec[idx - 1], j_prev - i_vec[idx]);
            cout << (p != pt_i ? "*" : "") << "j" << idx << " = sd" << ((p & 0x8000000000000000ULL)!=0) << "(" << i_vec[idx - 1] << ", " << cnt << ") = " << jk << endl;
            v   = m_tree.parent(v);
            j_prev = jk;
        }
        assert(jk == m_wt.select(m_wt.rank(i + 1, c) + cnt, c));
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
    /*
    cout << "('n', 8, 1)" << endl << wwt.select_at_dist('n', 8, 1) << endl << endl;
    cout << "('c', 6, 1)" << endl << wwt.select_at_dist('c', 6, 1) << endl << endl;
    cout << "('a', 8, 1)" << endl << wwt.select_at_dist('a', 8, 1) << endl << endl;
    cout << "('a', 11, 1)" << endl << wwt.select_at_dist('a', 11, 1) << endl << endl;
    */
    //cout << "('b', 16, 0)" << endl << wwt.select_at_dist('b', 16, 0) << endl << endl;
    char c = 'd';
    for(uint k = 0; k < 8; k++){
        for(uint i = 0; i < 14 - k; i++){
            if(k + i == 0)
                continue;
            cout << "('" << c << "', "<< i << ", " << k << ")" << endl << wwt.select_at_dist(c, i, k) << endl << endl;
        }
    }
}

int main(int argc, const char * argv[]) {
    ab();
    return 0;
}
