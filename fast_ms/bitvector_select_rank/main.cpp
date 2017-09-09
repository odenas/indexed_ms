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
#include <sdsl/bits.hpp>

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
        m_bv = wt.get_bv();
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
    size_type select_at_dist(const sdsl::wt_huff<>::value_type c, const size_type i, const size_type cnt) const
    {
        //for(int i=0; i < m_bv.size(); i++)
        //    cout << "m_bv[" << i << "] = " << m_bv[i] << endl;

        uint64_t p = m_tree.bit_path(c), pt_i = m_tree.bit_path(m_wt[i]);
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

        size_type jk = 0, j_prev = 0;
        if(!equal_prefix[p_len - 1] && i_vec[p_len - 1] > 0)
            i_vec[p_len - 1] -= 1;
        if((p & 0x8000000000000000ULL) == 0)
            j_prev = bit_select_at_dist0(v, i_vec[p_len - 1], cnt);
        else
            j_prev = bit_select_at_dist1(v, i_vec[p_len - 1], cnt);
        cout << (equal_prefix[p_len - 1] ? "" : "*") << "j" << p_len << " = sd" << ((p & 0x8000000000000000ULL)!=0) << "(" << i_vec[p_len - 1] << ", " << cnt << ") = " << j_prev << endl;
        v  = m_tree.parent(v);
        p <<= 1;
        
        for(uint32_t idx = p_len - 1; idx > 0; idx--, p <<= 1){
            if(!equal_prefix[idx - 1] && i_vec[idx - 1] > 0)
                i_vec[idx - 1] -= 1;
            if ((p & 0x8000000000000000ULL)==0)
                jk = bit_select_at_dist0(v, i_vec[idx - 1], j_prev - i_vec[idx]);
            else
                jk = bit_select_at_dist1(v, i_vec[idx - 1], j_prev - i_vec[idx]);
            cout << (equal_prefix[idx - 1] ? "" : "*") << "j" << idx << " = sd" << ((p & 0x8000000000000000ULL)!=0) << "(" << i_vec[idx - 1] << ", " << cnt << ") = " << jk << endl;
            v = m_tree.parent(v);
            j_prev = jk;
        }
        assert(j_prev == m_wt.select(m_wt.rank(i + 1, c) + cnt, c));
        return j_prev;
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
    cout << "('c', 1, 1)" << endl << wwt.select_at_dist('\x0', 1, 1) << endl << endl;
    /*
    char c = 'a';
    for(uint i = 1; i < 20; i++){
        for(uint k = 1; k < 50; k++){
            if(k + i == 0)
                continue;
            cout << "('" << c << "', "<< i << ", " << k << ")" << endl << wwt.select_at_dist(c, i, k) << endl << endl;
        }
    }
    c = 'b';
    for(uint i = 1; i < 20; i++){
        for(uint k = 1; k < 35; k++){
            if(k + i == 0)
                continue;
            cout << "('" << c << "', "<< i << ", " << k << ")" << endl << wwt.select_at_dist(c, i, k) << endl << endl;
        }
    }
    */
}

int main(int argc, const char * argv[]) {
    ab();
    return 0;
}
