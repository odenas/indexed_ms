//
//  maxrep_vector.hpp
//  fast_ms
//
//  Created by denas on 10/21/17.
//  Copyright © 2017 denas. All rights reserved.
//

#ifndef maxrep_vector_h
#define maxrep_vector_h


#include <iostream>
#include <string>
#include <vector>

#include <sdsl/int_vector.hpp>
#include "stree_sct3.hpp"

namespace fdms {
    using timer = std::chrono::high_resolution_clock;

    template<typename cst_t, typename vector_t>
    class Maxrep{
        //typedef StreeOhleb<> cst_t;
        typedef typename cst_t::node_type node_type;
        typedef typename cst_t::size_type size_type;
        //typedef sdsl::bit_vector vector_t;
        
    private:
        // strategies for determining maximality
        typedef bool (*maximal_test) (const cst_t& m_cst, const node_type& v) ;
        
        vector_t m_vec;

        void build_maxrep(const cst_t& m_st, maximal_test maximal_test_f_ptr){
            node_type currnode = m_st.root(), nextnode = m_st.root();
            bool direction_down = true, currnode_maximal = false;
            
            do{
                if(direction_down){
                    if(!m_st.is_leaf(currnode)){ // process currnode
                        if(maximal_test_f_ptr(m_st, currnode)){
                            m_vec[currnode.i] = m_vec[currnode.j] = currnode_maximal = 1;
                            
                            // try going down the subtree
                            nextnode = m_st.first_child(currnode);
                            if(m_st.is_root(nextnode) or !currnode_maximal)
                                direction_down = false;
                            else
                                currnode = nextnode;
                        } else { // no node in subtree is maximal
                            direction_down = false;
                        }
                    } else {
                        direction_down = false;
                    }
                } else {
                    nextnode = m_st.sibling(currnode);
                    if(m_st.is_root(nextnode)) {
                        currnode = m_st.parent(currnode);
                    } else {
                        currnode = nextnode;
                        direction_down = true;
                    }
                }
            } while(!m_st.is_root(currnode));
        }

        static void load_vec(const std::string fname, vector_t& vec){
            sdsl::load_from_file(vec, fname);
        }
        
        static void load_cst(const std::string fname, cst_t& st){
            sdsl::load_from_file(st, fname);
        }
        
    public:
        Maxrep(const cst_t& st, const bool verbose = false) {
            m_vec = vector_t(st.size()  + 1);
            sdsl::util::set_to_value(m_vec, 0);
            maximal_test mtest{Maxrep::rank_maximal_test};

            if(verbose){
                auto start = timer::now();
                std::cerr << " * computing MAXREP ";
                build_maxrep(st, mtest);
                int t = (int) std::chrono::duration_cast<std::chrono::milliseconds>(timer::now() - start).count();
                std::cerr << "DONE (" << t << " milliseconds)" << std::endl;
            } else {
                build_maxrep(st, mtest);
            }
        }

        Maxrep(){}

        Maxrep(const vector_t& v) {m_vec = v;}

        Maxrep(const Maxrep& m) {m_vec = m.m_vec;}
        
        static bool wl_cnt_maximal_test(const cst_t& m_st, const node_type& v) {
            size_type wl_count = 0;
            for(uint16_t i=0; i<m_st.csa.sigma; i++)
                wl_count += (m_st.is_root(m_st.single_rank_wl(v, m_st.csa.comp2char[i])) ? 0 : 1);
            return wl_count > 1;
        }
        
        static bool rank_maximal_test(const cst_t& st, const node_type& v) {
            char c = st.csa.bwt[v.j];
            auto ni = st.csa.bwt.double_rank(v.i, v.j + 1, c);
            size_type count = (ni.second - ni.first);
            return count != (v.j - v.i + 1);
        }

        bool operator[](size_type i) const {
            return m_vec[i];
        }

        inline bool is_maximal(const node_type v) const{
            return (v.i != v.j && m_vec[v.i] == 1 && m_vec[v.j] == 1);
        }
        
        inline bool is_intnode_maximal(const node_type v) const {
            return (m_vec[v.i] == 1 && m_vec[v.j] == 1);
        }
        
        inline size_type size() const { return m_vec.size(); }
        
        inline string desc(node_type v) const {
            return "(" + string(std::to_string(m_vec[v.i]) + ", " + std::to_string(m_vec[v.j]) + ")");
        }
        
        void set_to_one(const size_type size){
            m_vec.resize(size);
            sdsl::util::set_to_value(m_vec, 1);
        }

        void dump_vec(const std::string fname){
            sdsl::store_to_file(m_vec, fname);
        }

        static vector_t load_vec(const std::string fname){
            vector_t vec;
            sdsl::load_from_file(vec, fname);
            return vec;
        }

        static Maxrep load(const std::string vec_fname){
            return Maxrep(load_vec(vec_fname));
        }
        
        static size_type load_or_build(Maxrep& mr, const cst_t& st, const string& load_fname, const bool load){
            auto start = timer::now();
            if(load){
                cerr << " * loading MAXREP from " << load_fname << " ";
                vector_t vec;
                sdsl::load_from_file(vec, load_fname);
                mr = Maxrep(vec);
            } else {
                cerr << " * building MAXREP of length " << st.size() << " ";
                mr = Maxrep(st, false);
            }
            auto stop = timer::now();
            return std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
        }
    };
}


#endif /* maxrep_vector_h */
