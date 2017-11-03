//
//  parent_depth_list.hpp
//  fast_ms
//
//  Created by denas on 11/1/17.
//  Copyright © 2017 denas. All rights reserved.
//

#ifndef parent_depth_list_h
#define parent_depth_list_h

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>


#include "stree_sct3.hpp"
//#include "maxrep_vector.hpp"

using namespace std;


namespace fdms {
    typedef StreeOhleb<> cst_t;
    typedef cst_t::node_type node_type;
    typedef cst_t::size_type size_type;
    typedef cst_t::char_type char_type;
    typedef uint8_t char_value;
    using timer = std::chrono::high_resolution_clock;
    
    /*
     basically a pair (node, # parents, c)
     
     c parent calls to node will give the lowest ancestor with a weiner link with c
     */
    class node_with_depth{
    public:
        node_type m_node;
        size_type m_depth;
        char_type m_c;
        
        node_with_depth(){
            m_node = {0, 0, 0, 0, 0};
            m_depth = 0;
            m_c = '\0';
        }
        
        node_with_depth(const node_with_depth& other) : m_node{other.m_node}, m_depth{other.m_depth}, m_c{other.m_c} {}

        node_with_depth(const node_type n, const size_type d, const char_type c) : m_node{n}, m_depth{d}, m_c{c} {}

        friend std::ostream& operator<<(std::ostream& out, const node_with_depth& nwd){
            return out << std::setw(10) << nwd.m_node.i << ' '
            << std::setw(10) << nwd.m_node.j << ' '
            << std::setw(10) << nwd.m_node.ipos << ' '
            << std::setw(10) << nwd.m_node.cipos << ' '
            << std::setw(10) << nwd.m_node.jp1pos << ' '
            << std::setw(10) << nwd.m_c << ' '
            << std::setw(10) << nwd.m_depth << '\n';
        }
        friend std::istream& operator>>(std::istream& in, node_with_depth& nwd){
            return in >> nwd.m_node.i >>
            nwd.m_node.j >>
            nwd.m_node.ipos >>
            nwd.m_node.cipos >>
            nwd.m_node.jp1pos >>
            nwd.m_c >>
            nwd.m_depth;
        }
        
        friend std::ifstream& read(std::ifstream& in, node_with_depth& e){
            in.read(reinterpret_cast<char*>(&e.m_node.i), sizeof(e.m_node.i));
            in.read(reinterpret_cast<char*>(&e.m_node.j), sizeof(e.m_node.j));
            in.read(reinterpret_cast<char*>(&e.m_node.ipos), sizeof(e.m_node.ipos));
            in.read(reinterpret_cast<char*>(&e.m_node.cipos), sizeof(e.m_node.cipos));
            in.read(reinterpret_cast<char*>(&e.m_node.jp1pos), sizeof(e.m_node.jp1pos));
            in.read(reinterpret_cast<char*>(&e.m_c), sizeof(e.m_c));
            in.read(reinterpret_cast<char*>(&e.m_depth), sizeof(e.m_depth));
            return in;
        };
        
        friend std::ofstream& write(std::ofstream& out, node_with_depth& e){
            out.write(reinterpret_cast<char*>(&e.m_node.i), sizeof(e.m_node.i));
            out.write(reinterpret_cast<char*>(&e.m_node.j), sizeof(e.m_node.j));
            out.write(reinterpret_cast<char*>(&e.m_node.ipos), sizeof(e.m_node.ipos));
            out.write(reinterpret_cast<char*>(&e.m_node.cipos), sizeof(e.m_node.cipos));
            out.write(reinterpret_cast<char*>(&e.m_node.jp1pos), sizeof(e.m_node.jp1pos));
            out.write(reinterpret_cast<char*>(&e.m_c), sizeof(e.m_c));
            out.write(reinterpret_cast<char*>(&e.m_depth), sizeof(e.m_depth));
            return out;
        };
        
        bool operator==(const node_with_depth x){
            return (m_node.i == x.m_node.i &&
                    m_node.j == x.m_node.j &&
                    m_node.ipos == x.m_node.ipos &&
                    m_node.cipos == x.m_node.cipos &&
                    m_node.jp1pos == x.m_node.jp1pos &&
                    m_c == x.m_c &&
                    m_depth == x.m_depth);
        }
        
        string repr() const {
            return string("<(" +
                          std::to_string(m_node.i) + "," +
                          std::to_string(m_node.j) + "," +
                          std::to_string(m_node.ipos) + "," +
                          std::to_string(m_node.cipos) + "," +
                          std::to_string(m_node.jp1pos) + ") " +
                          std::to_string(m_c) + ", " +
                          std::to_string(m_depth) + ">");
        }

    };

    class NwdList {
    private:
        void report_progress(timer::time_point start_time, size_type curr_idx, size_type total){
            timer::time_point stop_time = timer::now();
            size_type elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time).count() + 1;
            cerr << endl << "[" << elapsed / 1000 << " s] " << 100.0 * curr_idx / total << "% @ " << (1.0 * curr_idx / elapsed) << " KHz";
        }

        void write_bin_v(const vector<node_with_depth> v, const string fname, const string suffix){
            std::ofstream out{fname + "." + suffix, std::ios::binary};
            for(auto e : v)
                write(out, e);
        }

        /*
         count the number of parents from `v` to the lowest ancestor with a weiner link
         */
        size_type parent_depth(const cst_t& st, const node_type start_node, const char_type c){
            size_type d = 0;
            node_type u = start_node, v = st.double_rank_fail_wl(u, c);
            
            while((!st.is_root(u)) and st.is_root(v)){
                u = st.parent(u);
                v = st.double_rank_fail_wl(u, c);
                d += 1;
            }
            return d;
        }
        
        /*
         For every traversed internal node `alpha` compute the number of ancestor nodes `p` (`alpha`, `p`) l2
         For every traversed leaf `alpha` compute the number of ancestor nodes `p` (`alpha`, `p`) l1
         */
        size_type fill_vectors(const cst_t& st, vector<node_with_depth>& l1, vector<node_with_depth>& l2, const size_t sample_freq, size_type& max_depth1, size_type& max_depth2){
            auto start_time = timer::now();
            size_type nodes_visited = 0;
            
            max_depth1 = 0;
            max_depth2 = 0;

            node_type currnode = st.first_child(st.root()), nextnode = st.first_child(st.root());
            bool direction_down = true;
            size_type report_freq = (st.size() > 10 ? st.size() / 10 : st.size());

            do{
                if(direction_down){
                    if(!st.is_root(currnode)){
                        if(nodes_visited++ % report_freq == 0) // report progress
                            report_progress(start_time, currnode.i, st.size());
                        
                        if(static_cast<size_t>(sample_freq*static_cast<unsigned long>(std::rand())/(RAND_MAX+1UL)) == 0){ // sample
                            for(size_type i = 0; i < st.csa.sigma; i++){
                                char_type c = st.csa.comp2char[i];
                                bool is_leaf = st.is_leaf(currnode);
                                node_with_depth n = node_with_depth(currnode, parent_depth(st, currnode, c), c);
                                if(is_leaf){
                                    l2.push_back(n);
                                    max_depth2 = (max_depth2 > n.m_depth ? max_depth2 : n.m_depth);
                                } else {
                                    l1.push_back(n);
                                    max_depth1 = (max_depth1 > n.m_depth ? max_depth1 : n.m_depth);
                                }
                            }
                        }

                        // move on
                        nextnode = st.first_child(currnode);
                        if(st.is_root(nextnode))
                            direction_down = false;
                        else
                            currnode = nextnode;
                    } else {
                        direction_down = false;
                    }
                } else {
                    nextnode = st.sibling(currnode);
                    if(st.is_root(nextnode)) {
                        currnode = st.parent(currnode);
                    } else {
                        currnode = nextnode;
                        direction_down = true;
                    }
                }
            } while(!st.is_root(currnode));
            return std::chrono::duration_cast<std::chrono::milliseconds>(timer::now() - start_time).count();
        }

        size_type find_max_depth(vector<node_with_depth> v){
            size_type d = 0;
            for(auto n : v)
                d = (d > n.m_depth ? d : n.m_depth);
            return d;
        }

    public:
        vector<node_with_depth> m_vec1, m_vec2;
        size_type max_depth1, max_depth2;
        
        NwdList(const vector<node_with_depth>& v1, const vector<node_with_depth>& v2) : m_vec1{v1}, m_vec2{v2} {
            max_depth1 = find_max_depth(m_vec1);
            max_depth2 = find_max_depth(m_vec2);
        };

        NwdList(const cst_t& st, size_t sample_freq){
            vector<node_with_depth> l1, l2;
            size_type d1 = 0, d2 = 0;
            cerr << " * NwdList from tree of " << st.size() << " nodes, alphabet of size " << st.csa.sigma;
            cerr << ", sample 1 / " << sample_freq << " nodes";
            size_type duration = fill_vectors(st, l1, l2, sample_freq, d1, d2);
            cerr << " DONE (" << duration / 1000 << "seconds)" << endl;
            m_vec1 = l1;
            m_vec2 = l2;
            max_depth2 = d2;
            max_depth1 = d1;
        }
        
        string repr() const {
            return string("NwdLists(" +
                          std::to_string(m_vec1.size()) + ", " +
                          std::to_string(m_vec2.size()) + ")");
        }
        
        void shuffle_vectors(){
            std::random_shuffle(m_vec1.begin(), m_vec1.end());
            std::random_shuffle(m_vec2.begin(), m_vec2.end());
        }
        
        void write_bin(const string fname){
            cerr << "writing to " << fname;
            auto start_time = timer::now();
            write_bin_v(m_vec1, fname, "i.bin");
            write_bin_v(m_vec2, fname, "l.bin");
            size_type duration = std::chrono::duration_cast<std::chrono::milliseconds>(timer::now() - start_time).count();
            cerr << " DONE (" << duration / 1000 << " seconds)" << endl;
        }

        void check(const cst_t& st, const string& s, const size_type sample_freq){
            cerr << " ** checking coverage" << endl;
            size_type alp_size = st.csa.sigma;
            
            if(alp_size * s.size() != sample_freq * m_vec2.size()){
                if(sample_freq == 1)
                    throw string("ERROR: expecting " +
                                 to_string(alp_size * s.size()) + " leaves with, got " +
                                 to_string(m_vec2.size()));
                else
                    cerr << string("WARNING: expecting ") +
                            to_string(alp_size * s.size()) + string(" leaves with, got ") +
                            to_string(m_vec2.size()) << endl;
            }
            cerr << "OK" << endl;
        }

        bool operator==(const NwdList other){
            bool sizes = (m_vec1.size() == other.m_vec1.size() &&
                          m_vec2.size() == other.m_vec2.size());
            
            bool elements = true;
            for(int i=0; i<m_vec1.size(); i++)
                elements = (m_vec1[i] == other.m_vec1[i] ? elements : false);
            for(int i=0; i<m_vec2.size(); i++)
                elements = (m_vec2[i] == other.m_vec2[i] ? elements : false);
            
            return (sizes && elements);
        }

    };

    vector<node_with_depth> load_nwd_vector_bin(const string fname){
        std::ifstream in {fname, std::ios::binary};
        vector<node_with_depth> v;
        node_with_depth e;
        
        while(true){
            if(!read(in, e)) break;
            v.push_back(e);
        }
        in.close();
        cerr << " ** " << fname << " --> " << v.size() << " elements " << endl;
        return v;
    };
    
    NwdList load_nwd_list_bin(const string fname){
        cerr << " * loading from " << fname << endl;
        NwdList e = NwdList(load_nwd_vector_bin(fname + ".i.bin"),
                            load_nwd_vector_bin(fname + ".l.bin"));
        return e;
    }
    
}
#endif /* parent_depth_list_h */
