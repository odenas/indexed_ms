//
//  edge_list.hpp
//  fast_ms
//
//  Created by denas on 10/17/17.
//  Copyright © 2017 denas. All rights reserved.
//

#ifndef edge_list_h
#define edge_list_h


#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>

//#include <sdsl/vectors.hpp>

#include "stree_sct3.hpp"
#include "cst_iterator.hpp"
#include "maxrep_vector.hpp"

using namespace std;


namespace fdms {
    typedef StreeOhleb<>::node_type node_type;
    typedef StreeOhleb<>::size_type size_type;
    typedef uint8_t char_value;
    using timer = std::chrono::high_resolution_clock;
    
    /*
     basically a pair (node, char)
     */
    
    class edge_type{
    public:
        node_type m_node;
        char m_c;
        
        edge_type() {
            m_node = {0, 0, 0, 0, 0};
            m_c = 0;
        }
        
        edge_type(node_type v, char c) : m_node{v}, m_c{c}{};
        
        friend std::ostream& operator<<(std::ostream& out, const edge_type& edge){
            return out << std::setw(10) << edge.m_node.i << ' '
            << std::setw(10) << edge.m_node.j << ' '
            << std::setw(10) << edge.m_node.ipos << ' '
            << std::setw(10) << edge.m_node.cipos << ' '
            << std::setw(10) << edge.m_node.jp1pos << ' '
            << std::setw(10) << edge.m_c << '\n';
        }
        friend std::istream& operator>>(std::istream& in, edge_type& edge){
            return in >> edge.m_node.i >>
            edge.m_node.j >>
            edge.m_node.ipos >>
            edge.m_node.cipos >>
            edge.m_node.jp1pos >>
            edge.m_c;
        }
        
        friend std::ifstream& read(std::ifstream& in, edge_type& e){
            in.read(reinterpret_cast<char*>(&e.m_node.i), sizeof(e.m_node.i));
            in.read(reinterpret_cast<char*>(&e.m_node.j), sizeof(e.m_node.j));
            in.read(reinterpret_cast<char*>(&e.m_node.ipos), sizeof(e.m_node.ipos));
            in.read(reinterpret_cast<char*>(&e.m_node.cipos), sizeof(e.m_node.cipos));
            in.read(reinterpret_cast<char*>(&e.m_node.jp1pos), sizeof(e.m_node.jp1pos));
            in.read(reinterpret_cast<char*>(&e.m_c), sizeof(e.m_c));
            return in;
        };
        
        friend std::ofstream& write(std::ofstream& out, edge_type& e){
            out.write(reinterpret_cast<char*>(&e.m_node.i), sizeof(e.m_node.i));
            out.write(reinterpret_cast<char*>(&e.m_node.j), sizeof(e.m_node.j));
            out.write(reinterpret_cast<char*>(&e.m_node.ipos), sizeof(e.m_node.ipos));
            out.write(reinterpret_cast<char*>(&e.m_node.cipos), sizeof(e.m_node.cipos));
            out.write(reinterpret_cast<char*>(&e.m_node.jp1pos), sizeof(e.m_node.jp1pos));
            out.write(reinterpret_cast<char*>(&e.m_c), sizeof(e.m_c));
            return out;
        };
        
        bool operator==(const edge_type x){
            return (m_node.i == x.m_node.i &&
                    m_node.j == x.m_node.j &&
                    m_node.ipos == x.m_node.ipos &&
                    m_node.cipos == x.m_node.cipos &&
                    m_node.jp1pos == x.m_node.jp1pos &&
                    m_c == x.m_c);
        }
        
        string repr() const {
            return string("<(" +
                          std::to_string(m_node.i) + "," +
                          std::to_string(m_node.j) + "," +
                          std::to_string(m_node.ipos) + "," +
                          std::to_string(m_node.cipos) + "," +
                          std::to_string(m_node.jp1pos) + ") " +
                          m_c + ">");
        }
    };
    
    
    class EdgeList {
    private:
        void report_progress(timer::time_point start_time, size_type curr_idx, size_type total){
            timer::time_point stop_time = timer::now();
            size_type elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time).count() + 1;
            cerr << endl << "[" << elapsed / 1000 << " s] " << 100.0 * curr_idx / total << "% @ " << (1.0 * curr_idx / elapsed) << " KHz";
        }

        /*
         For every traversed node `alpha` try all Weiner links for all characters `c`.
         
         If at least 2 Weiner links are successful, add pairs (alpha,c) to the list of:
         - successful Weiner links for maximal nodes, where c is a char generating a successful Weiner link
         - unsuccessful Weiner links for maximal nodes, where c' is a character that generated a failed Weiner link.
         If just one character `c` gave a successful Weiner link, add the pair
         - (alpha,c) to the list of successful Weiner links for non maximal nodes
         - (alpha,c') to the list unsuccessful Weiner links for non maximal nodes, for all characters c' != c
         */
        size_type fill_vectors(const StreeOhleb<>& st, vector<edge_type>& l1, vector<edge_type>& l2, vector<edge_type>& l3, vector<edge_type>& l4,
                               const size_t sample_freq){
            
            auto start_time = timer::now();
            size_type nodes_visited = 0;

            node_type currnode = st.root(), nextnode = st.root();
            bool direction_down = true;
            
            do{
                if(direction_down){
                    if(!st.is_leaf(currnode)){
                        if(nodes_visited++ % (st.size() / 10) == 0) // report progress
                            report_progress(start_time, currnode.i, st.size());

                        if(static_cast<size_t>(sample_freq*static_cast<unsigned long>(std::rand())/(RAND_MAX+1UL)) == 0){ // sample
                            // process current node
                            if(Maxrep::rank_maximal_test(st, currnode)){ // maximal node
                                for(size_type i = 0; i < st.csa.sigma; i++){
                                    bool is_max = !st.is_root(st.double_rank_fail_wl(currnode, st.csa.comp2char[i]));
                                    (is_max ? l1 : l2).push_back(edge_type(currnode, st.csa.comp2char[i]));
                                }
                            } else {
                                char wl_sym = st.csa.char2comp[st.csa.bwt[currnode.j]];
                                for(size_type i = 0; i < st.csa.sigma; i++){
                                    (wl_sym == i ? l3 : l4).push_back(edge_type(currnode, st.csa.comp2char[i]));
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
        
        void write_bin_v(const vector<edge_type> v, const string fname, const string suffix){
            std::ofstream out{fname + "." + suffix, std::ios::binary};
            for(auto e : v)
                write(out, e);
        }
        
        void write_txt_v(const vector<edge_type> v, const string fname, const string suffix){
            std::ofstream out{fname + "." + suffix};
            for(int i=0; i<v.size(); i++)
                out << v[i] << endl;
        }
        
    public:
        vector<edge_type> m_vec1, m_vec2, m_vec3, m_vec4;

        EdgeList(const vector<edge_type> l1,   // haswl, maximal
                 const vector<edge_type> l2,   // nowl, maximal
                 const vector<edge_type> l3,   // haswl, nonmaximal
                 const vector<edge_type> l4   // nowl, nonmaximal
                 ) : m_vec1{l1}, m_vec2{l2}, m_vec3{l3}, m_vec4{l4} {};

        EdgeList(const StreeOhleb<>& st, size_t sample_freq){
            vector<edge_type> l1, l2, l3, l4;

            cerr << " * EdgeLists from tree of " << st.size() << " nodes, alphabet of size " << st.csa.sigma;
            cerr << ", sample 1 / " << sample_freq << " nodes";
            size_type duration = fill_vectors(st, l1, l2, l3, l4, sample_freq);
            cerr << " DONE (" << duration / 1000 << "seconds)" << endl;
            m_vec1 = l1;
            m_vec2 = l2;
            m_vec3 = l3;
            m_vec4 = l4;
            
        }
        
        vector<edge_type> vec(const bool is_maximal, const bool has_wl){
            if(is_maximal)
                return (has_wl ? m_vec1 : m_vec2);
            else
                return (has_wl ? m_vec3 : m_vec4);
        }
        
        string repr() const {
            return string("EdgeLists(" +
                          std::to_string(m_vec1.size()) + ", " +
                          std::to_string(m_vec1.size()) + ", " +
                          std::to_string(m_vec1.size()) + ", " +
                          std::to_string(m_vec1.size()) + ")");
        }
        
        void shuffle_vectors(){
            std::random_shuffle(m_vec1.begin(), m_vec1.end());
            std::random_shuffle(m_vec2.begin(), m_vec2.end());
            std::random_shuffle(m_vec3.begin(), m_vec3.end());
            std::random_shuffle(m_vec4.begin(), m_vec4.end());
        }
        
        void write_bin(const string fname){
            cerr << "writing to " << fname;
            auto start_time = timer::now();
            write_bin_v(m_vec1, fname, "11.bin");
            write_bin_v(m_vec2, fname, "01.bin");
            write_bin_v(m_vec3, fname, "10.bin");
            write_bin_v(m_vec4, fname, "00.bin");
            size_type duration = std::chrono::duration_cast<std::chrono::milliseconds>(timer::now() - start_time).count();
            cerr << " DONE (" << duration / 1000 << " seconds)" << endl;
        }
        
        void write_txt(const string fname){
            write_txt_v(m_vec1, fname, "11.txt");
            write_txt_v(m_vec2, fname, "01.txt");
            write_txt_v(m_vec3, fname, "10.txt");
            write_txt_v(m_vec4, fname, "00.txt");
        }
        
        bool operator==(const EdgeList other){
            bool sizes = (m_vec1.size() == other.m_vec1.size() &&
                          m_vec2.size() == other.m_vec2.size() &&
                          m_vec3.size() == other.m_vec3.size() &&
                          m_vec4.size() == other.m_vec4.size());
            bool elements = true;
            for(int i=0; i<m_vec1.size(); i++)
                elements = (m_vec1[i] == other.m_vec1[i] ? elements : false);
            for(int i=0; i<m_vec2.size(); i++)
                elements = (m_vec2[i] == other.m_vec2[i] ? elements : false);
            for(int i=0; i<m_vec3.size(); i++)
                elements = (m_vec3[i] == other.m_vec3[i] ? elements : false);
            for(int i=0; i<m_vec4.size(); i++)
                elements = (m_vec4[i] == other.m_vec4[i] ? elements : false);
            return (sizes && elements);
        }
        
        void check_with_maxrep(Maxrep& maxrep){
            cerr << " ** checking maximal nodes";
            for(auto v: m_vec1){
                if(!maxrep.is_maximal(v.m_node))
                    throw std::ios::failure(string("ERROR: on edge " + v.repr() + " " + maxrep.desc(v.m_node)));
            }
            cerr << " ... ";
            for(auto v: m_vec2){
                if(!maxrep.is_maximal(v.m_node))
                    throw std::ios::failure(string("ERROR: on edge " + v.repr() + " " + maxrep.desc(v.m_node)));
            }
            cerr << "OK" << endl;

            cerr << " ** checking non-maximal nodes";
            for(auto v: m_vec3){
                if(maxrep.is_maximal(v.m_node))
                    throw std::ios::failure(string("ERROR: on edge " + v.repr() + " " + maxrep.desc(v.m_node)));
            }
            cerr << " ... ";
            for(auto v: m_vec4){
                if(maxrep.is_maximal(v.m_node))
                    throw std::ios::failure(string("ERROR: on edge " + v.repr() + " " + maxrep.desc(v.m_node)));
            }
            cerr << "OK" << endl;
        }

    };
    
    vector<edge_type> load_edge_vector_bin(const string fname){
        std::ifstream in {fname, std::ios::binary};
        vector<edge_type> v;
        edge_type e;
        
        while(true){
            if(!read(in, e)) break;
            v.push_back(e);
        }
        in.close();
        cerr << " ** " << fname << " --> " << v.size() << " elements " << endl;
        return v;
    };
    
    EdgeList load_edge_list_bin(const string fname){
        cerr << " * loading from " << fname << endl;
        EdgeList e = EdgeList(load_edge_vector_bin(fname + ".11.bin"),
                              load_edge_vector_bin(fname + ".01.bin"),
                              load_edge_vector_bin(fname + ".10.bin"),
                              load_edge_vector_bin(fname + ".00.bin"));
        return e;
    }
    
    vector<edge_type> load_edge_vector_txt(const string fname){
        vector<edge_type> l1;
        edge_type v{};
        std::ifstream in {fname};
        if(!in){
            throw std::ios::failure{string("failed to open" + fname)};
        }
        while(true){
            in >> v;
            if(!in) break;
            l1.push_back(v);
        }
        in.close();
        return l1;
    }
    
    EdgeList load_edge_list_txt(const string fname){
        return EdgeList(load_edge_vector_txt(fname + ".11.txt"),
                        load_edge_vector_txt(fname + ".01.txt"),
                        load_edge_vector_txt(fname + ".10.txt"),
                        load_edge_vector_txt(fname + ".00.txt"));
    }
}



#endif /* edge_list_h */
