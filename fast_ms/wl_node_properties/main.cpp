//
//  main.cpp
//  wl_node_properties
//
//  Created by denas on 10/15/17.
//  Copyright © 2017 denas. All rights reserved.
//

#include <iostream>
#include "sdsl/int_vector.hpp"
#include "utils.hpp"
#include "cmd_utils.hpp"
#include "cst_iterator.hpp"

using namespace std;
using namespace fdms;

typedef StreeOhleb<>::node_type node_type;
typedef uint8_t char_value;


class edge_type{
private:
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
};

class EdgeList {
private:
    /*
     For every traversed node `alpha` try all Weiner links for all characters `c`.
     
     If at least 2 Weiner links are successful, add pairs (alpha,c) to the list of:
     - successful Weiner links for maximal nodes, where c is a char generating a successful Weiner link
     - unsuccessful Weiner links for maximal nodes, where c' is a character that generated a failed Weiner link.
     If just one character `c` gave a successful Weiner link, add the pair
     - (alpha,c) to the list of successful Weiner links for non maximal nodes
     - (alpha,c') to the list unsuccessful Weiner links for non maximal nodes, for all characters c' != c
     */
    void fill_vectors(const StreeOhleb<>& st, vector<edge_type>& l1, vector<edge_type>& l2, vector<edge_type>& l3, vector<edge_type>& l4){
        size_type nodes_visited = 0;
        sdsl::bit_vector wl_presence(st.csa.sigma);
        for(cst_dfslr_iterator<StreeOhleb<>> pos(&st); !pos.end(); ++pos, nodes_visited++){
            node_type currnode = *pos;
            
            for(size_type i = 0; i < st.csa.sigma; i++)
                wl_presence[i] = st.has_wl(currnode, st.csa.comp2char[i]);
            size_type wl_count = sdsl::util::cnt_one_bits(wl_presence);
            
            if(wl_count >= 2){ // maximal node
                for(char_value i = 0; i < st.csa.sigma; i++){
                    (wl_presence[i] ? l1 : l2).push_back(edge_type(currnode, st.csa.comp2char[i]));
                }
            } else if (wl_count == 1){  // non-maximal node
                for(size_type i = 0; i < st.csa.sigma; i++){
                    (wl_presence[i] ? l3 : l4).push_back(edge_type(currnode, st.csa.comp2char[i]));
                }
            } else
                throw std::ios::failure(string("grrr. node with no WL"));
        }
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
             const vector<edge_type> l4) : // nowl, nonmaximal
    m_vec1{l1}, m_vec2{l2}, m_vec3{l3}, m_vec4{l4}{};
    
    EdgeList(const StreeOhleb<>& st){
        vector<edge_type> l1, l2, l3, l4;
        fill_vectors(st, l1, l2, l3, l4);
        m_vec1 = l1;
        m_vec2 = l2;
        m_vec3 = l3;
        m_vec4 = l4;
    }
    
    void shuffle_vectors(){
        std::random_shuffle(m_vec1.begin(), m_vec1.end());
        std::random_shuffle(m_vec2.begin(), m_vec2.end());
        std::random_shuffle(m_vec3.begin(), m_vec3.end());
        std::random_shuffle(m_vec4.begin(), m_vec4.end());
    }

    void write_bin(const string fname){
        write_bin_v(m_vec1, fname, "11.bin");
        write_bin_v(m_vec2, fname, "01.bin");
        write_bin_v(m_vec3, fname, "10.bin");
        write_bin_v(m_vec4, fname, "00.bin");
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
    return v;
};

EdgeList load_edge_list_bin(const string fname){
    return EdgeList(load_edge_vector_bin(fname + ".11.bin"),
                    load_edge_vector_bin(fname + ".01.bin"),
                    load_edge_vector_bin(fname + ".10.bin"),
                    load_edge_vector_bin(fname + ".00.bin"));
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

void comp(InputSpec& S_fwd, const string& out_path, InputFlags& flags){
    string s = S_fwd.load_s();
    /* build stree */
    StreeOhleb<> st;
    size_type load_cst_time = load_st(st, s, S_fwd.fwd_cst_fname, flags.load_stree);
    cerr << "DONE (" << load_cst_time / 1000 << "seconds)" << endl;

    EdgeList e{st};
    e.write_txt(out_path);
    EdgeList f {load_edge_list_txt(out_path)};
    cout << (e == f) << endl;

    e.write_bin(out_path);
    EdgeList g {load_edge_list_bin(out_path)};
    cout << (e == g) << endl;

}

int main(int argc, char **argv){
    OptParser input(argc, argv);
    if(argc == 1){
        const string base_dir = {"/Users/denas/projects/matching_statistics/indexed_ms/tests/datasets/testing/"};
        InputFlags flags(false, // lazy_wl
                         false,  // rank-and-fail
                         false,  // use maxrep
                         true,  // lca_parents
                         false, // space
                         false, // time
                         true,  // ans
                         false, // verbose
                         10,    // nr. progress messages for runs construction
                         10,    // nr. progress messages for ms construction
                         false, // load CST
                         false, // load MAXREP
                         1      // nthreads
                         );
        InputSpec sfwd_spec(base_dir + "rep_100000s_1000t_ab_1_5.s");
        const string out_path = "ciao";
        comp(sfwd_spec, base_dir + out_path, flags);
    } else {
        InputFlags flags(input);
        InputSpec sfwd_spec(input.getCmdOption("-s_path"));
        const string out_path = input.getCmdOption("-out_path");
        comp(sfwd_spec, out_path, flags);
    }
    
    return 0;
}
