/*
 * dump a CST to a file. to be used for later use
 */

#include <iostream>
#include <fstream>
#include <vector>

#include "fd_ms/input_spec.hpp"
#include "fd_ms/opt_parser.hpp"
#include "fd_ms/stree_sct3.hpp"


using namespace std;
using namespace fdms;
using timer = std::chrono::high_resolution_clock;

typedef StreeOhleb<>                cst_t;
typedef typename cst_t::size_type size_type;
typedef typename cst_t::char_type char_type;
typedef typename cst_t::node_type node_type;

void dfs(const cst_t& st, const node_type& v){
    cout << "(";
    node_type u = st.first_child(v);
    while(u != st.root()) {
        dfs(st, u);
        u = st.sibling(u);
    }
    cout << ")";
}

void comp(const InputSpec& s_spec){
    StreeOhleb<> st;
    StreeOhleb<>::size_type load_cst_time = cst_t::load_or_build(st, s_spec, false, false);
    cerr << "DONE (" << load_cst_time / 1000 << "seconds)" << endl;

    dfs(st, st.root());
    cout << endl;
}


int main(int argc, char **argv){
    OptParser input(argc, argv);
    InputSpec ispec;

    if(argc == 1){
        const string base_dir = {"/home/brt/code/matching_statistics/indexed_ms/fast_ms/tree_pruning/"};
        ispec = InputSpec(base_dir + "input.s", "");
    } else {
        ispec = InputSpec(input.getCmdOption("-s_path"), "");
    }
    comp(ispec);
    return 0;
}

