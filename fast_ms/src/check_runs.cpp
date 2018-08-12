/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   check_runs.cpp
 * Author: brt
 *
 * Created on May 15, 2018, 12:11 AM
 */

#include <cstdlib>

#include <sdsl/int_vector.hpp>

#include "fd_ms/opt_parser.hpp"
#include "fd_ms/input_spec.hpp"
#include "fd_ms/runs_ms.hpp"
#include "fd_ms/stree_sct3.hpp"
#include "fd_ms/query.hpp"


using namespace std;
using namespace fdms;

typedef StreeOhleb<>                cst_t;
typedef typename cst_t::node_type   node_type;
typedef typename cst_t::size_type   size_type;
typedef typename cst_t::char_type   char_type;
typedef sdsl::bit_vector            bitvec_t;
typedef MsVectors<cst_t, bitvec_t>  msvec_t;

cst_t st;


class InputFlags{
public:
    bool load_cst;

    InputFlags(){}

    InputFlags(const bool load_cst) : load_cst{load_cst} {}

    InputFlags(const InputFlags& i) : load_cst{i.load_cst} {}

    InputFlags(const OptParser args) : load_cst{args.getCmdOption("-load_cst") == "1"} {}
};

void load_or_build(const InputSpec& ispec, const bool load){
    string potential_stree_fname = ispec.fwd_cst_fname;
    if(load){
        std::cerr << " * loading the CST from " << potential_stree_fname << " ";
        sdsl::load_from_file(st, potential_stree_fname);
    } else {
        std::cerr << " * loadding  index string from " << ispec.s_fname << " " << endl;
        string s = ispec.load_s(false);
        std::cerr << " * building the CST of length " << s.size() << " ";
        sdsl::construct_im(st, s, 1);
    }
}

void comp(const InputSpec& tspec, const InputSpec& sspec, const InputFlags flags){
    size_type t_length = Query::query_length(tspec.s_fname);
    msvec_t ms_vec = msvec_t(t_length);

    load_or_build(sspec, flags.load_cst);
    ms_vec.fill_runs(tspec.s_fname, st, &cst_t::double_rank_fail_wl, &msvec_t::parent_sequence);
    sdsl::int_vector_buffer<1> runs(tspec.runs_fname, std::ios::in);
    for(int i=0; i<runs.size(); i++){
        cout << (runs[i] == ms_vec.runs[i] ? " " : "*") << i << "," << runs[i] << endl;
    }
}

int main(int argc, char** argv) {
    OptParser input(argc, argv);
    InputSpec s_spec, t_spec;
    InputFlags flags;

    if(argc == 1){
        const string base_dir = {"/home/brt/code/matching_statistics/indexed_ms/tests/datasets/testing/"};
        s_spec = InputSpec(base_dir + "mut_200s_64t_15.s");
        s_spec = InputSpec(base_dir + "mut_200s_64t_15.t");
        flags = InputFlags(false);
    } else {
        s_spec = InputSpec(input.getCmdOption("-s_path"));
        t_spec = InputSpec(input.getCmdOption("-t_path"));
        flags = InputFlags(input);
    }
    comp(t_spec, s_spec, flags);
    return 0;
}

