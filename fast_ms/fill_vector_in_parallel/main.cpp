//
//  main.cpp
//  fill_vector_in_parallel
//
//  Created by denas on 3/13/17.
//  Copyright © 2017 denas. All rights reserved.
//

#include <iostream>
#include <string>
#include <future>
#include <thread>

#include <sdsl/suffix_trees.hpp>

#include "basic.hpp"
#include "fd_ms.hpp"
#include "stree_sada.hpp"
#include "stree_sct3.hpp"

using namespace fdms;

sdsl::int_vector<8> runs(1);
size_type thread_id = 0;

int fill_runs(size_type i, size_type j){
    cout << "[" << i << ", " << j << ")" << endl;
    for(; i<j; i++)
        runs[i] = i;
    return 0;
}

void comp(InputSpec& T, InputSpec& S_fwd, const string& out_path, InputFlags& flags){
    const size_type nr_threads = flags.ms_progress;
    string t = T.load_s();
    runs.resize(t.size());
    if(runs.size() % nr_threads != 0)
        exit(1);

    for(size_type i=0; i<runs.size(); i++)
        runs[i] = 0;
    size_type chunk_size = runs.size() / nr_threads;

    std::vector<std::future<int>> results(nr_threads);

    for(size_type i=0; i<nr_threads; i++)
        results[i] = std::async(std::launch::async, fill_runs, i * chunk_size, (i + 1) * chunk_size);

    for(size_type i=0; i<nr_threads; i++)
        results[i].get();


    // print
    for(size_type i=0; i<runs.size(); i++)
        cout << "runs[" << i << "] = " << (int) runs[i] << endl;
}

int main(int argc, char **argv) {
    InputParser input(argc, argv);
    if(argc == 1){
        const string base_dir = {"/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/test_input_data/"};
        InputFlags flags(false, // lazy_wl
                         false, // sada cst
                         false, // space
                         true, // time
                         false, // ans
                         false, // verbose
                         0,     // nr. progress messages for runs construction
                         10,     // nr. progress messages for ms construction
                         false  // load CST
                         );
        InputSpec tspec(base_dir + "abcde200_16384t.txt");
        InputSpec sfwd_spec(base_dir + "abcde200_16384s.txt");
        const string out_path = "";
        comp(tspec, sfwd_spec, out_path, flags);
    } else {
        InputFlags flags(input);
        InputSpec tspec(input.getCmdOption("-t_path"));
        InputSpec sfwd_spec(input.getCmdOption("-s_path"));
        const string out_path = input.getCmdOption("-out_path");
        comp(tspec, sfwd_spec, out_path, flags);
    }
    return 0;
}
