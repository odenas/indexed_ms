//
//  main.cpp
//  test_select_support
//
//  Created by denas on 1/24/17.
//  Copyright © 2017 denas. All rights reserved.
//


/*
 building a select_support datastructure is costly. test this on given inputs
*/

#include <iostream>
#include <sdsl/suffix_trees.hpp>

#include "fd_ms.hpp"


using namespace std;
using namespace fdms;

using timer = std::chrono::high_resolution_clock;


monitor::size_dict do_test(const string& t){
    monitor::size_dict time_usage;
    bvector ms(t.size());
    unsigned long size_in_bytes_ms_select1 = 0;
    unsigned long k = 0;

    auto runs_start = timer::now();
    while(++k < t.size()){
        sdsl::select_support_mcl<1,1> ms_select1(&ms);
        size_in_bytes_ms_select1 = (size_in_bytes_ms_select1 < sdsl::size_in_bytes(ms_select1) ?
                                    sdsl::size_in_bytes(ms_select1) : size_in_bytes_ms_select1);
        switch(t[k]){
            case 'A': ms[k] = 0;
                break;
            case 'C': ms[k] = 1;
                break;
            case 'G': ms[k] = 0;
                break;
            case 'T': ms[k] = 1;
                break;
            default: ms[k] = 0;
        }
    }
    auto runs_stop = timer::now();
    cout << std::chrono::duration_cast<std::chrono::milliseconds>(runs_stop - runs_start).count() << endl;

    return time_usage;
}

int main(int argc, char **argv) {
    InputParser input(argc, argv);

    if(argc == 1){
        const string base_dir = {"/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/genome_tests/genome_data/"};
        InputFlags flags(false, // lazy_wl
                         false, // sada cst
                         false, // space
                         false, // time
                         true,  // ans
                         true,  // verbose
                         10,    //nr. progress messages for runs construction
                         5000,  //nr. progress messages for ms construction
                         true   // load CST
                         );
        InputSpec tspec(base_dir + "../Homo_sapiens.GRCh38.dna.chromosome.22.juststring");
        //InputSpec sfwd_spec(base_dir + "Mus_musculus.GRCm38.dna.chromosome.MT.juststring");
    } else {
        InputFlags flags(input);
        InputSpec tspec(input.getCmdOption("-t_path"));
        //InputSpec sfwd_spec(input.getCmdOption("-s_path"));
        do_test(tspec.load_s());
    }
    return 0;
}
