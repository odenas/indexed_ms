//
//  main.cpp
//  vec_agg
//
//  Created by denas on 1/26/17.
//  Copyright © 2017 denas. All rights reserved.
//

/*
 given a list of files of dumped vectors produce a vector
 that is the element-wise max of the given vectors
 */


#include <iostream>
#include <sdsl/vectors.hpp>


int main(int argc, const char **argv) {
    sdsl::int_vector<32> MS;

    sdsl::load_from_file(MS, argv[1]);
    std::cerr << "processing " << argc - 1 << " files of length " << MS.size() << " ... " << std::endl;
    for(int i=2; i<argc; i++){
        sdsl::int_vector<32> curr_ms;
        sdsl::load_from_file(curr_ms, argv[i]);
        for(unsigned long long j=0; j<MS.size(); j++){
            if(MS[j] < curr_ms[j])
                MS[j] = curr_ms[j];
        }
    }

    std::cerr << "dumping ..." << std::endl;


    for(unsigned long long j=0; j<MS.size(); j++)
        std::cout << MS[j] << " ";
    std::cout << std::endl;
}
