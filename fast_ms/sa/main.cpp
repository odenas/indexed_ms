//
//  main.cpp
//  sa
//
//  Created by denas on 12/12/16.
//  Copyright © 2016 denas. All rights reserved.
//

#include <iostream>
#include <string>

#include <sdsl/suffix_arrays.hpp>

using namespace sdsl;
using namespace std;



int main(int argc, const char * argv[]) {
    if (argc != 3){
        cout << "usage: " << argv[0] << " <input_file> <output_file>" << endl;
        exit(1);
    }

    csa_wt<wt_huff<>, 32, 32, sa_order_sa_sampling<>, isa_sampling<>, byte_alphabet> csa3;
    construct(csa3, argv[1], 1);
    store_to_file(csa3, argv[2]);
}
