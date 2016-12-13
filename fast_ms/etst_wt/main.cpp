
#include <iostream>
#include <string>
//#include <sdsl/suffix_trees.hpp>
//#include <sdsl/construct_sa.hpp>

#include "fd_ms.hpp"

using namespace std;
using namespace fdms;

//typedef csa_wt<> csa_t;

int main(int argc, char* argv[])
{
    string fname = "/Users/denas/Desktop/FabioImplementation/software/indexed_ms/fast_ms/etst_wt/main.cpp";

    //csa_t csa;
    //construct(csa, fname, 1);

    StreeOhleb<> cst;
    construct(cst, fname, 1);

    cout << cst.size() << endl;
    
}
