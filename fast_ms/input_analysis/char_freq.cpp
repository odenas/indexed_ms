#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <sdsl/vectors.hpp>
#include "fd_ms/opt_parser.hpp"
#include "fd_ms/query.hpp"

using namespace std;
typedef unsigned long size_type;
typedef map<char, size_type> counter_t;

int main(int argc, char **argv) {
    fdms::OptParser input(argc, argv);
    string path(input.getCmdOption("-path"));
    counter_t freq{};

    fdms::Query_fwd a{path, 1024*1024};

    cout << path << ": " << a.size() << endl;
    for(size_type j=0; j < a.size(); j++)
        freq[a[j]] += 1;

    cout << "**" << endl;

    for (auto item : freq)
        cout << item.first << ": " << item.second << endl;

    return 0;
}

