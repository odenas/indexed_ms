#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <future>
#include <thread>

#include "fd_ms/opt_parser.hpp"
#include "fd_ms/input_spec.hpp"
#include "fd_ms/help.hpp"

using namespace std;
using namespace fdms;

string load_s(const string s_fname) {
    string s;
    std::ifstream s_file{s_fname};
    while (s_file >> s)
        ;
    return s;
}

/**
 * longest prefix of t that occurs in s
 */
int comp(const string& s, const string& t){
    for(int l = t.size(); l > 0; l--){
        if(s.find(t.substr(0, l)) != string::npos)
            return l;
    }
    return 0;
}

int main(int argc, char **argv) {
    OptParser input(argc, argv);
    string s = load_s(input.getCmdOption("-s_path"));
    string t = load_s(input.getCmdOption("-t_path"));

    for(int i = 0; i < t.size(); i++){
        cout << comp(s, t.substr(i, t.size() - i));
        if(i < t.size() - 1)
            cout << " ";
    }
    cout << endl;
    return 0;
}
