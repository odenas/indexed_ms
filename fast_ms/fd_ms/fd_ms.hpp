//
//  fd_ms.hpp
//  fast_ms
//
//  Created by denas on 11/2/16.
//  Copyright © 2016 denas. All rights reserved.
//

#ifndef fd_ms_h
#define fd_ms_h

#include <iostream>
#include <string>
#include <vector>
#include <sdsl/suffix_trees.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "Bwt.hpp"
#include "stree_sada.hpp"
#include "stree_sct3.hpp"


using namespace std;
using namespace sdsl;

typedef unsigned long size_type;


namespace fdms{


class InputSpec{
private:
    bit_vector parse_bitstr(string& s){
        bit_vector b(s.size());

        for(size_type i = 0; i < s.size(); i++)
            b[i] = ((unsigned char)s[i] - 48);
        return b;
    }
public:
    string s_fname, sbp_fname;

    InputSpec(string s_fn, string sbp_fn) : s_fname(s_fn), sbp_fname(sbp_fn){}

    string load_s(){
        string s;
        std::ifstream s_file {s_fname};
        while(s_file >> s)
            ;
        return s;
    }

    bit_vector load_bps(){
        string sbp;
        std::ifstream sbp_file {sbp_fname};
        sbp_file >> sbp;
        return parse_bitstr(sbp);
    }
};


class Interval{
private:
    Bwt& bwt;

public:
    size_type lb, ub;

    Interval(Bwt& bwt_, char c) : bwt{bwt_} {
        lb = bwt.C[bwt.char2int[c]];
        ub = bwt.C[bwt.char2int[c] + 1] - 1;
    }

    void set(size_type l, size_type u){
        lb = l;
        ub = u;
    }

    void bstep(char c){
        int cc = bwt.char2int[c];
        lb = bwt.C[cc] + bwt.rank(lb, c);
        ub = bwt.C[cc] + bwt.rank(ub + 1, c) - 1;
    }

    bool is_empty(){
        return lb > ub;
    }

    void dump(){
        std::cout << "{";
        if(is_empty())
            std::cout << "-, -";
        else
            std::cout << lb << ", " << ub;
        std::cout << "}" << endl;
    }
};

// copied from http://stackoverflow.com/questions/865668/how-to-parse-command-line-arguments-in-c
class InputParser{
public:
    std::string empty = "";
    InputParser (int &argc, char **argv){
        for (int i=1; i < argc; ++i)
            this->tokens.push_back(std::string(argv[i]));
    }
    /// @author iain
    const std::string& getCmdOption(const std::string &option) const{
        std::vector<std::string>::const_iterator itr;
        itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
        if (itr != this->tokens.end() && ++itr != this->tokens.end()){
            return *itr;
        }
        return empty;
    }
    /// @author iain
    bool cmdOptionExists(const std::string &option) const{
        return std::find(this->tokens.begin(), this->tokens.end(), option)
        != this->tokens.end();
    }
private:
    std::vector <std::string> tokens;
};

class InputFlags{
public:
    bool space_usage, mach_space_usage;
    bool time_usage;
    bool answer;
    bool verbose;

    InputFlags(bool space, bool mach_space, bool time_, bool ans, bool v) :
        space_usage {space},
        mach_space_usage {mach_space},
        time_usage {time_},
        answer {ans},
        verbose{v}
    {}

    InputFlags (InputParser input) :
        space_usage {input.getCmdOption("-s") == "1"},      // space usage
        mach_space_usage {input.getCmdOption("-S") == "1"}, // resident/vm memory usage
        time_usage {input.getCmdOption("-t") == "1"},       // time usage
        answer {input.getCmdOption("-a") == "1"},           // answer
        verbose{input.getCmdOption("-v") == "1"}            // verbose
    {}
        
};

}

#endif /* fd_ms_h */
