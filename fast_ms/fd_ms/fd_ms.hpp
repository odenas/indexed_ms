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
#include <sdsl/suffix_trees.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "Bwt.hpp"
#include "stree.hpp"


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

}

#endif /* fd_ms_h */
