//
//  Bwt.hpp
//  fast_ms
//
//  Created by denas on 11/15/16.
//  Copyright © 2016 denas. All rights reserved.
//

#ifndef Bwt_h
#define Bwt_h

#include <string>
#include <iostream>
#include <sdsl/int_vector.hpp>
#include <sdsl/wavelet_trees.hpp>
extern "C" {
#include "dbwt.h"
}

using namespace std;

namespace fdms{
class Bwt{
private:
    unsigned char *bwt;
    sdsl::wt_huff<> wtree;
    typedef unsigned long size_type;

    void parse_alphabet(const unsigned char *S, size_type s_len){
        // char2int[char] = 1 iff char occurs in S
        for(int i = 0; i < s_len; i++){
            unsigned char c = S[i];
            if(char2int[c] == 0){
                char2int[c] = 1;
                sigma++;
            }
        }
        assert (char2int['#'] == 0);
        char2int['#'] = 1;

        // char2int[s] = rank of s
        for(int i = 0, j = 0; i < 128; i++)
            if(char2int[i] > 0)
                char2int[i] = j++;
    }

    void computeC(const unsigned char *S, size_type s_len){
        size_type cnt[sigma];
        for(int i = 0; i < sigma; i++)
            cnt[i] = 0;

        for(int i = 0; i < s_len; i++)
            cnt[char2int[S[i]]] += 1;
        cnt[char2int['#']] = 1;

        for(int i = 1; i <= sigma; i++)
            C[i] = C[i - 1] + cnt[i - 1];
        assert (C[sigma] == s_len + 1);
    }

public:
    size_type C[128], bwt_len;
    uint8_t char2int[128];
    uint8_t sigma = 1; // 0 reserved for '#'

    //Bwt(){}

    Bwt(const unsigned char *S){
        for(int i=0; i<128; i++)
            C[i] = char2int[i] = 0u;
        unsigned int last = 0;
        size_type s_len = std::strlen((const char *)S);
        bwt_len = s_len + 1;

        parse_alphabet(S, s_len);
        computeC(S, s_len);

        // bwt
        bwt = dbwt_bwt((unsigned char *)S, (long)s_len, &last, 0u);
        bwt[last] = '#';

        // construct the wavelet tree on the bwt char sequence
        string tmp_file = sdsl::ram_file_name(sdsl::util::to_string(sdsl::util::pid()) + "_" + sdsl::util::to_string(sdsl::util::id()));
        sdsl::store_to_file(*bwt, tmp_file);
        construct(wtree, tmp_file, 1);
        sdsl::ram_fs::remove(tmp_file);
    }

    size_type rrank(size_type i, char c){
        return wtree.rank(i, c);
    }

    size_type rank(size_type i, char c){
        size_type cnt = 0;
        for(int j=0; j<i; j++){
            if(bwt[j] == c)
                cnt++;
        }
        return cnt;
    }

    void output_partial_vec(sdsl::int_vector<> v, const int idx, const char name[], bool verbose){
        if (!verbose)
            return ;
        cout << name << ": ";
        for(int i=0; i<v.size(); i++)
            cout << v[i];
        cout << endl;
        for(int i=0; i<strlen(name) + 2; i++)
            cout << " ";
        for(int i=0; i<v.size(); i++)
            cout << (i == idx ? "*" : " ");
        cout << endl;
    }
        
    void show_bwt(const string alp){
        for(int i=0; i<wtree.size(); i++)
            cout << (char)wtree[i] << " ";
        cout << endl;
        sdsl::int_vector<> y(wtree.size() + 1);
        for(char s: alp){
            cout << s << ": " << endl;
            for(int i=0; i<wtree.size() + 1; i++)
                y[i] = wtree.rank(i, s);
            output_partial_vec(y, 0, "y", true);
        }
    }
};

}
#endif /* Bwt_h */
