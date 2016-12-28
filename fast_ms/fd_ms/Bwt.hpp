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

#include "basic.hpp"

using namespace std;

namespace fdms{
    class Bwt{
    private:
        sdsl::wt_huff<> wtree;
        unsigned char *bwt;

        void parse_alphabet(string& S, size_type s_len){
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
            alphabet = new char[sigma];

            // char2int[s] = rank of s
            for(int i = 0, j = 0; i < 128; i++)
                if(char2int[i] > 0){
                    alphabet[j] = (char)i;
                    char2int[i] = j++;
                }
        }

        void computeC(string& S, size_type s_len){
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
        size_type size_in_bytes__wtree, size_in_bytes__bwt, size_in_bytes__C, size_in_bytes__char2int, size_in_bytes__Sigma;
        uint8_t char2int[128];
        uint8_t sigma = 1; // 0 reserved for '#'
        char *alphabet;

        Bwt(string& S){
            for(int i=0; i<128; i++)
                C[i] = char2int[i] = 0u;
            unsigned int last = 0;
            size_type s_len = S.size();
            bwt_len = s_len + 1;

            parse_alphabet(S, s_len);
            computeC(S, s_len);

            // bwt
            bwt = dbwt_bwt((unsigned char *)S.c_str(), (long)s_len, &last, 0u);
            bwt[last] = '#';
            string bbwt((char *) bwt);

            // construct the wavelet tree on the bwt char sequence
            string tmp_file = sdsl::ram_file_name(sdsl::util::to_string(sdsl::util::pid()) + "_" + sdsl::util::to_string(sdsl::util::id()));
            sdsl::store_to_file(bbwt, tmp_file);
            construct(wtree, tmp_file, 1);
            sdsl::ram_fs::remove(tmp_file);

            size_in_bytes__wtree = sdsl::size_in_bytes(wtree);
            size_in_bytes__bwt = bwt_len;
            size_in_bytes__C = 128;
            size_in_bytes__char2int = 128;
            size_in_bytes__Sigma = sigma;
        }

        sdsl::wt_huff<>& get_wtree(){ return wtree; }

        std::pair<size_type, size_type> double_rank(size_type i, size_type j, char c){
            return wtree.double_rank(i, j, c);
        }

        /*
         * The number of occurrences of symbol c in the prefix [0..i-1]
         * for 0 <= i <= size()
         */
        size_type rank(size_type i, char c){
            return wtree.rank(i, c);
        }

        size_type rrank(size_type i, char c){
            size_type cnt = 0;
            for(int j = 0; j < i; j++){
                if(wtree[j] == c)
                    cnt++;
            }
            return cnt;
        }

        /*
         * largest j: c occurs i times in bwt[0:j-1] -- not so sure about this.
         *
         * (from sdsl) The index i in [0..size()-1] of the j-th occurrence of symbol c
         * 0 < j <= occurrence(c) in wtree
         */
        size_type select(size_type j, char c){
            return wtree.select(j, c);
        }

        size_type sselect(size_type i, char c){
            size_type cnt = 0;
            for(int j = 0; j < bwt_len; j++){
                if(wtree[j] == c)
                    cnt++;
                if(cnt > i)
                    return cnt;
            }
            return wtree.size();
        }

        size_type lf(size_type i){
            char c = wtree[i];
            return C[char2int[c]] + rank(i, c);
        }

        size_type lf_rev(size_type i){
            int j = 0;
            do{
                if(C[j] <= i && i < C[j + 1])
                    break;
                j++;
            } while(j < sigma);
            return wtree.select(i - C[j] + 1, alphabet[j]);
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

        char operator[](size_type i) { return wtree[i]; }
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

    void test(){
        typedef unsigned long size_type;

        string Sfwd {"aabbaba"};
        Bwt Bwtfwd(Sfwd);

        for(size_type i=0; i < Bwtfwd.bwt_len; i++)
            cout << Bwtfwd[i] << endl;

        cout << "rank, 'a': ";
        for(size_type i=0; i <= Bwtfwd.bwt_len; i++)
            cout << (int) Bwtfwd.rank(i, 'a');
        cout << endl << "select, 'a': ";
        for(size_type i=1; i <= Bwtfwd.rank(Bwtfwd.bwt_len, 'a'); i++)
            cout << (int) Bwtfwd.select(i, 'a');
        cout << endl;

        cout << "rank, 'b': ";
        for(size_type i=0; i <= Bwtfwd.bwt_len; i++)
            cout << (int) Bwtfwd.rank(i, 'b');
        cout << endl << "select, 'b': ";
        for(size_type i=1; i <= Bwtfwd.rank(Bwtfwd.bwt_len, 'b'); i++)
            cout << (int) Bwtfwd.select(i, 'b');
        cout << endl;


        for(size_type i=0; i<Bwtfwd.bwt_len; i++)
            cout << "LF(" << i << ") = " << Bwtfwd.lf(i) << endl;

        for(size_type i=0; i<Bwtfwd.bwt_len; i++)
            cout << "LF^-1(" << i << ") = " << Bwtfwd.lf_rev(i) << endl;
        
        for(size_type i=0; i<Bwtfwd.bwt_len; i++)
            cout << "LF(LF^-1(" << i << ")) = " << Bwtfwd.lf(Bwtfwd.lf_rev(i)) << endl;
    }

}
#endif /* Bwt_h */
