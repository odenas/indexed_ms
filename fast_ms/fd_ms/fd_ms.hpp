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


extern "C" {
#include "dbwt.h"
}

using namespace std;
using namespace sdsl;

typedef unsigned long size_type;

class Bwt{
private:
    unsigned char *bwt;
    wt_huff<> wtree;

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

        string tmp_file = ram_file_name(util::to_string(util::pid()) + "_" + util::to_string(util::id()));
        store_to_file(*bwt, tmp_file);
        construct(wtree, tmp_file, 1);
        ram_fs::remove(tmp_file);
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

    void output_partial_vec(int_vector<> v, const int idx, const char name[], bool verbose){
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
        int_vector<> y(wtree.size() + 1);
        for(char s: alp){
            cout << s << ": " << endl;
            for(int i=0; i<wtree.size() + 1; i++)
                y[i] = wtree.rank(i, s);
            output_partial_vec(y, 0, "y", true);
        }

    }
};

class Interval{
private:
    Bwt *bwt;

public:
    size_type lb, ub;

    Interval(Bwt *bwt_, char c){
        bwt = bwt_;
        lb = bwt->C[bwt->char2int[c]];
        ub = bwt->C[bwt->char2int[c] + 1] - 1;
    }

    void set(size_type l, size_type u){
        lb = l;
        ub = u;
    }

    void bstep(char c){
        int cc = bwt->char2int[c];
        lb = bwt->C[cc] + bwt->rank(lb, c);
        ub = bwt->C[cc] + bwt->rank(ub + 1, c) - 1;
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

class Mstat{
private:
    bit_vector ms, runs;
    size_type ms_size;

    void output_partial_vec(const bit_vector v, size_type idx, const char name[], bool verbose){
        if (!verbose)
            return ;
        cout << name << ": ";
        for(size_type i = 0; i < (size_type)v.size(); i++)
            cout << v[i];
        cout << endl;
        for(size_type i = 0; i < (size_type)strlen(name) + 2; i++)
            cout << " ";
        for(size_type i = 0; i < (size_type)v.size(); i++)
            cout << (i == idx ? "*" : " ");
        cout << endl;
    }

    size_type find_k_prim(bit_vector runs, size_type k){
        size_t total_zeros = rank_support_v<0>(&runs)(runs.size());
        size_t zeros = rank_support_v<0>(&runs)(k + 1);

        if(total_zeros == zeros)
            return runs.size();

        bit_vector::select_0_type b_sel(&runs);
        return b_sel(zeros + 1);
    }

    bit_vector build_runs(string t, string s, Bwt bwt, const bool verbose){
        cst_sct3<> st_of_s;
        construct_im(st_of_s, s, 1);

        bit_vector runs(ms_size);
        size_type k = ms_size, c = t[k - 1];
        Interval I{&bwt, static_cast<char>(c)};

        cst_sct3<>::node_type v = st_of_s.child(st_of_s.root(), c); // stree node
        while(--k > 0){
            c = t[k-1];
            I.bstep(c);
            if(I.is_empty()){
                runs[k] = 0;
                // update I to the parent of the proper locus of w until we can extend by 'c'
                do{
                    v = st_of_s.parent(v);
                    I.set(v.i, v.j);
                    I.bstep(c);
                } while(I.is_empty());
                v = st_of_s.wl(v, c); // update v
            } else {
                v = st_of_s.wl(v, c); // update v
                runs[k] = 1;
            }
            output_partial_vec(runs, k, "runs", verbose);
        }
        return runs;
    }

    size_type get_ms(bit_vector ms_, size_type k){
        if(k == -1)
            return (size_type) 1;
        return bit_vector::select_1_type (&ms_)(k + 1) - (2 * k);
    }

    bit_vector build_ms(string t, string s_rev, Bwt bwt, const bool verbose){
        bit_vector ms(ms_size * 2);
        cst_sct3<> st_of_s;
        construct_im(st_of_s, s_rev, 1);
        size_type k = 0, h_star = k + 1, k_prim, ms_idx = 0;
        uint8_t c = t[k];
        Interval I{&bwt, static_cast<char>(c)};
        cst_sct3<>::node_type v = st_of_s.child(st_of_s.root(), c); // stree node

        while(k < ms_size){
            output_partial_vec(ms, ms_idx, "ms", verbose);
            for(; !I.is_empty() && h_star < ms_size; ){
                c = t[h_star];
                I.bstep(c);
                if(!I.is_empty()){
                    v = st_of_s.wl(v, c);
                    h_star ++;
                }
            }
            for(int i = 0; i < h_star - k - get_ms(ms, k - 1) + 1; i++)
                ms[ms_idx++] = 0;
            if(h_star - k - get_ms(ms, k - 1) + 1 > 0)
                ms[ms_idx++] = 1;

            if(h_star < ms_size){
                do {
                    v = st_of_s.parent(v);
                    I.set((int)v.i, (int)v.j);
                    I.bstep(t[h_star]);
                } while(I.is_empty());
                h_star = h_star + 1;
            }
            // k_prim: index of the first zero to the right of k in runs
            k_prim = find_k_prim(runs, k);

            for(size_type i = k + 1; i <= k_prim - 1; i++)
                ms[ms_idx++] = 1;

            // update v
            v = st_of_s.wl(v, c);
            k = k_prim;
        }
        output_partial_vec(ms, ms_idx, "ms", verbose);
        return ms;
    }

public:
    Mstat(string T, string S, Bwt& bwt_fw, string Srev, Bwt& bwt_rev, const bool verbose){
        //ms_size = std::strlen((char *)T);
        ms_size = T.size();
        runs = build_runs(T, S, bwt_fw, verbose);
        ms = build_ms(T, Srev, bwt_rev, verbose);
    }

    size_type operator[](size_type k){ return get_ms(ms, k); }

    void dump(){
        for (int i = 0; i < ms_size; i++)
            cout << (*this)[i];
        cout << endl;
    }
};



#endif /* fd_ms_h */
