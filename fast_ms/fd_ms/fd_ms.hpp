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
#include "fd_ms.hpp"


using namespace std;
using namespace sdsl;


class Bwt{
private:
    std::string *bwt;
    wt_huff<> wtree;

    /**
     * parse a string of the type: 4-12-31 into an array [4, 12, 31]
     */
    int *parse_C(std::string Cstr){
        int *C = (int *)malloc(sizeof(int) * (Cstr.size() + 1) / 2 );
        int i = 0, j = 0, k = 0;

        while((j = (int)Cstr.find('-', i)) > 0){
            C[k++] = atoi(Cstr.substr(i, j - i).c_str());
            i = ++j;
        }
        C[k] = atoi(Cstr.substr(i, Cstr.size() - i).c_str());
        return C;
    }

    int *parse_alphabet(const std::string A){
        int *char2int = (int *)malloc(sizeof(int) * 128);
        for(int i=0; i<A.size(); i++)
            char2int[A[i]] = i;
        return char2int;
    }

public:
    const int *C;
    const int *char2int;

    Bwt(std::string *bwt, const int *C, const int *c2i){
        this->bwt = bwt;
        this->C = C;
        char2int = c2i;
    }
    Bwt(std::string *bwt, const string *Cstr, const string *A){
        this->bwt = bwt;
        C = parse_C(*Cstr);
        char2int = parse_alphabet(*A);

        string tmp_file = ram_file_name(util::to_string(util::pid()) + "_" + util::to_string(util::id()));
        //cout << tmp_file << endl;
        store_to_file(*bwt, tmp_file);
        construct(wtree, tmp_file, 1);
        ram_fs::remove(tmp_file);
    }

    int rank(int i, char c){
        return (int) wtree.rank(i, c);
    }

    int rrank(int i, char c){
        int cnt = 0;
        for(char cc: bwt->substr(0, i)){
            if(cc == c)
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
    int lb, ub;

    Interval(Bwt *bwt_, char c){
        bwt = bwt_;
        lb = bwt->C[bwt->char2int[c]];
        ub = bwt->C[bwt->char2int[c] + 1] - 1;
    }

    void set(int l, int u){
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
    uint ms_size;

    void output_partial_vec(const bit_vector v, const int idx, const char name[], bool verbose){
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

    int find_k_prim(bit_vector runs, int k){
        size_t total_zeros = rank_support_v<0>(&runs)(runs.size());
        size_t zeros = rank_support_v<0>(&runs)(k+1);
        if(total_zeros == zeros)
            return (int) (runs.size());
        bit_vector::select_0_type b_sel(&runs);
        return (int) b_sel(zeros + 1);
    }

    bit_vector build_runs(string t, string s, Bwt bwt, const bool verbose){
        cst_sct3<> st_of_s;
        construct_im(st_of_s, s, 1);

        bit_vector runs(t.size());
        uint8_t k = t.size(), c = t[k - 1];
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
                    I.set((int) v.i, (int) v.j);
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

    uint get_ms(bit_vector ms_, uint k){
        if(k == -1)
            return (uint) 1;
        return ((uint) bit_vector::select_1_type (&ms_)(k + 1)) - (2 * k);
    }

    bit_vector build_ms(string t, string s_rev, Bwt bwt, const bool verbose){
        bit_vector ms(t.size() * 2);
        cst_sct3<> st_of_s;
        construct_im(st_of_s, s_rev, 1);
        unsigned int k = 0, h_star = k + 1, k_prim, ms_idx = 0;
        uint8_t c = t[k];
        Interval I{&bwt, static_cast<char>(c)};
        cst_sct3<>::node_type v = st_of_s.child(st_of_s.root(), c); // stree node

        while(k < t.size()){
            output_partial_vec(ms, ms_idx, "ms", verbose);
            for(; !I.is_empty() && h_star < t.size(); ){
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

            if(h_star < t.size()){
                do {
                    v = st_of_s.parent(v);
                    I.set((int)v.i, (int)v.j);
                    I.bstep(t[h_star]);
                } while(I.is_empty());
                h_star = h_star + 1;
            }
            // k_prim: index of the first zero to the right of k in runs
            k_prim = find_k_prim(runs, k);

            for(int i=k+1; i<=k_prim - 1; i++)
                ms[ms_idx++] = 1;

            // update v
            v = st_of_s.wl(v, c);
            k = k_prim;
        }
        output_partial_vec(ms, ms_idx, "ms", verbose);
        return ms;
    }

public:
    Mstat(string T, string S, Bwt bwt_fw, string Srev, Bwt bwt_rev, const bool verbose){
        runs = build_runs(T, S, bwt_fw, verbose);
        ms = build_ms(T, Srev, bwt_rev, verbose);
        ms_size = (uint)T.length();
    }

    uint operator[](int k){ return get_ms(ms, k); }

    void dump(){
        for (int i = 0; i < ms_size; i++)
            cout << (*this)[i];
        cout << endl;
    }
};



#endif /* fd_ms_h */
