/*
 * fabio_djamal_ms.cpp
 *
 *  Created on: Oct 13, 2016
 *      Author: denas
 */

#include <iostream>
#include <string>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_trees.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/bp_support.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/bit_vectors.hpp>

using namespace std;
using namespace sdsl;

class Bwt{
private:
    string bwt;

    /**
     * parse a string of the type: 4-12-31 into an array [4, 12, 31]
     */
    int *parse_C(string Cstr){
        int *C = (int *)malloc(sizeof(int) * (Cstr.size() + 1) / 2 );
        int i = 0, j = 0, k = 0;

        while((j = (int)Cstr.find('-', i)) > 0){
            C[k++] = atoi(Cstr.substr(i, j - i).c_str());
            i = ++j;
        }
        C[k] = atoi(Cstr.substr(i, Cstr.size() - i).c_str());
        return C;
    }

    int *parse_alphabet(const string A){
        int *char2int = (int *)malloc(sizeof(int) * 128);
        for(int i=0; i<A.size(); i++)
            char2int[A[i]] = i;
        return char2int;
    }

public:
    Bwt(string bwt, const int *C, const int *c2i){
        this->bwt = bwt;
        this->C = C;
        char2int = c2i;
    }
    Bwt(const string bwt, const string Cstr, const string A){
        this->bwt = bwt;
        C = parse_C(Cstr);
        char2int = parse_alphabet(A);
    }

    const int *C;
    const int *char2int;

    int rank(int i, char c){
        int cnt = 0;
        for(char cc: bwt.substr(0, i)){
            if(cc == c)
                cnt++;
        }
        return cnt;
    }
};


class Intervall{
private:
    Bwt *bwt;

public:
    int lb, ub;

    Intervall(Bwt *bwt_, char c){
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
        cout << "{";
        if(is_empty())
            cout << "-, -";
        else
            cout << lb << ", " << ub;
        cout << "}" << endl;
    }
};



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

uint MS(bit_vector ms, uint k){
    if(k == -1)
        return (uint) 1;
    return ((uint) bit_vector::select_1_type (&ms)(k + 1)) - (2 * k);
}

int find_k_prim(bit_vector runs, int k){
    size_t total_zeros = rank_support_v<0>(&runs)(runs.size());
    size_t zeros = rank_support_v<0>(&runs)(k+1);
    if(total_zeros == zeros)
        return (int) (runs.size());
    bit_vector::select_0_type b_sel(&runs);
    return (int) b_sel(zeros + 1);
}

bit_vector phase1(cst_sct3<> *st_of_s_ptr, string t, Bwt *bwt, bool verbose){
    cst_sct3<> st_of_s = *st_of_s_ptr;
    bit_vector runs(t.size());
    uint8_t k = t.size(), c = t[k - 1];
    Intervall II{bwt, static_cast<char>(c)};
    cst_sct3<>::node_type v = st_of_s.child(st_of_s.root(), c); // stree node

    while(--k > 0){
        c = t[k-1];
        II.bstep(c);
        if(II.is_empty()){
            runs[k] = 0;
            // update I to the parent of the proper locus of w until we can extend by 'c'
            do{
                v = st_of_s.parent(v);
                II.set((int) v.i, (int) v.j);
                II.bstep(c);
            } while(II.is_empty());
            v = st_of_s.wl(v, c); // update v
        } else {
            v = st_of_s.wl(v, c); // update v
            runs[k] = 1;
        }
        output_partial_vec(runs, k, "runs", verbose);
    }
    return runs;
}

void phase2(bit_vector runs, bit_vector *ms, cst_sct3<> *st_of_s_ptr, string t, Bwt *bwt, bool verbose){
    cst_sct3<> st_of_s = (*st_of_s_ptr);
    uint8_t k = 0, c = t[k], h_star = k + 1, k_prim, ms_idx = 0;
    Intervall II{bwt, static_cast<char>(c)};
    cst_sct3<>::node_type v = st_of_s.child(st_of_s.root(), c); // stree node


    while(k < t.size()){
        output_partial_vec(*ms, ms_idx, "ms", verbose);
        for(; !II.is_empty() && h_star < t.size(); ){
            c = t[h_star];
            II.bstep(c);
            if(!II.is_empty()){
                v = st_of_s.wl(v, c);
                h_star ++;
            }
        }
        for(int i = 0; i < h_star - k - MS(*ms, k - 1) + 1; i++)
            (*ms)[ms_idx++] = 0;
        if(h_star - k - MS(*ms, k - 1) + 1 > 0)
            (*ms)[ms_idx++] = 1;

        if(h_star < t.size()){
            do {
                v = st_of_s.parent(v);
                II.set((int)v.i, (int)v.j);
                II.bstep(t[h_star]);
            } while(II.is_empty());
            h_star = h_star + 1;
        }
        // k_prim: index of the first zero to the right of k in runs
        k_prim = find_k_prim(runs, k);

        for(int i=k+1; i<=k_prim - 1; i++)
            (*ms)[ms_idx++] = 1;

        // update v
        v = st_of_s.wl(v, c);
        k = k_prim;
    }
    output_partial_vec(*ms, ms_idx, "ms", verbose);
}


bit_vector compute_ms(string T, string S, string BWTfw, string BWTrev, string A, string Cstr, bool verbose){
    string Srev {S};
    for(int i=0; i<S.length(); i++)
        Srev[S.length() - i - 1] = S[i];

    bit_vector ms(2 * T.length());
    cst_sct3<> st_of_s, st_of_s_rev;
    construct_im(st_of_s, S, 1);
    construct_im(st_of_s_rev, Srev, 1);

    if(verbose){
        cout << " i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << endl;
        csXprintf(cout, "%2I %2S %3s %3P %2p %3B   %:3T", st_of_s.csa);
        cout << " i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << endl;
        csXprintf(cout, "%2I %2S %3s %3P %2p %3B   %:3T", st_of_s_rev.csa);
    }

    Bwt bwt(BWTfw, Cstr, A);
    bit_vector runs = phase1(&st_of_s, T, &bwt, verbose);
    output_partial_vec(runs, 0, "runs", verbose);
    if(verbose)
        cout << endl;

    Bwt bwt_rev(BWTrev, Cstr, A);
    phase2(runs, &ms, &st_of_s_rev, T, &bwt_rev, verbose);
    if(verbose){
        for(int i=0; i<T.length(); i++)
            cout << "MS[" << i << "] = " << (int) MS(ms, i) << endl;
    }
    return ms;
}


int main(int argc, char **argv){
    bit_vector ms;
    string T, S, BWTfw, BWTrev, A, Cstr;

    if(argc == 2){ // process file
        std::ifstream in_file {argv[1]};
        if(!in_file){
            cout << "could not open file " << argv[1] << endl;
            return 1;
        }

        while (in_file >> T >> S >> BWTfw >> BWTrev >> A >> Cstr){
            ms = compute_ms(T, S, BWTfw, BWTrev, A, Cstr, 0);
            cout << T + " " + S + " ";
            for (int i = 0; i < T.length(); i++)
                cout << MS(ms, i);
            cout << endl;
        }
    } else {
        //bbacbcacbbacbabcbbcc bcacaccaabacacbb bcabcccbac#cbaaaa bbccccaca#babcaaa #abc 0-1-7-11
        /*
        string T {"babcbacba"};
        string S {"bcbabaacacb"};
        string BWTfw {"bbbaccac#aab"};
        string BWTrev {"bcabcca#aabb"};
        string A {"#abc"};
        string Cstr {"0-1-7-11"};
        */

        std::fstream in_file {"/Users/denas/Desktop/FabioImplementation/software/indexed_ms/tests/a.txt"};
        in_file >> T >> S >> BWTfw >> BWTrev >> A >> Cstr;

        ms = compute_ms(T, S, BWTfw, BWTrev, A, Cstr, 1);

        cout << T + " " + S + " ";
        for (int i = 0; i < T.length(); i++)
            cout << MS(ms, i);
        cout << endl;

        //cout << "usage: " << argv[0] << "<file of strings>" << endl;
        //cout << "or" << endl;
        //cout << "usage: " << argv[0] << "<T> <S>" << endl;
        //return 1;
    }

    return 0;
}

