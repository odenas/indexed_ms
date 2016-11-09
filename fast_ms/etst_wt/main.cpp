//
//  test_wt.cpp
//  fast_ms
//
//  Created by denas on 11/4/16.
//  Copyright © 2016 denas. All rights reserved.
//

#include <stdio.h>
#include <string.h>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/bit_vectors.hpp>
using namespace std;
using namespace sdsl;



void output_partial_vec(int_vector<> v, const int idx, const char name[]){
    //string name {"x"
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

int r_rank(wt_huff<> wtree, int i, char c){
    return (int) wtree.rank(i, c);
}

int r_rrank(string bwt, int i, char c){
    int cnt = 0;
    for(char cc: bwt.substr(0, i)){
        if(cc == c)
            cnt++;
    }
    return cnt;
}

int main_rank(int argc, char **argv){
    char xx[11] = {'a', 'b', 'b', 'a', 'b', 'a', 'a', 'b', 'b', 'a', '#'};
    char x[3] = {'#', 'a', 'b'};
    wt_huff<> wtree;

    string tmp_file = ram_file_name(util::to_string(util::pid()) + "_" + util::to_string(util::id()));
    //cout << tmp_file << endl;
    store_to_file(xx, tmp_file);
    construct(wtree, tmp_file, 1);
    ram_fs::remove(tmp_file);

    //construct_im(wtree, xx, 1);
    //construct(wtree, "/Users/denas/Desktop/FabioImplementation/software/indexed_ms/aa", 1);

    for(int i=0; i<wtree.size(); i++)
        cout << (char)wtree[i] << " ";
    cout << endl;
    int_vector<> y(wtree.size() + 1), yy(wtree.size() + 1);


    for(char s: x){
        cout << s << ": " << endl;
        for(int i=0; i<=wtree.size(); i++){
            y[i] = r_rank(wtree, i, s);
            yy[i] = r_rrank(xx, i, s);
        }

        output_partial_vec(y, 0, " y");
        output_partial_vec(yy, 0, "yy");
        
    }
    return 0;
}

int main(int argc, char **argv){
    string t {"111111000110011100011100110100111000011000111111010010100011100011100001111100110000011101011010001100100110001111010100110011110001110100101000111001100000111111000011100100111100110010100011100110001100011000111110101000101010000111110110011000000000001100000000000000000000000000000000000000000000"};

    bit_vector ms(t.size());
    for(int i=0; i<t.size(); i++){
        int j = stoi(t.substr(i, 1));
        //cout << i << ": " << t[i] << " " << j << endl;
        ms[i] = j;
    }
    cout << t  << endl;
    cout << ms << endl;

    for(int k=1; k<127; k++)
        cout << k << ":" << ((uint) bit_vector::select_1_type (&ms)(k + 1)) << endl;
    return 0;
}

