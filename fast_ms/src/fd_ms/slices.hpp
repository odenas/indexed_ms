//
//  slices.hpp
//  fast_ms
//
//  Created by denas on 10/27/17.
//  Copyright © 2017 denas. All rights reserved.
//

#ifndef slices_h
#define slices_h

#include <cassert>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>


using namespace std;


namespace fdms {

    template<typename size_type>
    class Slices {
        typedef pair<size_type, size_type> pair_t;

    public:
        size_t input_size, nslices;
        vector<pair_t> slices;

        Slices() : input_size{0}, nslices{0}, slices{vector<pair_t>(0)}
        {
        }

        /**
         * split 14 positions into 4 chunks
         * 14 = s * 4 + e (s = 3, e = 2)
         * so, allocate 4 chunks of size `s` and increase the size of the first `e` chunks by 1
         **/
        Slices(const size_type input_size, const size_type nslices) :
        input_size{input_size}, nslices{nslices}, slices{vector<pair_t>(nslices)}
        {
            size_type chunk = input_size / nslices;
            size_type extra = input_size % nslices;
            size_type step = 0;

            for (size_type i = 0, from = 0; i < nslices; i++) {
                // increase the size of the the first `extra` slices by 1
                step = chunk + (i < extra ? 1 : 0);
                assert(step > 0);
                slices[i] = std::make_pair(from, from + step);
                from += step;
            }
        }

        Slices(const Slices& other) :
        input_size{other.input_size}, nslices{other.nslices}, slices{vector<pair_t>(other.nslices)}
        {
            for (int i = 0; i < other.nslices; i++)
                slices[i] = other.slices[i];
        }
        
        Slices& operator=(const Slices& other){
            slices.resize(other.nslices);
            for (int i = 0; i < other.nslices; i++)
                slices[i] = other.slices[i];
            input_size = other.input_size;
            nslices = other.nslices;
            return *this;
        }

        pair_t operator[](size_type i) const {
            return slices[i];
        }
        
        size_type slice_length(size_type i) const {
            return slices[i].second - slices[i].first;
        }
        
        /*
         * Given an index in the domain of slices, infer the slice that contains it.
         * E.g., if a slice object has 4 slices over 14 positions
         * |....|....|....|..|
         * with slices [0, 4), [4, 8), [8, 12), [12, 14)
         * positions 0, 1, 2, 3 -> 0
         * positions 4, 5, 6, 7 -> 1
         * positions 8, 9, 10, 11 -> 2
         * positions 12, 13 -> 3
         */
        size_type slice_idx(const size_type i) const {
            if(i >= input_size)
            	throw "i > " + to_string(input_size);
            assert(i < input_size);
            size_type s_idx = 0; 
            while(slices[s_idx].second <= i)
                s_idx++;
            assert(slices[s_idx].first <= i && i < slices[s_idx].second);
            assert(s_idx < nslices);
            return s_idx;
        }

        string repr(size_type i) const {
            return string("[" + to_string(slices[i].first) + " .. " + to_string(slices[i].second) + ")");
        }
    };
};

#endif /* slices_h */
