#ifndef slices_h
#define slices_h

#include <iostream>
#include <fstream>
#include <vector>
#include <string>


using namespace std;


namespace fdms {
    template<typename size_type>
    class Slices{
        typedef pair<size_type, size_type> Interval;

    public:
        size_t input_size, nslices;
        vector<Interval> slices;
        
        Slices() : input_size{0}, nslices{0}, slices{vector<Interval>(0)} {}

        Slices(const size_type input_size, const size_type nslices) :
        input_size{input_size}, nslices{nslices}, slices{vector<Interval>(nslices)}
        {
            size_type chunk = input_size / nslices;
            size_type extra = input_size % nslices;
            size_type step = 0;
            
            //std::vector<Interval> slices(nslices);
            for(size_type i=0, from = 0; i<nslices; i++){
                step = chunk + (i < extra ? 1 : 0);
                slices[i] = std::make_pair(from, from + step);
                from += step;
            }
        }

        Slices(const Slices& other) :
        input_size{other.input_size}, nslices{other.nslices}, slices{vector<Interval>(other.nslices)}
        {
            for(int i=0; i<other.nslices; i++)
                slices[i] = other.slices[i];
        }
        
        Interval operator[](size_type i) const { return slices[i]; }

        string repr(size_type i) {
            return string("[" + to_string(slices[i].first) + " .. " + to_string(slices[i].second) + ")");
        }
    };
};

#endif /* slices_h */
