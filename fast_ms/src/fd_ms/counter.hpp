/**
 * helper data structure for keeping counts of things
 */

#ifndef counter_h
#define counter_h

#include <map>
#include <string>
#include <chrono>


using namespace std;
using timer = std::chrono::high_resolution_clock;

namespace fdms{
    template<typename size_type>
    class Counter{
    public:
        map<string, size_type> reg;

        void register_now(const string& key,
            const std::chrono::time_point<std::chrono::_V2::system_clock, std::chrono::nanoseconds> start,
            const bool print_time=false)
        {
            auto end = timer::now();
            reg[key] = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            if(print_time)
                cerr << "DONE (" << reg[key] / 1000000 << " seconds)" << endl;
        }

        void register_add(const string& key,
            const std::chrono::time_point<std::chrono::_V2::system_clock, std::chrono::nanoseconds> start,
            const bool print_time=false)
        {
            auto end = timer::now();
            reg[key] += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            if(print_time)
                cerr << "DONE (" << reg[key] / 1000000 << " seconds)" << endl;
        }

    };
};
#endif
