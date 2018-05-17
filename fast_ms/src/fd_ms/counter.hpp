/**
 * helper data structure for keeping counts of things
 */

#ifndef map_structs_h
#define map_structs_h

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

	    Counter(){}

	    void register_now(const string key, const std::chrono::time_point<std::chrono::_V2::system_clock, std::chrono::nanoseconds> start){
		    auto end = timer::now();
		    reg[key] = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		    cerr << "DONE (" << reg[key] / 1000 << " seconds)" << endl;
	    }
    };
};

#endif