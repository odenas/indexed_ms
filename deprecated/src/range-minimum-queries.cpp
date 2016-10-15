/*
 * range-minimum-queries.cpp
 *
 *  Created on: Oct 13, 2016
 *      Author: denas
 */


#include <sdsl/rmq_support.hpp>
#include <iostream>
#include <cstdlib>


using namespace sdsl;
using namespace std;


int main(int argc, char** argv){
	uint64_t len = 50;
	if(argc > 1){
		len = stoull(argv[1]);
	}

	int_vector<> v(len, 7); // a vector of length 50 (by default)
	for(uint64_t i=1; i<len; ++i){
		uint64_t x = v[i-1];
		v[i] = (x*7-1) % 97;
	}


	cout << "size_in_bytes(v) = " << size_in_bytes(v) << endl;
	rmq_succinct_sct<> rmq(&v);
	cout << "size_in_bytes(rmq) / v.size() = " << size_in_bytes(rmq) / v.size() << endl;

	if(v.size() < 100){
		cout << "v = " << v << endl;
	}

	uint64_t left = 0, right = v.size() - 1;

	// no need for v to answer queries
	util::clear(v);
	uint64_t count = 0;

	while(left < right){
		uint64_t min_pos = rmq(left, right);
		if(++count < 10){
			cout << "minimum of v[" << left << "..." << right << "] at index ";
			cout << min_pos << endl;
		}
		if(min_pos - left > right - min_pos){
			right = min_pos - 1;
		} else {
			left = min_pos + 1;
		}
	}

	cout << endl;
	cout << "write_structure<JSON_FORMAT>:" << endl;
    cout << "----------------------------- " << endl;
    write_structure<JSON_FORMAT>(rmq, cout);
    cout << endl;
}

