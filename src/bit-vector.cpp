/*
 * bit-vector.cpp
 *
 *  Created on: Oct 10, 2016
 *      Author: denas
 */




#include <sdsl/bit_vectors.hpp>
#include <iostream>

using namespace std;
using namespace sdsl;

int main (){
	bit_vector b(10000000, 0);
	b[8] = 1;
	rank_support_v<> rb(&b);

	cout << rb(8) << endl;
	cout << rb(9) << endl;

	cout << "Size of b in MB: " << size_in_mega_bytes(b) << endl;
	cout << "Size of rb in MB: " << size_in_mega_bytes(rb) << endl;

	// a rrr-vector and the associated rank support vector
	rrr_vector<127> rrrb(b);
	rrr_vector<127>::rank_1_type rank_rrrb(&rrrb);

	cout << rank_rrrb(8) << endl;
	cout << rank_rrrb(9) << endl;

	cout << "Size of rrrb in MB: " << size_in_mega_bytes(rrrb) << endl;
	cout << "Size of rank_rrrb in MB: " << size_in_mega_bytes(rank_rrrb) << endl;

	rrr_vector<127>::select_1_type select_rrrb(&rrrb);
	cout << "position of the first 1 in b: " << select_rrrb(1) << endl;

	bit_vector x;
	util::assign(x, bit_vector(10000000, 1));

	int_vector<> v(100, 2, 7);

	cout << "v[5]=" << v[5] << endl;
	v[5] = 120;
	cout << "v[5]=" << v[5] << endl;

	int_vector<32> w(100, 4);
	write_structure<JSON_FORMAT>(rrrb, cout);
	cout << endl;
}
