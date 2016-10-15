/*
 * fabio_djamal_ms.cpp
 *
 *  Created on: Oct 13, 2016
 *      Author: denas
 */

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/csa_wt.hpp>

using namespace std;
using namespace sdsl;

int main(int argc, char **argv){

	csa_wt<> my_sa;
	construct_im(my_sa, "abracadabra", 1);

	for (size_t i=0; i<my_sa.size(); ++i){
		if(0 == my_sa.bwt[i]){
			cout << '#';
		} else {
			cout << my_sa.bwt[i];
		}
	}
	cout << endl;
		//cout << util::to_string(my_sa.bwt[i], 1) << endl;

	cout << " i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << endl;
	csXprintf(cout, "%2I %2S %3s %3P %2p %3B   %:3T", my_sa);
}

