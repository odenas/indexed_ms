#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <sdsl/vectors.hpp>

//#define V_SIZE 500000

int main(int argc, char **argv){
	sdsl::bit_vector a(V_SIZE);
	a[a.size() - 1] = 1;
    return 0;
}

