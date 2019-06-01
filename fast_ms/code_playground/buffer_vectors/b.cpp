#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <sdsl/vectors.hpp>

//#define V_SIZE 500000
#define B_SIZE (V_SIZE > 200 ? V_SIZE / 100 : V_SIZE)

int main(int argc, char **argv){
    {
        sdsl::int_vector_buffer<1>a("aa.bin", std::ios::out, B_SIZE);
        a[V_SIZE - 1] = 1;
    }
    sdsl::int_vector_buffer<1>a("aa.bin", std::ios::in, B_SIZE);
    return 0;
}


