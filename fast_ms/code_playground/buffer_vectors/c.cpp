#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <sdsl/vectors.hpp>

//#define V_SIZE 500000
#define B_SIZE (V_SIZE > 200 ? V_SIZE / 100 : V_SIZE)
using namespace std;

int main(int argc, char **argv)
{
    for(int i=1; i<3; i++){
        sdsl::int_vector_buffer<1>a("cc.bin." + std::to_string(i), std::ios::out, 32);
        for(int j=0; j<10; j++)
            a[j] = j % 2;
    }
    sdsl::int_vector_buffer<1>b("cc.bin", std::ios::out, 32);
    int k = 0;
    for(int i=1; i<3; i++){
        sdsl::int_vector_buffer<1>a("cc.bin." + std::to_string(i), std::ios::in, 32);
        for(int j=0; j<10; j++){
            b[k + j] = a[j];
            cout << "b[" << k + j << "] <- a["<< j << "] = " << a[j] << endl;
        }
        k += a.size();
    }
    return 0;
}



