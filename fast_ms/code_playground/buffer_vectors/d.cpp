#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <sdsl/vectors.hpp>

using namespace std;

int main(int argc, char **argv)
{
    {
        sdsl::int_vector_buffer<1>a("dd.bin", std::ios::out, 32);
        for(int j=0; j<10; j++)
            a[j] = j % 2;
    }
    sdsl::int_vector_buffer<1>b("dd.bin", std::ios::out, 32);
    for(int j=0; j<10; j++)
        cout << "d[" << j << "] = "<< b[j] << endl;
    return 0;
}




