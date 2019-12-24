#ifndef mt19937ar_h
#define mt19937ar_h
void init_genrand(unsigned long long s);
unsigned long long genrand_int32(void);

void init_by_array(unsigned long long init_key[], int key_length);
long genrand_int31(void);
double genrand_real1(void);
double genrand_real2(void);
double genrand_real3(void);
double genrand_res53(void);

#endif
