#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include"mt19937ar.h"
#include"naive_ms_range_sum.h"
#include"range_ms_sum.h"

#define nqueries 4096
#define nbits 2048
#define max_query_range 256
#define max_first 256

int main()
{
	uint64_t * bitvec;
	uint64_t nwords=(nbits+63)/64;
	uint64_t i;
	uint64_t first_idx;
	uint64_t last_idx;
	uint64_t first_ms;
	uint64_t first_select_res;
	uint64_t last_select_res;
	uint64_t range_sum1,range_sum2;
	bitvec=(uint64_t *) malloc(nwords*8);
	
	init_genrand(time(0));
	for(i=0;i<nwords/4;i++)
		bitvec[i]=0;
	for(;i<nwords;i++)
		bitvec[i]=((uint64_t)(genrand_int32())<<32)+genrand_int32();
	for(i=0;i<nqueries;i++)
	{
		first_idx=genrand_int32()%max_first;
		last_idx=first_idx+genrand_int32()%max_query_range;
		first_select_res=naive_select64(first_idx,nbits,bitvec);
		last_select_res=naive_select64(last_idx,nbits,bitvec);
		first_ms=first_select_res-2*first_idx;
		range_sum1=range_ms_sum_fast64(first_ms,first_select_res,
			last_select_res,bitvec);
		range_sum2=naive_range_ms64(first_idx,last_idx,nbits,bitvec);
		if(range_sum1!=range_sum2)
			break;
	};
	if(i<nqueries)
	{
		printf("Failure at query number %d\n",(int)i);
		printf("Start of range is %d\n",(int)first_idx);
		printf("End of range is %d\n",(int)last_idx);
		printf("naive ms range sum is %u\n",(int)range_sum1);
		printf("byte-based ms range sum is %u\n",(int)range_sum2);
	}
	else
		printf("Success\n");
	free(bitvec);
	return 0;
};
