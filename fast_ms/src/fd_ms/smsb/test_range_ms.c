#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include"mt19937ar.h"
#include"naive_ms_range_sum.h"
#include"range_ms_sum.h"

//#define nqueries 1024
//#define nbits 65536
//#define max_query_range 4096
//#define max_first 16384

#define nqueries 4096
#define nbits 2048
#define max_query_range 256
#define max_first 256

int main()
{
	uint32_t * bitvec;
	uint32_t nwords=(nbits+31)/32;
	uint32_t i;
	uint32_t first_idx;
	uint32_t last_idx;
	uint32_t first_ms;
	uint32_t first_select_res;
	uint32_t last_select_res;
	uint32_t range_sum1,range_sum2;
	bitvec=(uint32_t *) malloc(nwords*4);
	
	init_genrand(time(0));
	for(i=0;i<nwords/4;i++)
		bitvec[i]=0;
	for(;i<nwords;i++)
		bitvec[i]=genrand_int32();
	for(i=0;i<nqueries;i++)
	{
		first_idx=genrand_int32()%max_first;
		last_idx=first_idx+genrand_int32()%max_query_range;
		first_select_res=naive_select(first_idx,nbits,bitvec);
		last_select_res=naive_select(last_idx,nbits,bitvec);
		first_ms=first_select_res-2*first_idx;
		/**
		 * cmpute the sum of ms values in range [i, j]
		 * first_ms: ms[i]
		 * first_select_res: select1(i)
		 * last_select_res: select1(j)
		 * bitvec: ms
		 */
		range_sum1=range_ms_sum_fast(first_ms,first_select_res,
			last_select_res,bitvec);
		range_sum2=naive_range_ms(first_idx,last_idx,nbits,bitvec);
		if(range_sum1!=range_sum2)
			break;
	};
	if(i<nqueries)
	{
		printf("Failure at query number %d\n",i);
		printf("Start of range is %d\n",first_idx);
		printf("End of range is %d\n",last_idx);
		printf("naive ms range sum is %u\n",range_sum1);
		printf("byte-based ms range sum is %u\n",range_sum2);
	}
	else
		printf("Success\n");
	free(bitvec);
	return 0;
};
