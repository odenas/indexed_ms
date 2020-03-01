#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include"mt19937ar.h"
#include"naive_ms_range_max.h"
#include"ms_range_max.h"

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
	uint32_t first_select_res;
	uint32_t last_select_res;
	uint32_t ms_range_max_val1,ms_range_max_val2;
//	uint32_t ms_range_max_val3,ms_range_max_val4;
	uint32_t ms_range_max_pos1,ms_range_max_pos2;
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
		ms_range_max_val1=naive_ms_range_max(first_idx,last_idx,nbits,bitvec,&ms_range_max_pos1);
		ms_range_max_val2=ms_range_max_fast(first_idx,first_select_res,
			last_select_res,bitvec,&ms_range_max_pos2);
//		ms_range_max_val3=ms_at_idx(ms_range_max_pos1,nbits,bitvec);
//		ms_range_max_val4=ms_at_idx(ms_range_max_pos2,nbits,bitvec);
		if(ms_range_max_val1!=ms_range_max_val2 || ms_range_max_pos1!=ms_range_max_pos2)
			break;
	};
	if(i<nqueries)
	{
		printf("Failure at query number %d\n",i);
		printf("Start of range is %d\n",first_idx);
		printf("End of range is %d\n",last_idx);
		printf("naive ms range max val is %u\n",ms_range_max_val1);
		printf("naive ms range max pos is %u\n",ms_range_max_pos1);
		printf("byte-based ms range max val is %u\n",ms_range_max_val2);
		printf("byte-based ms range max pos is %u\n",ms_range_max_pos2);
	}
	else
		printf("Success\n");
	free(bitvec);
	return 0;
};
