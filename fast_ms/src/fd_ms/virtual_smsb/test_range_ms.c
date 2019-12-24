#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include"mt19937ar.h"
#include"naive_ms_range_sum.h"
#include"range_ms_sum.h"

static uint32_t range_ms_sum_fast_in_blocks(uint32_t first_ms, uint32_t first_bit_idx,
			uint32_t last_bit_idx,uint32_t * bitvec)
{
	uint32_t ms_sum=first_ms;
	uint32_t curr_ms_sum;
	uint32_t virtual_ms=first_ms;
	uint32_t first_word_idx;
	uint32_t last_word_idx;
	uint32_t i;
	if(last_bit_idx<=first_bit_idx)
		return virtual_ms;
	first_word_idx=(first_bit_idx+1)/32;
	first_bit_idx=(first_bit_idx+1)%32;
	last_word_idx=last_bit_idx/32;
	last_bit_idx=last_bit_idx%32;
	if(last_word_idx==first_word_idx)
	{
		range_ms_sum_fast_ext(virtual_ms,first_bit_idx,last_bit_idx,
			&bitvec[first_word_idx],&virtual_ms,&curr_ms_sum);
		ms_sum+=curr_ms_sum;
		return ms_sum;		
	}
	range_ms_sum_fast_ext(virtual_ms,first_bit_idx,31,
		&bitvec[first_word_idx],&virtual_ms,&curr_ms_sum);
	ms_sum+=curr_ms_sum;
	for(i=first_word_idx+1;i<last_word_idx;i++)
	{
		range_ms_sum_fast_ext_word(virtual_ms,bitvec[i],&virtual_ms,&curr_ms_sum);
		ms_sum+=curr_ms_sum;
	};
	range_ms_sum_fast_ext(virtual_ms,0,last_bit_idx,
		&bitvec[last_word_idx],&virtual_ms,&curr_ms_sum);
	ms_sum+=curr_ms_sum;
	return ms_sum;
};

#define nqueries 4096
#define nbits 2048
#define max_query_range 256
#define max_first 256
//#define max_first 64

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
	uint32_t range_sum1,range_sum2,range_sum3,range_sum4;
	uint32_t virtual_ms;
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
		range_sum1=range_ms_sum_fast(first_ms,first_select_res,
			last_select_res,bitvec);
		range_sum2=naive_range_ms(first_idx,last_idx,nbits,bitvec);
		range_ms_sum_fast_ext(first_ms,first_select_res+1,
			last_select_res,bitvec,&virtual_ms,&range_sum3);
		range_sum3+=first_ms;
		range_sum4=range_ms_sum_fast_in_blocks(first_ms,first_select_res,
			last_select_res,bitvec);
		if(range_sum1!=range_sum2 || range_sum2!=range_sum3 || range_sum3!=range_sum4)
			break;
	};
	if(i<nqueries)
	{
		printf("Failure at query number %d\n",i);
		printf("Start of range is %d\n",first_idx);
		printf("End of range is %d\n",last_idx);
		printf("First bit is %d\n",first_select_res);
		printf("Last bit is %d\n",last_select_res);
		printf("naive ms range sum is %u\n",range_sum1);
		printf("byte-based ms range sum is %u\n",range_sum2);
		printf("virtual based ms range sum is %u\n",range_sum3);
		printf("virtual block-based ms range sum is %u\n",range_sum4);
	}
	else
		printf("Success\n");
	free(bitvec);
	return 0;
};
