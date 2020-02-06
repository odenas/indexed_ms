#include"naive_ms.h"

uint32_t naive_select(uint32_t bit_idx,uint32_t bitvec_size,const uint32_t * bitvec)
{
	uint32_t idx;
	uint32_t next_bit;
	for(idx=0;idx<bitvec_size;idx++)
	{
		next_bit=read_bit(idx,bitvec);
		if(next_bit)
			if((bit_idx--)==0)
				break;
	};
	return idx;


};
uint64_t naive_select64(uint64_t bit_idx,uint64_t bitvec_size,const uint64_t * bitvec)
{
	uint64_t idx;
	uint64_t next_bit;
	for(idx=0;idx<bitvec_size;idx++)
	{
		next_bit=read_bit64(idx,bitvec);
		if(next_bit)
			if((bit_idx--)==0)
				break;
	};
	return idx;


};

int8_t naive_count1_byte(uint8_t byte)
{
	int8_t count=0;
	uint8_t idx;

	for(idx=0;idx<8;idx++)
		count+=(byte>>idx)&1;
	return count;
};

int8_t naive_ms_incr_byte(uint8_t byte)
{
	int8_t res=0;
	uint8_t idx;
	uint8_t next_bit;
	for(idx=0;idx<8;idx++)
	{
		next_bit=(byte>>idx)&1;
		if(next_bit==0)
			res++;
		else
			res--;

	};
	return res;
};
