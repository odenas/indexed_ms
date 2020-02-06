#include"naive_ms_range_sum.h"


uint32_t naive_range_ms(uint32_t first,uint32_t last,uint32_t bitvec_size,uint32_t * bitvec)
{
	uint32_t i;
	uint32_t sum=0;
	for(i=first;i<=last;i++)
		sum+=ms_at_idx(i,bitvec_size,bitvec);
	return sum;
};
uint64_t naive_range_ms64(uint64_t first,uint64_t last,uint64_t bitvec_size,const uint64_t * bitvec)
{
	uint64_t i;
	uint64_t sum=0;
	for(i=first;i<=last;i++)
		sum+=ms_at_idx64(i,bitvec_size,bitvec);
	return sum;
};


int32_t naive_range_ms_byte(int32_t virt_last_ms,uint8_t byte)
{
	int32_t sum=0;
	int32_t curr_ms=virt_last_ms;
	uint8_t idx;
	uint8_t next_bit;
	for(idx=0;idx<8;idx++)
	{
		next_bit=(byte>>idx)&1;
		if(next_bit==0)
			curr_ms++;
		else
		{
			curr_ms--;
			sum+=curr_ms;
		}
	};
	return sum;
};
