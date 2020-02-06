#include"naive_ms_range_max.h"


uint32_t naive_ms_range_max(uint32_t first,uint32_t last,uint32_t bitvec_size,uint32_t * bitvec, uint32_t * _pos)
{
	uint32_t i;
	uint32_t max_ms=0;
	uint32_t curr_ms;
	uint32_t pos=first;
	for(i=first;i<=last;i++)
	{
		curr_ms=ms_at_idx(i,bitvec_size,bitvec);
		if(i==first || curr_ms>max_ms)
		{
			max_ms=curr_ms;
			pos=i;
		};
	};
	(*_pos)=pos;
	return max_ms;
};

uint64_t naive_ms_range_max64(uint64_t first,uint64_t last,uint64_t bitvec_size,uint64_t * bitvec, uint64_t *_pos)
{
	uint64_t i;
	uint64_t max_ms=0;
	uint64_t curr_ms;
	uint64_t pos=0;
	for(i=first;i<=last;i++)
	{
		curr_ms=ms_at_idx64(i,bitvec_size,bitvec);
		if(i==first || curr_ms>max_ms)
		{
			max_ms=curr_ms;
			pos=i;
		};
	};
	(*_pos)=pos;
	return max_ms;

};


int8_t naive_ms_range_max_byte_incr(uint8_t byte,int8_t *_pos)
{
	int8_t max_ms=-9;
	int8_t pos=-1;
	int8_t curr_ms=0;
	int8_t idx;
	int8_t next_bit;
	int8_t ms_pos=0;

	for(idx=0;idx<8;idx++)
	{
		next_bit=(byte>>idx)&1;
		if(next_bit==0)
			curr_ms++;
		else
		{
			curr_ms--;
			if(curr_ms>max_ms)
			{
				pos=ms_pos;
				max_ms=curr_ms;				
			};
			ms_pos++;
		}
	};
	(*_pos)=pos;
	if(max_ms==-9) max_ms=0;
	return max_ms;
};
