#include"ms_range_sum.h"

uint32_t range_ms_sum_fast(uint32_t first_ms, uint32_t first_bit_idx,
			uint32_t last_bit_idx,uint32_t * bitvec)
// Compute range_ms given the first ms, the bit that corresponds to that first ms (obtained with select)
// and the bit corresponding to the last ms (again obtained using select)
{
	uint32_t word_idx;
	uint32_t last_word_idx;
	uint32_t byte_idx;
	uint32_t last_byte_idx;
	uint32_t curr_word;
	uint8_t curr_byte;
	uint32_t range_sum;
	range_ms_params_t range_ms_params;

	word_idx=first_bit_idx/32;
	byte_idx=(first_bit_idx%32)/8;
	first_bit_idx=first_bit_idx%8;
	last_word_idx=last_bit_idx/32;
	last_byte_idx=(last_bit_idx%32)/8;
	last_bit_idx=last_bit_idx%8;

	curr_word=bitvec[word_idx];
	curr_byte=(curr_word>>(8*byte_idx))&0xff;
	curr_byte&=(0xff<<first_bit_idx);

	if(word_idx==last_word_idx)
	{
		if(byte_idx==last_byte_idx)
		{
			if(first_bit_idx==last_bit_idx)
				return first_ms;
			curr_byte&=(0xff>>(7-last_bit_idx));
			range_ms_params=range_ms_lookup_table[curr_byte];
			return (first_ms-first_bit_idx+1)*range_ms_params.sum_mul+range_ms_params.sum_add;
		}
		first_ms-=first_bit_idx;
		first_ms++;
		range_ms_params=range_ms_lookup_table[curr_byte];
		range_sum=first_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
		first_ms+=range_ms_params.ms_incr;
		for(byte_idx++;byte_idx<last_byte_idx;byte_idx++)
		{
			curr_byte=(curr_word>>(8*byte_idx))&0xff;
			range_ms_params=range_ms_lookup_table[curr_byte];
			range_sum+=first_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
			first_ms+=range_ms_params.ms_incr;
		};
		curr_byte=(curr_word>>(8*byte_idx))&0xff;
		curr_byte&=(0xff>>(7-last_bit_idx));
		range_ms_params=range_ms_lookup_table[curr_byte];
		range_sum+=first_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
		return range_sum;
	} 
	first_ms-=first_bit_idx;
	first_ms++;
	range_ms_params=range_ms_lookup_table[curr_byte];
	range_sum=first_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
	first_ms+=range_ms_params.ms_incr;
	for(byte_idx++;byte_idx<4;byte_idx++)
	{
		curr_byte=(curr_word>>(8*byte_idx))&0xff;
		range_ms_params=range_ms_lookup_table[curr_byte];
		range_sum+=first_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
		first_ms+=range_ms_params.ms_incr;
	};
	for(word_idx++;word_idx<last_word_idx;word_idx++)
	{
		curr_word=bitvec[word_idx];	
		for(byte_idx=0;byte_idx<4;byte_idx++)
		{
			curr_byte=(curr_word>>(8*byte_idx))&0xff;
			range_ms_params=range_ms_lookup_table[curr_byte];
			range_sum+=first_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
			first_ms+=range_ms_params.ms_incr;
		};
	};
	curr_word=bitvec[word_idx];
	for(byte_idx=0;byte_idx<last_byte_idx;byte_idx++)
	{
		curr_byte=(curr_word>>(8*byte_idx))&0xff;
		range_ms_params=range_ms_lookup_table[curr_byte];
		range_sum+=first_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
		first_ms+=range_ms_params.ms_incr;
	};
	curr_byte=(curr_word>>(8*byte_idx));
	curr_byte&=(0xff>>(7-last_bit_idx));
	range_ms_params=range_ms_lookup_table[curr_byte];
	range_sum+=first_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
	return range_sum;
};
uint64_t range_ms_sum_fast64(uint64_t first_ms, uint64_t first_bit_idx,
			uint64_t last_bit_idx, const uint64_t * bitvec)
{
	uint64_t word_idx;
	uint64_t last_word_idx;
	uint32_t byte_idx;
	uint32_t last_byte_idx;
	uint64_t curr_word;
	uint8_t curr_byte;
	uint64_t range_sum;
	range_ms_params_t range_ms_params;

	word_idx=first_bit_idx/64;
	byte_idx=(first_bit_idx%64)/8;
	first_bit_idx=first_bit_idx%8;
	last_word_idx=last_bit_idx/64;
	last_byte_idx=(last_bit_idx%64)/8;
	last_bit_idx=last_bit_idx%8;

	curr_word=bitvec[word_idx];
	curr_byte=(curr_word>>(8*byte_idx))&0xff;
	curr_byte&=(0xff<<first_bit_idx);

	if(word_idx==last_word_idx)
	{
		if(byte_idx==last_byte_idx)
		{
			if(first_bit_idx==last_bit_idx)
				return first_ms;
			curr_byte&=(0xff>>(7-last_bit_idx));
			range_ms_params=range_ms_lookup_table[curr_byte];
			return (first_ms-first_bit_idx+1)*range_ms_params.sum_mul+range_ms_params.sum_add;
		}
		first_ms-=first_bit_idx;
		first_ms++;
		range_ms_params=range_ms_lookup_table[curr_byte];
		range_sum=first_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
		first_ms+=range_ms_params.ms_incr;
		for(byte_idx++;byte_idx<last_byte_idx;byte_idx++)
		{
			curr_byte=(curr_word>>(8*byte_idx))&0xff;
			range_ms_params=range_ms_lookup_table[curr_byte];
			range_sum+=first_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
			first_ms+=range_ms_params.ms_incr;
		};
		curr_byte=(curr_word>>(8*byte_idx))&0xff;
		curr_byte&=(0xff>>(7-last_bit_idx));
		range_ms_params=range_ms_lookup_table[curr_byte];
		range_sum+=first_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
		return range_sum;
	} 
	first_ms-=first_bit_idx;
	first_ms++;
	range_ms_params=range_ms_lookup_table[curr_byte];
	range_sum=first_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
	first_ms+=range_ms_params.ms_incr;
	for(byte_idx++;byte_idx<8;byte_idx++)
	{
		curr_byte=(curr_word>>(8*byte_idx))&0xff;
		range_ms_params=range_ms_lookup_table[curr_byte];
		range_sum+=first_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
		first_ms+=range_ms_params.ms_incr;
	};
	for(word_idx++;word_idx<last_word_idx;word_idx++)
	{
		curr_word=bitvec[word_idx];	
		for(byte_idx=0;byte_idx<8;byte_idx++)
		{
			curr_byte=(curr_word>>(8*byte_idx))&0xff;
			range_ms_params=range_ms_lookup_table[curr_byte];
			range_sum+=first_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
			first_ms+=range_ms_params.ms_incr;
		};
	};
	curr_word=bitvec[word_idx];
	for(byte_idx=0;byte_idx<last_byte_idx;byte_idx++)
	{
		curr_byte=(curr_word>>(8*byte_idx))&0xff;
		range_ms_params=range_ms_lookup_table[curr_byte];
		range_sum+=first_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
		first_ms+=range_ms_params.ms_incr;
	};
	curr_byte=(curr_word>>(8*byte_idx));
	curr_byte&=(0xff>>(7-last_bit_idx));
	range_ms_params=range_ms_lookup_table[curr_byte];
	range_sum+=first_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
	return range_sum;

};

void range_ms_sum_fast_ext(uint32_t in_virtual_ms, uint32_t first_bit_idx,
			uint32_t last_bit_idx,uint32_t * bitvec,
			uint32_t * out_virtual_ms, uint32_t * ms_range_sum)
{

	uint32_t word_idx;
	uint32_t last_word_idx;
	uint32_t byte_idx;
	uint32_t last_byte_idx;
	uint32_t curr_word;
	uint8_t curr_byte;
	uint32_t range_sum;
	uint32_t virtual_ms=in_virtual_ms;
	range_ms_params_t range_ms_params;
	if(first_bit_idx>last_bit_idx)
	{
		range_sum=0;
		goto set_output;
	};
	word_idx=first_bit_idx/32;
	byte_idx=(first_bit_idx%32)/8;
	first_bit_idx=first_bit_idx%8;
	last_word_idx=last_bit_idx/32;
	last_byte_idx=(last_bit_idx%32)/8;
	last_bit_idx=last_bit_idx%8;

	curr_word=bitvec[word_idx];
	curr_byte=(curr_word>>(8*byte_idx))&0xff;
	curr_byte&=(0xff<<first_bit_idx);

	if(word_idx==last_word_idx)
	{
		if(byte_idx==last_byte_idx)
		{
			curr_byte&=(0xff>>(7-last_bit_idx));
			range_ms_params=range_ms_lookup_table[curr_byte];
			range_sum=(virtual_ms-first_bit_idx)*range_ms_params.sum_mul+range_ms_params.sum_add;
			virtual_ms+=range_ms_params.ms_incr-(7-last_bit_idx)-first_bit_idx;
			goto set_output;
		}
		virtual_ms-=first_bit_idx;
		range_ms_params=range_ms_lookup_table[curr_byte];
		range_sum=virtual_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
		virtual_ms+=range_ms_params.ms_incr;
		for(byte_idx++;byte_idx<last_byte_idx;byte_idx++)
		{
			curr_byte=(curr_word>>(8*byte_idx))&0xff;
			range_ms_params=range_ms_lookup_table[curr_byte];
			range_sum+=virtual_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
			virtual_ms+=range_ms_params.ms_incr;
		};
		curr_byte=(curr_word>>(8*byte_idx))&0xff;
		curr_byte&=(0xff>>(7-last_bit_idx));
		range_ms_params=range_ms_lookup_table[curr_byte];
		range_sum+=virtual_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
		virtual_ms+=range_ms_params.ms_incr-(7-last_bit_idx);
		goto set_output;
	} 
	virtual_ms-=first_bit_idx;
	range_ms_params=range_ms_lookup_table[curr_byte];
	range_sum=virtual_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
	virtual_ms+=range_ms_params.ms_incr;
	for(byte_idx++;byte_idx<4;byte_idx++)
	{
		curr_byte=(curr_word>>(8*byte_idx))&0xff;
		range_ms_params=range_ms_lookup_table[curr_byte];
		range_sum+=virtual_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
		virtual_ms+=range_ms_params.ms_incr;
	};
	for(word_idx++;word_idx<last_word_idx;word_idx++)
	{
		curr_word=bitvec[word_idx];	
		for(byte_idx=0;byte_idx<4;byte_idx++)
		{
			curr_byte=(curr_word>>(8*byte_idx))&0xff;
			range_ms_params=range_ms_lookup_table[curr_byte];
			range_sum+=virtual_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
			virtual_ms+=range_ms_params.ms_incr;
		};
	};
	curr_word=bitvec[word_idx];
	for(byte_idx=0;byte_idx<last_byte_idx;byte_idx++)
	{
		curr_byte=(curr_word>>(8*byte_idx))&0xff;
		range_ms_params=range_ms_lookup_table[curr_byte];
		range_sum+=virtual_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
		virtual_ms+=range_ms_params.ms_incr;
	};
	curr_byte=(curr_word>>(8*byte_idx));
	curr_byte&=(0xff>>(7-last_bit_idx));
	range_ms_params=range_ms_lookup_table[curr_byte];
	range_sum+=virtual_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
	virtual_ms+=range_ms_params.ms_incr-(7-last_bit_idx);
set_output:
	(*ms_range_sum)=range_sum;
	(*out_virtual_ms)=virtual_ms;

};
// Compute the contribution to sum of ms values and updates the virtual ms sum
void range_ms_sum_fast_ext64(uint64_t in_virtual_ms, uint64_t first_bit_idx,
			uint64_t last_bit_idx,uint64_t * bitvec,
			uint64_t * out_virtual_ms, uint64_t * ms_range_sum)
{

	uint64_t word_idx;
	uint64_t last_word_idx;
	uint64_t byte_idx;
	uint64_t last_byte_idx;
	uint64_t curr_word;
	uint8_t curr_byte;
	uint64_t range_sum;
	uint64_t virtual_ms=in_virtual_ms;
	range_ms_params_t range_ms_params;
	if(first_bit_idx>last_bit_idx)
	{
		range_sum=0;
		goto set_output;
	};
	word_idx=first_bit_idx/64;
	byte_idx=(first_bit_idx%64)/8;
	first_bit_idx=first_bit_idx%8;
	last_word_idx=last_bit_idx/64;
	last_byte_idx=(last_bit_idx%64)/8;
	last_bit_idx=last_bit_idx%8;

	curr_word=bitvec[word_idx];
	curr_byte=(curr_word>>(8*byte_idx))&0xff;
	curr_byte&=(0xff<<first_bit_idx);

	if(word_idx==last_word_idx)
	{
		if(byte_idx==last_byte_idx)
		{
			curr_byte&=(0xff>>(7-last_bit_idx));
			range_ms_params=range_ms_lookup_table[curr_byte];
			range_sum=(virtual_ms-first_bit_idx)*range_ms_params.sum_mul+range_ms_params.sum_add;
			virtual_ms+=range_ms_params.ms_incr-(7-last_bit_idx)-first_bit_idx;
			goto set_output;
		}
		virtual_ms-=first_bit_idx;
		range_ms_params=range_ms_lookup_table[curr_byte];
		range_sum=virtual_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
		virtual_ms+=range_ms_params.ms_incr;
		for(byte_idx++;byte_idx<last_byte_idx;byte_idx++)
		{
			curr_byte=(curr_word>>(8*byte_idx))&0xff;
			range_ms_params=range_ms_lookup_table[curr_byte];
			range_sum+=virtual_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
			virtual_ms+=range_ms_params.ms_incr;
		};
		curr_byte=(curr_word>>(8*byte_idx))&0xff;
		curr_byte&=(0xff>>(7-last_bit_idx));
		range_ms_params=range_ms_lookup_table[curr_byte];
		range_sum+=virtual_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
		virtual_ms+=range_ms_params.ms_incr-(7-last_bit_idx);
		goto set_output;
	} 
	virtual_ms-=first_bit_idx;
	range_ms_params=range_ms_lookup_table[curr_byte];
	range_sum=virtual_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
	virtual_ms+=range_ms_params.ms_incr;
	for(byte_idx++;byte_idx<8;byte_idx++)
	{
		curr_byte=(curr_word>>(8*byte_idx))&0xff;
		range_ms_params=range_ms_lookup_table[curr_byte];
		range_sum+=virtual_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
		virtual_ms+=range_ms_params.ms_incr;
	};
	for(word_idx++;word_idx<last_word_idx;word_idx++)
	{
		curr_word=bitvec[word_idx];	
		for(byte_idx=0;byte_idx<8;byte_idx++)
		{
			curr_byte=(curr_word>>(8*byte_idx))&0xff;
			range_ms_params=range_ms_lookup_table[curr_byte];
			range_sum+=virtual_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
			virtual_ms+=range_ms_params.ms_incr;
		};
	};
	curr_word=bitvec[word_idx];
	for(byte_idx=0;byte_idx<last_byte_idx;byte_idx++)
	{
		curr_byte=(curr_word>>(8*byte_idx))&0xff;
		range_ms_params=range_ms_lookup_table[curr_byte];
		range_sum+=virtual_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
		virtual_ms+=range_ms_params.ms_incr;
	};
	curr_byte=(curr_word>>(8*byte_idx));
	curr_byte&=(0xff>>(7-last_bit_idx));
	range_ms_params=range_ms_lookup_table[curr_byte];
	range_sum+=virtual_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
	virtual_ms+=range_ms_params.ms_incr-(7-last_bit_idx);
set_output:
	(*ms_range_sum)=range_sum;
	(*out_virtual_ms)=virtual_ms;

};

