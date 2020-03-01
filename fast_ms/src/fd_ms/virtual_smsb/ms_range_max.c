#include"ms_range_max.h"

uint32_t ms_range_max_fast(uint32_t first_ms_idx, uint32_t first_bit_idx,
			uint32_t last_bit_idx,uint32_t * bitvec, uint32_t * _pos)
{
	uint32_t word_idx;
	uint32_t last_word_idx;
	uint32_t byte_idx;
	uint32_t last_byte_idx;
	uint32_t curr_word;
	uint8_t curr_byte;
	uint32_t first_ms=first_bit_idx-2*first_ms_idx;
	uint32_t max_ms_val=first_ms;
	uint32_t max_ms_pos=first_ms_idx;
	uint32_t curr_ms_pos;
	uint32_t curr_ms;
	ms_range_max_params_t ms_range_max_params;

	word_idx=first_bit_idx/32;
	byte_idx=(first_bit_idx%32)/8;
	first_bit_idx=first_bit_idx%8;
	last_word_idx=last_bit_idx/32;
	last_byte_idx=(last_bit_idx%32)/8;
	last_bit_idx=last_bit_idx%8;

	curr_word=bitvec[word_idx];
	curr_byte=(curr_word>>(8*byte_idx))&0xff;
	curr_byte&=(0xff<<first_bit_idx);
	curr_ms=first_ms-(first_bit_idx-1);
	curr_ms_pos=first_ms_idx;

	if(word_idx==last_word_idx)
	{
		if(byte_idx==last_byte_idx)
		{
			if(first_bit_idx==last_bit_idx)
			{
				(*_pos)=first_ms_idx;
				return first_ms;
			};
			curr_byte&=(0xff>>(7-last_bit_idx));
			ms_range_max_params=ms_range_max_lookup_table[curr_byte];
			max_ms_val=curr_ms+ms_range_max_params.max_incr;
			max_ms_pos=curr_ms_pos+ms_range_max_params.max_pos;
			(*_pos)=max_ms_pos;
			return max_ms_val;
		}
		ms_range_max_params=ms_range_max_lookup_table[curr_byte];
		max_ms_val=curr_ms+ms_range_max_params.max_incr;
		max_ms_pos=curr_ms_pos+ms_range_max_params.max_pos;
		curr_ms+=ms_range_max_params.ms_incr;
		curr_ms_pos+=ms_range_max_params.count1;
		for(byte_idx++;byte_idx<last_byte_idx;byte_idx++)
		{
			curr_byte=(curr_word>>(8*byte_idx))&0xff;
			ms_range_max_params=ms_range_max_lookup_table[curr_byte];
			if(curr_ms+ms_range_max_params.max_incr>max_ms_val)
			{
				max_ms_val=curr_ms+ms_range_max_params.max_incr;
				max_ms_pos=curr_ms_pos+ms_range_max_params.max_pos;
			};
			curr_ms+=ms_range_max_params.ms_incr;
			curr_ms_pos+=ms_range_max_params.count1;
		};
		curr_byte=(curr_word>>(8*byte_idx))&0xff;
		curr_byte&=(0xff>>(7-last_bit_idx));
		ms_range_max_params=ms_range_max_lookup_table[curr_byte];
		if(curr_ms+ms_range_max_params.max_incr>max_ms_val)
		{
			max_ms_val=curr_ms+ms_range_max_params.max_incr;
			max_ms_pos=curr_ms_pos+ms_range_max_params.max_pos;
		};
		(*_pos)=max_ms_pos;
		return max_ms_val;
	} 
	ms_range_max_params=ms_range_max_lookup_table[curr_byte];
	max_ms_val=curr_ms+ms_range_max_params.max_incr;
	max_ms_pos=curr_ms_pos+ms_range_max_params.max_pos;
	curr_ms+=ms_range_max_params.ms_incr;
	curr_ms_pos+=ms_range_max_params.count1;
	for(byte_idx++;byte_idx<4;byte_idx++)
	{
		curr_byte=(curr_word>>(8*byte_idx))&0xff;
		ms_range_max_params=ms_range_max_lookup_table[curr_byte];
		if(curr_ms+ms_range_max_params.max_incr>max_ms_val)
		{
			max_ms_val=curr_ms+ms_range_max_params.max_incr;
			max_ms_pos=curr_ms_pos+ms_range_max_params.max_pos;
		};
		curr_ms+=ms_range_max_params.ms_incr;
		curr_ms_pos+=ms_range_max_params.count1;
	};
	for(word_idx++;word_idx<last_word_idx;word_idx++)
	{
		curr_word=bitvec[word_idx];	
		for(byte_idx=0;byte_idx<4;byte_idx++)
		{
			curr_byte=(curr_word>>(8*byte_idx))&0xff;
			ms_range_max_params=ms_range_max_lookup_table[curr_byte];
			if(curr_ms+ms_range_max_params.max_incr>max_ms_val)
			{
				max_ms_val=curr_ms+ms_range_max_params.max_incr;
				max_ms_pos=curr_ms_pos+ms_range_max_params.max_pos;
			};
			curr_ms+=ms_range_max_params.ms_incr;
			curr_ms_pos+=ms_range_max_params.count1;
		};
	};
	curr_word=bitvec[word_idx];
	for(byte_idx=0;byte_idx<last_byte_idx;byte_idx++)
	{
		curr_byte=(curr_word>>(8*byte_idx))&0xff;
		ms_range_max_params=ms_range_max_lookup_table[curr_byte];
		if(curr_ms+ms_range_max_params.max_incr>max_ms_val)
		{
			max_ms_val=curr_ms+ms_range_max_params.max_incr;
			max_ms_pos=curr_ms_pos+ms_range_max_params.max_pos;
		};
		curr_ms+=ms_range_max_params.ms_incr;
		curr_ms_pos+=ms_range_max_params.count1;
	};
	curr_byte=(curr_word>>(8*byte_idx));
	curr_byte&=(0xff>>(7-last_bit_idx));
	ms_range_max_params=ms_range_max_lookup_table[curr_byte];
	if(curr_ms+ms_range_max_params.max_incr>max_ms_val)
	{
		max_ms_val=curr_ms+ms_range_max_params.max_incr;
		max_ms_pos=curr_ms_pos+ms_range_max_params.max_pos;
	};
	(*_pos)=max_ms_pos;
	return max_ms_val;
};

uint64_t ms_range_max_fast64(uint64_t first_ms_idx, uint64_t first_bit_idx,
			uint64_t last_bit_idx, const uint64_t * bitvec, uint64_t * _pos)
{
	uint64_t word_idx;
	uint64_t last_word_idx;
	uint64_t byte_idx;
	uint64_t last_byte_idx;
	uint64_t curr_word;
	uint8_t curr_byte;
	uint64_t first_ms=first_bit_idx-2*first_ms_idx;
	uint64_t max_ms_val=first_ms;
	uint64_t max_ms_pos=first_ms_idx;
	uint64_t curr_ms_pos;
	uint64_t curr_ms;
	ms_range_max_params_t ms_range_max_params;
	
	word_idx=first_bit_idx/64;
	byte_idx=(first_bit_idx%64)/8;
	first_bit_idx=first_bit_idx%8;
	last_word_idx=last_bit_idx/64;
	last_byte_idx=(last_bit_idx%64)/8;
	last_bit_idx=last_bit_idx%8;

	curr_word=bitvec[word_idx];
	curr_byte=(curr_word>>(8*byte_idx))&0xff;
	curr_byte&=(0xff<<first_bit_idx);
	curr_ms=first_ms-(first_bit_idx-1);
	curr_ms_pos=first_ms_idx;

	if(word_idx==last_word_idx)
	{
		if(byte_idx==last_byte_idx)
		{
			if(first_bit_idx==last_bit_idx)
			{
				(*_pos)=first_ms_idx;
				return first_ms;
			};
			curr_byte&=(0xff>>(7-last_bit_idx));
			ms_range_max_params=ms_range_max_lookup_table[curr_byte];
			max_ms_val=curr_ms+ms_range_max_params.max_incr;
			max_ms_pos=curr_ms_pos+ms_range_max_params.max_pos;
			(*_pos)=max_ms_pos;
			return max_ms_val;
		}
		ms_range_max_params=ms_range_max_lookup_table[curr_byte];
		max_ms_val=curr_ms+ms_range_max_params.max_incr;
		max_ms_pos=curr_ms_pos+ms_range_max_params.max_pos;
		curr_ms+=ms_range_max_params.ms_incr;
		curr_ms_pos+=ms_range_max_params.count1;
		for(byte_idx++;byte_idx<last_byte_idx;byte_idx++)
		{
			curr_byte=(curr_word>>(8*byte_idx))&0xff;
			ms_range_max_params=ms_range_max_lookup_table[curr_byte];
			if(curr_ms+ms_range_max_params.max_incr>max_ms_val)
			{
				max_ms_val=curr_ms+ms_range_max_params.max_incr;
				max_ms_pos=curr_ms_pos+ms_range_max_params.max_pos;
			};
			curr_ms+=ms_range_max_params.ms_incr;
			curr_ms_pos+=ms_range_max_params.count1;
		};
		curr_byte=(curr_word>>(8*byte_idx))&0xff;
		curr_byte&=(0xff>>(7-last_bit_idx));
		ms_range_max_params=ms_range_max_lookup_table[curr_byte];
		if(curr_ms+ms_range_max_params.max_incr>max_ms_val)
		{
			max_ms_val=curr_ms+ms_range_max_params.max_incr;
			max_ms_pos=curr_ms_pos+ms_range_max_params.max_pos;
		};
		(*_pos)=max_ms_pos;
		return max_ms_val;
	} 
	ms_range_max_params=ms_range_max_lookup_table[curr_byte];
	max_ms_val=curr_ms+ms_range_max_params.max_incr;
	max_ms_pos=curr_ms_pos+ms_range_max_params.max_pos;
	curr_ms+=ms_range_max_params.ms_incr;
	curr_ms_pos+=ms_range_max_params.count1;
	for(byte_idx++;byte_idx<8;byte_idx++)
	{
		curr_byte=(curr_word>>(8*byte_idx))&0xff;
		ms_range_max_params=ms_range_max_lookup_table[curr_byte];
		if(curr_ms+ms_range_max_params.max_incr>max_ms_val)
		{
			max_ms_val=curr_ms+ms_range_max_params.max_incr;
			max_ms_pos=curr_ms_pos+ms_range_max_params.max_pos;
		};
		curr_ms+=ms_range_max_params.ms_incr;
		curr_ms_pos+=ms_range_max_params.count1;
	};
	for(word_idx++;word_idx<last_word_idx;word_idx++)
	{
		curr_word=bitvec[word_idx];	
		for(byte_idx=0;byte_idx<8;byte_idx++)
		{
			curr_byte=(curr_word>>(8*byte_idx))&0xff;
			ms_range_max_params=ms_range_max_lookup_table[curr_byte];
			if(curr_ms+ms_range_max_params.max_incr>max_ms_val)
			{
				max_ms_val=curr_ms+ms_range_max_params.max_incr;
				max_ms_pos=curr_ms_pos+ms_range_max_params.max_pos;
			};
			curr_ms+=ms_range_max_params.ms_incr;
			curr_ms_pos+=ms_range_max_params.count1;
		};
	};
	curr_word=bitvec[word_idx];
	for(byte_idx=0;byte_idx<last_byte_idx;byte_idx++)
	{
		curr_byte=(curr_word>>(8*byte_idx))&0xff;
		ms_range_max_params=ms_range_max_lookup_table[curr_byte];
		if(curr_ms+ms_range_max_params.max_incr>max_ms_val)
		{
			max_ms_val=curr_ms+ms_range_max_params.max_incr;
			max_ms_pos=curr_ms_pos+ms_range_max_params.max_pos;
		};
		curr_ms+=ms_range_max_params.ms_incr;
		curr_ms_pos+=ms_range_max_params.count1;
	};
	curr_byte=(curr_word>>(8*byte_idx));
	curr_byte&=(0xff>>(7-last_bit_idx));
	ms_range_max_params=ms_range_max_lookup_table[curr_byte];
	if(curr_ms+ms_range_max_params.max_incr>max_ms_val)
	{
		max_ms_val=curr_ms+ms_range_max_params.max_incr;
		max_ms_pos=curr_ms_pos+ms_range_max_params.max_pos;
	};
	(*_pos)=max_ms_pos;
	return max_ms_val;
};

